#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "Netlist.h"
#include "Solver.h"
#include "Waveform.h"

#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <cmath>
#include <complex>
#include <unordered_map>

#define PYBIND11_LIMITED_API


namespace py = pybind11;

// ==============================================================================
//  ESTRUTURAS AUXILIARES
// ==============================================================================
struct FastPWMController {
    Netlist::Component* sw;
    double Ts;
    double on_time;
    double phase_offset;
    double dead_time;
};

static std::string upper(std::string s) {
    for (auto& c : s) c = (char)std::toupper((unsigned char)c);
    return s;
}

// ==============================================================================
//  CLASSE WRAPPER (PyCircuit)
// ==============================================================================
class PyCircuit {
public:
    Netlist nl;
    std::map<std::string, double> pwm_config;

    // Persistência de estado (Warm Start)
    Vector last_state;
    bool has_state = false;

    // Parâmetros do Solver
    int solver_max_iter = 100;
    double solver_tol = 1e-4;
    double solver_eps = 1e-6; 

    // Construtor
    PyCircuit() {}

    // Configuração dos parâmetros do solver
    void set_params(int max_iter, double tol, double eps = 1e-6) {
        this->solver_max_iter = max_iter;
        this->solver_tol = tol;
        this->solver_eps = eps;
    }

    // Reseta o estado (para começar simulação do zero)
    void reset_state() {
        has_state = false;
        last_state = Vector();
    }

    // Atualiza valor de componente existente
    void update_value(std::string name, double val) {
        auto* c = nl.find_component(name);
        if (c) c->value = val;
    }

    // --- Métodos de Construção da Netlist ---
    void add_resistor(std::string name, std::string n1, std::string n2, double val) {
        nl.add_resistor(name, val, { n1, n2 });
    }
    void add_capacitor(std::string name, std::string n1, std::string n2, double val) {
        nl.add_capacitor(name, val, { n1, n2 });
    }
    void add_inductor(std::string name, std::string n1, std::string n2, double val) {
        nl.add_inductor(name, val, { n1, n2 });
    }
    void add_dc_source(std::string name, std::string p, std::string n, double val) {
        nl.add_voltage_source(name, val, { p, n }, "DC");
    }
    void add_mosfet(std::string name, std::string d, std::string s) {
        nl.add_controlled_switch(name, { d, s }, false);
        // Adiciona diodo antiparalelo automaticamente
        nl.add_switch("D_" + name, { s, d }, false);
    }

    void set_pwm(std::string sw_name, double freq, double duty, double phase, double dead_time) {
        std::string key = upper(sw_name);
        pwm_config[key + "_FREQ"] = freq;
        pwm_config[key + "_DUTY"] = duty;
        pwm_config[key + "_PHASE"] = phase;
        pwm_config[key + "_DEAD"] = dead_time;
    }

    void add_diode(std::string name, std::string a, std::string k) {
        nl.add_switch(name, { a, k }, false);
    }
    void add_transformer(std::string name, std::string p1, std::string p2, std::string s1, std::string s2, double ratio) {
        nl.add_transformer(name, ratio, { p1, p2, s1, s2 });
    }

    // ==============================================================================
    //  MÉTODO RUN (Principal)
    // ==============================================================================
    py::dict run(std::string mode, double dt, double t_end, double period,
        std::vector<std::string> probes, bool export_matrices = false) {

        Netlist work_nl = this->nl; // Copia a topologia base para não alterar a original

        // Carrega o estado anterior se existir (Warm Start)
        if (has_state && last_state.size() == work_nl.total_mna_size()) {
            Solver::set_state_vector(work_nl, last_state);
        }

        // --- 1. Configuração de Probes e Instrumentação ---
        std::map<std::string, std::string> probe_map;
        std::set<std::string> components_to_probe;
        for (const auto& p : probes) {
            std::string U = upper(p);
            // Se for corrente I(Comp), marca o componente
            if (U.size() > 3 && U[0] == 'I' && U[1] == '(')
                components_to_probe.insert(p.substr(2, p.size() - 3));
        }

        // Insere fontes de tensão 0V em série para medir corrente onde necessário
        struct InstrumentTask { std::string comp_name, src_name, new_node, old_node; };
        std::vector<InstrumentTask> tasks;
        for (auto& comp : work_nl.get_components()) {
            if (components_to_probe.count(comp.name)) {
                // Fontes de tensão e Trafos já medem corrente nativamente no MNA
                if (comp.ctype == "V" || comp.ctype == "TF") {
                    probe_map[comp.name] = comp.name;
                    continue;
                }
                if (comp.nodes.size() < 2) continue;

                // Para R, L, C, Switch: insere V_sense em série
                InstrumentTask task = { comp.name, "V_sense_" + comp.name, comp.name + "_sense", comp.nodes[1] };
                tasks.push_back(task);
                probe_map[comp.name] = task.src_name;
            }
        }
        for (const auto& task : tasks) {
            auto* c = work_nl.find_component(task.comp_name);
            if (c) {
                c->nodes[1] = task.new_node; // Desconecta nó original
                work_nl.add_voltage_source(task.src_name, 0.0, { task.new_node, task.old_node }, "DC");
            }
        }

        // Mapeamento de índices para leitura rápida
        auto node_map = work_nl.node_index();
        std::unordered_map<std::string, int> v_src_idx;
        int nc = work_nl.count_nodes_excluding_gnd();
        int v_i = 0;
        for (const auto& c : work_nl.get_components()) {
            if (c.ctype == "V" || c.ctype == "TF") v_src_idx[c.name] = nc + v_i++;
        }

        struct FastProbe {
            std::string name; int idx; std::vector<double> buf; Waveform stats;
            FastProbe(std::string n, int i) : name(n), idx(i), stats(n) {}
        };

        std::vector<FastProbe> fast_probes;
        size_t expected_points = (t_end > 0 && dt > 0) ? (size_t)(t_end / dt) + 10 : 100;

        for (const auto& p : probes) {
            int idx = -1;
            std::string U = upper(p);
            if (U.find("V(") == 0) { // Tensão Nodal V(no)
                std::string n = p.substr(2, p.size() - 3);
                std::transform(n.begin(), n.end(), n.begin(), ::tolower);
                if (n == "ground") n = "gnd";
                if (node_map.count(n)) idx = node_map.at(n) - 1;
            }
            else if (U.find("I(") == 0) { // Corrente no componente I(Comp)
                std::string c = p.substr(2, p.size() - 3);
                std::string src = probe_map.count(c) ? probe_map[c] : c;
                if (v_src_idx.count(src)) idx = v_src_idx[src];
            }
            FastProbe fp(p, idx);
            if (expected_points > 0) fp.buf.reserve(expected_points);
            fast_probes.push_back(std::move(fp));
        }

        // --- 2. Configuração PWM Otimizada ---
        std::vector<FastPWMController> active_pwms;
        auto sws = work_nl.get_switches();
        for (auto* sw : sws) {
            std::string k = upper(sw->name);
            if (sw->is_controlled && pwm_config.count(k + "_FREQ")) {
                double f = pwm_config[k + "_FREQ"];
                if (f > 0) {
                    FastPWMController ctrl;
                    ctrl.sw = sw;
                    ctrl.Ts = 1.0 / f;
                    double duty = pwm_config[k + "_DUTY"];
                    double phase = pwm_config[k + "_PHASE"];
                    double dead = pwm_config[k + "_DEAD"];

                    ctrl.on_time = duty * ctrl.Ts - dead;
                    ctrl.phase_offset = std::fmod(phase * ctrl.Ts, ctrl.Ts);
                    ctrl.dead_time = dead;
                    active_pwms.push_back(ctrl);
                }
            }
        }

        // Lambda de controle que será chamado a cada passo de tempo pelo Solver
        auto control_optimized = [&](double t, Netlist& /*unused*/) {
            for (const auto& ctrl : active_pwms) {
                double local_t = std::fmod(t + ctrl.Ts - ctrl.phase_offset, ctrl.Ts);
                if (local_t < 0) local_t += ctrl.Ts;
                ctrl.sw->gate_signal = (local_t < ctrl.on_time);
            }
            };

        // --- 3. Execução da Simulação ---
        std::vector<double> time_vec;
        if (expected_points > 0) time_vec.reserve(expected_points);
        Vector state;

        // Estrutura para capturar matrizes se solicitado
        SolverDebugInfo debug_info;
        SolverDebugInfo* p_debug = export_matrices ? &debug_info : nullptr;

        {
            // Libera o GIL do Python para permitir multithreading se houver
            py::gil_scoped_release release;

            if (mode == "shooting") {
                // No Shooting Method, usamos um dt fino para precisão, mas limitado pelo periodo
                double shoot_dt = (dt > 0) ? dt : period / 200.0;

                state = Solver::find_steady_state_with_switches(
                    work_nl, shoot_dt, period, control_optimized,
                    this->solver_max_iter, this->solver_tol,
                    this->solver_eps, // Passa o EPS configurado no Python
                    p_debug,          // Passa o ponteiro para extrair matrizes
                    false             // Verbose off
                );
            }
            else { // Modo Transient (default)
                if (!has_state) {
                    // Inicialização a frio
                    work_nl.set_time(0.0);
                    control_optimized(0.0, work_nl);
                    state = Solver::get_state_vector(work_nl);
                }
                else {
                    state = last_state;
                }
            }

            // Atualiza estado persistente
            this->last_state = state;
            this->has_state = true;

            // Se t_end > 0, roda uma simulação no tempo para gerar formas de onda
            if (t_end > 0) {
                // Reseta tempo para 0 para facilitar plotagem 
                work_nl.set_time(0.0);
                Solver::set_state_vector(work_nl, state);

                Solver::simulate_time(work_nl, state, dt, t_end, [&](double t, const Vector& sol) {
                    // 1. Executa controle
                    control_optimized(t, work_nl);

                    // 2. Salva dados
                    time_vec.push_back(t);
                    for (auto& fp : fast_probes) {
                        double v = 0.0;
                        if (fp.idx >= 0 && fp.idx < sol.size()) v = sol(fp.idx);
                        fp.buf.push_back(v);
                        fp.stats.add(t, v);
                    }
                    });

                // Atualiza o last_state com o final da transiente também
                this->last_state = Solver::get_state_vector(work_nl);
            }
        }

        // --- 4. Empacotamento dos Resultados ---
        py::dict res, sig, sta;

        // Vetor de tempo
        if (!time_vec.empty()) res["t"] = time_vec;

        // Sinais e Estatísticas
        for (const auto& fp : fast_probes) {
            if (!fp.buf.empty()) sig[fp.name.c_str()] = fp.buf;

            py::dict s;
            s["mean"] = fp.stats.mean();
            s["rms"] = fp.stats.rms();
            s["max"] = fp.stats.max();
            s["min"] = fp.stats.min();
            s["ripple"] = fp.stats.ripple();
            sta[fp.name.c_str()] = s;
        }
        res["signals"] = sig;
        res["stats"] = sta;

        // Exportação das matrizes (Se solicitado e disponível)
        if (export_matrices && debug_info.captured) {
            res["jacobian"] = debug_info.J; // Converte Eigen::Matrix -> Numpy Array automaticamente
            res["inverse_jacobian"] = debug_info.B;
        }

        return res;
    }
};

// ==============================================================================
//  DEFINIÇÃO DO MÓDULO PYTHON
// ==============================================================================
PYBIND11_MODULE(_core, m) {
    m.doc() = "Fast SPICE Core Backend";

    py::class_<PyCircuit>(m, "Circuit")
        .def(py::init<>())
        .def("set_params", &PyCircuit::set_params,
            py::arg("max_iter"), py::arg("tol"), py::arg("eps") = 1e-6)

        .def("reset_state", &PyCircuit::reset_state)
        .def("update_value", &PyCircuit::update_value)

        // variáveis internas 
        .def_readwrite("last_state", &PyCircuit::last_state)
        .def_readwrite("has_state", &PyCircuit::has_state)
        // ----------------------------------------------------

        // Componentes Passivos
        .def("add_resistor", &PyCircuit::add_resistor)
        .def("add_capacitor", &PyCircuit::add_capacitor)
        .def("add_inductor", &PyCircuit::add_inductor)

        // Fontes e Chaves
        .def("add_dc_source", &PyCircuit::add_dc_source)
        .def("add_mosfet", &PyCircuit::add_mosfet)
        .def("add_diode", &PyCircuit::add_diode)
        .def("add_transformer", &PyCircuit::add_transformer)

        // PWM
        .def("set_pwm", &PyCircuit::set_pwm,
            py::arg("sw_name"), py::arg("freq"), py::arg("duty"),
            py::arg("phase"), py::arg("dead_time"))

        // Run
        .def("run", &PyCircuit::run,
            py::arg("mode"), py::arg("dt"), py::arg("t_end"),
            py::arg("period") = 0.0,
            py::arg("probes") = std::vector<std::string>(),
            py::arg("export_matrices") = false);
}