// Solver.cpp

#define EIGEN_NO_DEBUG
#define NDEBUG

#define _USE_MATH_DEFINES

#include "Solver.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <omp.h> 




#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

// ==============================================================================
//  ESTRUTURAS DE DADOS OTIMIZADAS (PODs) - ESCOPO GLOBAL
// ==============================================================================

struct FastCap {
    int node1, node2;
    double g_eq;
    double* v_prev_ptr;
};

struct FastInd {
    int node1, node2;
    double val, dt_div_L;
    double* i_prev_ptr;
};

struct FastSource {
    int matrix_idx, type;
    double val, freq;
};

struct SwitchInfo {
    int n1_idx, n2_idx;
    bool is_controlled;
    bool* gate_ptr;
};

struct FastLutEntry {
    Matrix inv;
    bool valid = false;
};

struct FastSimulationContext {
    std::vector<FastCap> caps;
    std::vector<FastInd> inds;
    std::vector<FastSource> sources;
    Vector Z, sol;
};

// ==============================================================================
//  AUXILIARES E HELPERS (ESCOPO GLOBAL)
// ==============================================================================

static int get_switch_mask(const std::vector<Netlist::Component*>& switches) {
    int mask = 0;
    for (size_t i = 0; i < switches.size(); i++) {
        if (switches[i]->is_on) mask |= (1 << i);
    }
    return mask;
}

inline bool get_next_switch_state(bool was_on, double vd) {
    constexpr double Ron = 1e-3;
    constexpr double V_THRESHOLD_ON = 0.5;
    constexpr double I_THRESHOLD_OFF = -10e-3;
    if (was_on) return (vd / Ron) > I_THRESHOLD_OFF;
    else return vd > V_THRESHOLD_ON;
}

// Helpers de leitura rápida
inline double get_v_direct(int idx, const Vector& sol) {
    if (idx < 0) return 0.0;
    return sol(idx);
}

// Compila a netlist para estruturas rápidas (Definida globalmente)
FastSimulationContext compile_netlist(Netlist& netlist, double dt) {
    FastSimulationContext ctx;
    auto node_map = netlist.node_index();
    int node_eqs = netlist.count_nodes_excluding_gnd();
    int size = netlist.total_mna_size();

    ctx.Z = Vector::Zero(size);
    ctx.sol = Vector::Zero(size);

    int vsrc_counter = 0;
    for (auto& c : netlist.get_components()) {
        int n1 = -1, n2 = -1;
        if (c.nodes.size() >= 2) {
            if (node_map.count(c.nodes[0]) && c.nodes[0] != "gnd" && c.nodes[0] != "0")
                n1 = node_map.at(c.nodes[0]) - 1;
            if (node_map.count(c.nodes[1]) && c.nodes[1] != "gnd" && c.nodes[1] != "0")
                n2 = node_map.at(c.nodes[1]) - 1;
        }

        if (c.ctype == "C") ctx.caps.push_back({ n1, n2, c.value / dt, &c.v_prev });
        else if (c.ctype == "L") ctx.inds.push_back({ n1, n2, c.value, dt / c.value, &c.i_prev });
        else if (c.ctype == "V") {
            int k = node_eqs + vsrc_counter;
            int t = (c.type == "SQUARE") ? 1 : 0;
            ctx.sources.push_back({ k, t, c.value, c.freq });
            vsrc_counter++;
        }
        else if (c.ctype == "TF") vsrc_counter++;
    }
    return ctx;
}

// Preenche vetor Z (Definida globalmente)
inline void fast_fill_Z(FastSimulationContext& ctx, double t) {
    ctx.Z.setZero();
    for (const auto& src : ctx.sources) {
        double val = src.val;
        if (src.type == 1 && src.freq > 0) {
            double period = 1.0 / src.freq;
            double phase = fmod(t, period);
            val = (phase < period / 2) ? src.val : 0.0;
        }
        ctx.Z(src.matrix_idx) = val;
    }
    for (const auto& c : ctx.caps) {
        double i_hist = c.g_eq * (*c.v_prev_ptr);
        if (c.node1 >= 0) ctx.Z(c.node1) += i_hist;
        if (c.node2 >= 0) ctx.Z(c.node2) -= i_hist;
    }
    for (const auto& l : ctx.inds) {
        double i_hist = (*l.i_prev_ptr);
        if (l.node1 >= 0) ctx.Z(l.node1) -= i_hist;
        if (l.node2 >= 0) ctx.Z(l.node2) += i_hist;
    }
}

// ==============================================================================
//  PRE-COMPUTAÇÃO
// ==============================================================================

Matrix Solver::build_G_for_state(Netlist& netlist, double dt,
    const std::map<std::string, bool>& switch_state_map) {

    int size = netlist.total_mna_size();
    Matrix G = Matrix::Zero(size, size);
    auto node_index = netlist.node_index();
    int node_equations = netlist.count_nodes_excluding_gnd();
    int vsrc_counter = 0;
    const double gmin = 1e-9;

    for (int r = 0; r < node_equations; r++) G(r, r) += gmin;

    for (const auto& c : netlist.get_components()) {
        int n1 = -1, n2 = -1, r1 = -1, r2 = -1;
        if (c.nodes.size() >= 2) {
            if (node_index.count(c.nodes[0])) n1 = node_index.at(c.nodes[0]);
            if (node_index.count(c.nodes[1])) n2 = node_index.at(c.nodes[1]);
            r1 = (n1 == 0 ? -1 : n1 - 1);
            r2 = (n2 == 0 ? -1 : n2 - 1);
        }

        double g = 0.0;
        bool is_conductive = false;

        if (c.ctype == "S") {
            bool is_on = false;
            auto it = switch_state_map.find(c.name);
            if (it != switch_state_map.end()) is_on = it->second;
            g = is_on ? 1000.0 : 1e-8;
            is_conductive = true;
        }
        else if (c.ctype == "R") { g = 1.0 / c.value; is_conductive = true; }
        else if (c.ctype == "C") { g = c.value / dt; is_conductive = true; }
        else if (c.ctype == "L") { g = dt / c.value; is_conductive = true; }

        if (is_conductive) {
            if (r1 >= 0) G(r1, r1) += g;
            if (r2 >= 0) G(r2, r2) += g;
            if (r1 >= 0 && r2 >= 0) { G(r1, r2) -= g; G(r2, r1) -= g; }
        }
        else if (c.ctype == "V") {
            int k = node_equations + vsrc_counter;
            if (r1 >= 0) { G(k, r1) = 1; G(r1, k) = 1; }
            if (r2 >= 0) { G(k, r2) = -1; G(r2, k) = -1; }
            vsrc_counter++;
        }
        else if (c.ctype == "TF") {
            int np_p = node_index.at(c.nodes[0]), np_m = node_index.at(c.nodes[1]);
            int ns_p = node_index.at(c.nodes[2]), ns_m = node_index.at(c.nodes[3]);
            int rp_p = (np_p == 0 ? -1 : np_p - 1), rp_m = (np_m == 0 ? -1 : np_m - 1);
            int rs_p = (ns_p == 0 ? -1 : ns_p - 1), rs_m = (ns_m == 0 ? -1 : ns_m - 1);
            int k = node_equations + vsrc_counter;
            double n = c.value;

            if (rp_p >= 0) G(k, rp_p) = 1; if (rp_m >= 0) G(k, rp_m) = -1; if (rs_p >= 0) G(k, rs_p) = -n; if (rs_m >= 0) G(k, rs_m) = n;
            if (rp_p >= 0) G(rp_p, k) = 1; if (rp_m >= 0) G(rp_m, k) = -1; if (rs_p >= 0) G(rs_p, k) = -n; if (rs_m >= 0) G(rs_m, k) = n;
            vsrc_counter++;
        }
    }
    return G;
}

static std::vector<FastLutEntry> precompute_matrices_inverse(Netlist& netlist, double dt,
    std::vector<Netlist::Component*>& switches) {

        
    int num_switches = static_cast<int>(switches.size()); 
    int total_states = 1 << num_switches;
    std::vector<FastLutEntry> lut(total_states);

#pragma omp parallel for 
    for (int i = 0; i < total_states; i++) {
        std::map<std::string, bool> state_map;
        for (int j = 0; j < num_switches; j++) {
            state_map[switches[j]->name] = (i >> j) & 1;
        }
        Matrix G = Solver::build_G_for_state(netlist, dt, state_map);
        Eigen::FullPivLU<Matrix> lu(G);
        if (lu.isInvertible()) {
            lut[i].inv = lu.inverse();
            lut[i].valid = true;
        }
        else {
            lut[i].valid = false;
        }
    }
    return lut;
}


std::pair<std::map<SwitchState, LutEntry>,
std::vector<Netlist::Component*>>
Solver::precompute_matrices(Netlist& n, double dt) { return {}; }

Vector Solver::transient_Z_matrix(Netlist& n, double dt, double t) { return Vector(); }
void Solver::update_state(Netlist& n, const Vector& s, const std::unordered_map<std::string, int>& m, double dt) {}
Vector Solver::find_steady_state(Netlist& n, double dt, double p, int i, double t, bool v) {
    return find_steady_state_with_switches(n, dt, p, nullptr, i, t, 1e-6, nullptr, v);
}

// Get/Set States
Vector Solver::get_state_vector(Netlist& netlist) {
    std::vector<double> s;
    for (auto& c : netlist.get_components()) {
        if (c.ctype == "C") s.push_back(c.v_prev);
        else if (c.ctype == "L") s.push_back(c.i_prev);
    }
    Vector res(s.size());
    for (size_t i = 0; i < s.size(); i++) res(i) = s[i];
    return res;
}

void Solver::set_state_vector(Netlist& netlist, const Vector& states) {
    int idx = 0;
    for (auto& c : netlist.get_components()) {
        if (c.ctype == "C") c.v_prev = states(idx++);
        else if (c.ctype == "L") c.i_prev = states(idx++);
    }
}

double Solver::get_v(const std::string& node_name, const Vector& sol,
    const std::unordered_map<std::string, int>& node_map) {
    if (node_name == "gnd" || node_name == "0") return 0.0;
    auto it = node_map.find(node_name);
    if (it == node_map.end()) return 0.0;
    int idx = it->second;
    if (idx - 1 >= sol.size()) return 0.0;
    return sol(idx - 1);
}

// ==============================================================================
//  SHOOTING METHOD OTIMIZADO COM CONTROLE E DEBUG
// ==============================================================================

Vector Solver::find_steady_state_with_switches(
    Netlist& netlist,
    double dt,
    double period,
    std::function<void(double t, Netlist& nl)> control_logic,
    int max_iters,
    double tol,
    double eps,                 // Parâmetro de sensibilidade
    SolverDebugInfo* debug_out, // Ponteiro opcional de Debug
    bool verbose)
{
    if (verbose) std::cout << "Iniciando Fast Shooting Method..." << std::endl;

    auto switches = netlist.get_switches();
    auto lut_vec = precompute_matrices_inverse(netlist, dt, switches);
    auto node_map = netlist.node_index();

    std::vector<SwitchInfo> sw_infos(switches.size());
    for (size_t i = 0; i < switches.size(); i++) {
        std::string n1 = switches[i]->nodes[0];
        std::string n2 = switches[i]->nodes[1];
        int idx1 = (n1 != "gnd" && n1 != "0" && node_map.count(n1)) ? node_map.at(n1) - 1 : -1;
        int idx2 = (n2 != "gnd" && n2 != "0" && node_map.count(n2)) ? node_map.at(n2) - 1 : -1;

        sw_infos[i].n1_idx = idx1;
        sw_infos[i].n2_idx = idx2;
        sw_infos[i].is_controlled = switches[i]->is_controlled;
        sw_infos[i].gate_ptr = &switches[i]->gate_signal;
    }

    FastSimulationContext ctx = compile_netlist(netlist, dt);

    int steps = static_cast<int>(period / dt);
    Vector x_k = get_state_vector(netlist);
    int dum_mask = 0;

    // --- ENGINE INTERNA: Roda 1 Período ---
    auto simulate_period = [&](const Vector& start_x, int& out_mask) -> Vector {
        set_state_vector(netlist, start_x);
        int curr_mask = 0;
        // Recalcula mask inicial baseada no estado atual das chaves no netlist
        for (size_t i = 0; i < switches.size(); i++) {
            if (switches[i]->is_on) curr_mask |= (1 << i);
        }

        double t_base = 0.0; // Shooting assume regime permanente periódico (t relativo)

        for (int k = 0; k < steps; k++) {
            double t = t_base + k * dt;

            // 1. Injeta Lógica de Controle
            if (control_logic) control_logic(t, netlist);

            // 2. Preenche Z (Chama função global)
            fast_fill_Z(ctx, t);

            // 3. Resolve Topologia
            int switch_iter = 0;
            bool topology_changed = true;
            while (topology_changed && switch_iter < 8) {
                if (curr_mask < lut_vec.size() && lut_vec[curr_mask].valid) {
                    ctx.sol.noalias() = lut_vec[curr_mask].inv * ctx.Z;
                }
                else {
                    ctx.sol.setZero();
                }

                int next_mask = 0;
                for (size_t idx = 0; idx < switches.size(); idx++) {
                    if (sw_infos[idx].is_controlled) {
                        // Segue o Gate
                        if (*(sw_infos[idx].gate_ptr)) next_mask |= (1 << idx);
                    }
                    else {
                        // Segue o Diodo (V/I) - Usa helpers globais
                        double v1 = (sw_infos[idx].n1_idx >= 0) ? ctx.sol(sw_infos[idx].n1_idx) : 0.0;
                        double v2 = (sw_infos[idx].n2_idx >= 0) ? ctx.sol(sw_infos[idx].n2_idx) : 0.0;

                        bool was_on = (curr_mask >> idx) & 1;
                        // Lógica de histerese do diodo
                        bool next_on = false;
                        if (was_on) next_on = ((v1 - v2) / 1e-3) > -10e-3;
                        else next_on = (v1 - v2) > 0.5;

                        if (next_on) next_mask |= (1 << idx);
                    }
                }
                if (next_mask != curr_mask) { curr_mask = next_mask; topology_changed = true; }
                else { topology_changed = false; }
                switch_iter++;
            }

            // 4. Update de elementos reativos
            for (auto& c : ctx.caps) {
                double v1 = (c.node1 >= 0) ? ctx.sol(c.node1) : 0.0;
                double v2 = (c.node2 >= 0) ? ctx.sol(c.node2) : 0.0;
                *c.v_prev_ptr = v1 - v2;
            }
            for (auto& l : ctx.inds) {
                double v1 = (l.node1 >= 0) ? ctx.sol(l.node1) : 0.0;
                double v2 = (l.node2 >= 0) ? ctx.sol(l.node2) : 0.0;
                *l.i_prev_ptr += (v1 - v2) * l.dt_div_L;
            }
        }
        out_mask = curr_mask;
        return get_state_vector(netlist);
        };

    // --- NEWTON-RAPHSON ---
    if (verbose) std::cout << "Warm-up..." << std::endl;
    int dum;
    for (int i = 0; i < 10; i++) x_k = simulate_period(x_k, dum);

    auto F = [&](const Vector& x) -> Vector { return (simulate_period(x, dum) - x); };

    Vector Fx = F(x_k);
    if (Fx.norm() < tol) return x_k;

    int n_states = static_cast<int>(x_k.size());
    Matrix J = Matrix::Zero(n_states, n_states);

    // Jacobiano Inicial
    for (int j = 0; j < n_states; j++) {

        double current_val = std::abs(x_k(j));
        double scaled_eps = std::max(current_val * eps, 1e-7);

        Vector xp = x_k; xp(j) += scaled_eps;
        J.col(j) = (F(xp) - Fx) / scaled_eps;
    }
    Matrix B = J.fullPivLu().inverse();

    // Exporta Debug se solicitado
    if (debug_out) {
        debug_out->J = J;
        debug_out->B = B;
        debug_out->captured = true;
    }

    for (int i = 0; i < max_iters; i++) {
        Vector dx = -B * Fx;
        Vector x_new = x_k + dx;
        Vector Fx_new = F(x_new);
        bool divergence = (Fx_new.norm() > Fx.norm());

        // Line Search
        if (divergence) {
            double alpha = 1;
            for (int k = 0; k < 15; k++) {
                Vector trial = x_k + alpha * dx;
                Vector Ft = F(trial);
                if (Ft.norm() < Fx.norm()) { x_new = trial; Fx_new = Ft; divergence = false; break; }
                alpha *= 0.5;
            }
        }

        // Broyden Update ou Restart
        if (divergence) {
            if (verbose) std::cout << "  Divergencia. Recalc J..." << std::endl;
            for (int j = 0; j < n_states; j++) {

                double current_val = std::abs(x_k(j));
                double scaled_eps = std::max(current_val * eps, 1e-7);

                Vector xp = x_k; xp(j) += scaled_eps;
                J.col(j) = (F(xp) - Fx) / scaled_eps;
            }
            B = J.fullPivLu().inverse();

            // Atualiza Debug
            if (debug_out) { debug_out->J = J; debug_out->B = B; }

            dx = -B * Fx;
            x_new = x_k + 0.25 * dx;
            Fx_new = F(x_new);
        }
        else {
            Vector s = x_new - x_k;
            Vector y = Fx_new - Fx;
            Vector By = B * y;
            double den = s.dot(By);
            if (std::abs(den) > 1e-12) B = B + ((s - By) * (s.transpose() * B)) / den;
        }

        x_k = x_new; Fx = Fx_new;
        if (verbose) std::cout << "Iter " << i + 1 << " Err: " << Fx.norm() << std::endl;
        if (Fx.norm() < tol) break;
    }

    // Se saiu do loop, salva a última B aproximada se ainda não salvou a exata
    if (debug_out && !debug_out->captured) {
        debug_out->B = B;
    }

    Vector final_x = simulate_period(x_k, dum_mask);
    set_state_vector(netlist, final_x);
    for (size_t i = 0; i < switches.size(); i++) switches[i]->is_on = (dum_mask >> i) & 1;
    return final_x;
}

// ==============================================================================
//  SIMULATE TIME (TRANSIENTE COMUM)
// ==============================================================================

void Solver::simulate_time(
    Netlist& netlist, const Vector& start_state, double dt, double total_time,
    std::function<void(double, const Vector&)> callback) {

    auto switches = netlist.get_switches();
    auto lut_vec = precompute_matrices_inverse(netlist, dt, switches);
    auto node_map = netlist.node_index();

    // Chama compile_netlist global
    FastSimulationContext ctx = compile_netlist(netlist, dt);

    std::vector<SwitchInfo> sw_infos(switches.size());
    for (size_t i = 0; i < switches.size(); i++) {
        std::string n1 = switches[i]->nodes[0];
        std::string n2 = switches[i]->nodes[1];
        int idx1 = (n1 != "gnd" && n1 != "0" && node_map.count(n1)) ? node_map.at(n1) - 1 : -1;
        int idx2 = (n2 != "gnd" && n2 != "0" && node_map.count(n2)) ? node_map.at(n2) - 1 : -1;
        sw_infos[i].n1_idx = idx1; sw_infos[i].n2_idx = idx2;
        sw_infos[i].is_controlled = switches[i]->is_controlled;
        sw_infos[i].gate_ptr = &switches[i]->gate_signal;
    }

    set_state_vector(netlist, start_state);
    int curr_mask = get_switch_mask(switches);
    double t0 = netlist.get_time();
    int steps = static_cast<int>(total_time / dt);

    for (int k = 0; k < steps; k++) {
        double t = t0 + k * dt;
        fast_fill_Z(ctx, t); // Chama global

        int switch_iter = 0;
        bool topology_changed = true;
        while (topology_changed && switch_iter < 8) {
            if (curr_mask < lut_vec.size() && lut_vec[curr_mask].valid)
                ctx.sol.noalias() = lut_vec[curr_mask].inv * ctx.Z;
            else ctx.sol.setZero();

            int next_mask = 0;
            for (size_t idx = 0; idx < switches.size(); idx++) {
                if (sw_infos[idx].is_controlled) {
                    if (*(sw_infos[idx].gate_ptr)) next_mask |= (1 << idx);
                }
                else {
                    double v1 = get_v_direct(sw_infos[idx].n1_idx, ctx.sol);
                    double v2 = get_v_direct(sw_infos[idx].n2_idx, ctx.sol);
                    if (get_next_switch_state((curr_mask >> idx) & 1, v1 - v2)) next_mask |= (1 << idx);
                }
            }
            if (next_mask != curr_mask) { curr_mask = next_mask; topology_changed = true; }
            else topology_changed = false;
            switch_iter++;
        }

        callback(t, ctx.sol);

        for (auto& c : ctx.caps) {
            double v1 = get_v_direct(c.node1, ctx.sol);
            double v2 = get_v_direct(c.node2, ctx.sol);
            *c.v_prev_ptr = v1 - v2;
        }
        for (auto& l : ctx.inds) {
            double v1 = get_v_direct(l.node1, ctx.sol);
            double v2 = get_v_direct(l.node2, ctx.sol);
            *l.i_prev_ptr += (v1 - v2) * l.dt_div_L;
        }
    }
    for (size_t i = 0; i < switches.size(); i++) switches[i]->is_on = (curr_mask >> i) & 1;
    netlist.set_time(t0 + steps * dt);
}