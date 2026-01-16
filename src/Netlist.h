#ifndef NETLIST_H
#define NETLIST_H

#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>

class Netlist {

public:
    double get_time() const { return sim_time; }
    void set_time(double t) { sim_time = t; }
private:
    double sim_time = 0.0;

public:
    struct Component {
        std::string ctype;
        std::string name;
        double value;
        std::vector<std::string> nodes;

        // Parâmetros adicionais
        std::string type;  // DC, SQUARE, etc
        double freq;

        // Memória para simulação transiente
        double v_prev;
        double i_prev;

        // Estado lógico (switches)
        bool is_on;

        // --- NOVOS CAMPOS PARA CONTROLE ---
        bool is_controlled; // true = Controlada (MOSFET/IGBT), false = Diodo
        bool gate_signal;   // Sinal de comando externo (PWM)

        Component(const std::string& ct, const std::string& n, double v,
            const std::vector<std::string>& nds)
            : ctype(ct), name(n), value(v), nodes(nds),
            type("DC"), freq(0), v_prev(0), i_prev(0), is_on(false),
            is_controlled(false), gate_signal(false) {
        }
    };

private:
    std::vector<Component> components;
    std::set<std::string> nodes;
    int voltage_source_count;

    void register_nodes(const std::vector<std::string>& nds);

public:
    Netlist();

    void add_resistor(const std::string& name, double value,
        const std::vector<std::string>& nodes);

    void add_capacitor(const std::string& name, double value,
        const std::vector<std::string>& nodes);

    void add_inductor(const std::string& name, double value,
        const std::vector<std::string>& nodes);

    void add_voltage_source(const std::string& name, double value,
        const std::vector<std::string>& nodes,
        const std::string& type = "DC", double freq = 0);

    void add_current_source(const std::string& name, double value,
        const std::vector<std::string>& nodes);

    // Adiciona chave não controlada (Diodo ideal)
    void add_switch(const std::string& name,
        const std::vector<std::string>& nodes,
        bool initial_state = false);

    // Adiciona chave controlada (MOSFET/IGBT ideal) 
    void add_controlled_switch(const std::string& name,
        const std::vector<std::string>& nodes,
        bool initial_state = false);

    void add_transformer(const std::string& name, double turn_ratio,
        const std::vector<std::string>& nodes);

    Component* find_component(const std::string& name) {
        for (auto& c : components)
            if (c.name == name)
                return &c;
        return nullptr;
    }

    std::vector<Component*> get_switches();
    int count_nodes_excluding_gnd() const;
    int total_mna_size() const;
    std::unordered_map<std::string, int> node_index() const;

    const std::vector<Component>& get_components() const { return components; }
    std::vector<Component>& get_components() { return components; }
};

#endif // NETLIST_H