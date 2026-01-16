#include "Netlist.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

Netlist::Netlist() : voltage_source_count(0) {}

void Netlist::register_nodes(const std::vector<std::string>& nds) {
    for (const auto& node : nds) {
        nodes.insert(node);
    }
}

void Netlist::add_resistor(const std::string& name, double value,
    const std::vector<std::string>& nds) {
    components.emplace_back("R", name, value, nds);
    register_nodes(nds);
}

void Netlist::add_capacitor(const std::string& name, double value,
    const std::vector<std::string>& nds) {
    components.emplace_back("C", name, value, nds);
    register_nodes(nds);
}

void Netlist::add_inductor(const std::string& name, double value,
    const std::vector<std::string>& nds) {
    components.emplace_back("L", name, value, nds);
    register_nodes(nds);
}

void Netlist::add_voltage_source(const std::string& name, double value,
    const std::vector<std::string>& nds,
    const std::string& type, double freq) {
    Component comp("V", name, value, nds);
    comp.type = type;
    comp.freq = freq;
    components.push_back(comp);
    register_nodes(nds);
    voltage_source_count++;
}

void Netlist::add_current_source(const std::string& name, double value,
    const std::vector<std::string>& nds) {
    components.emplace_back("I", name, value, nds);
    register_nodes(nds);
}

// Chave convencional (Diodo)
void Netlist::add_switch(const std::string& name,
    const std::vector<std::string>& nds,
    bool initial_state) {
    Component comp("S", name, 0, nds);
    comp.is_on = initial_state;
    comp.is_controlled = false; // Diodo
    components.push_back(comp);
    register_nodes(nds);
}

// Chave controlada (MOSFET/IGBT)
void Netlist::add_controlled_switch(const std::string& name,
    const std::vector<std::string>& nds,
    bool initial_state) {
    Component comp("S", name, 0, nds);
    comp.is_on = initial_state;
    comp.is_controlled = true;       // Controlada
    comp.gate_signal = initial_state;
    components.push_back(comp);
    register_nodes(nds);
}

void Netlist::add_transformer(const std::string& name, double turn_ratio,
    const std::vector<std::string>& nds) {
    components.emplace_back("TF", name, turn_ratio, nds);
    register_nodes(nds);
    voltage_source_count++;
}

std::vector<Netlist::Component*> Netlist::get_switches() {
    std::vector<Component*> switches;
    for (auto& c : components) {
        if (c.ctype == "S") {
            switches.push_back(&c);
        }
    }
    return switches;
}

int Netlist::count_nodes_excluding_gnd() const {
    int count = 0;
    for (const auto& n : nodes) {
        std::string lower = n;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != "gnd" && lower != "0") {
            count++;
        }
    }
    return count;
}

int Netlist::total_mna_size() const {
    return count_nodes_excluding_gnd() + voltage_source_count;
}

std::unordered_map<std::string, int> Netlist::node_index() const {
    std::unordered_map<std::string, int> mapping;
    std::vector<std::string> sorted_nodes(nodes.begin(), nodes.end());
    std::sort(sorted_nodes.begin(), sorted_nodes.end());

    int idx = 1;
    for (const auto& n : sorted_nodes) {
        std::string lower = n;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower != "gnd" && lower != "0") {
            mapping[n] = idx++;
        }
        else {
            mapping[n] = 0;
        }
    }
    return mapping;
}


static std::string trim(std::string s) {
    auto notspace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
    return s;
}

static std::string upper(std::string s) {
    for (auto& c : s) c = (char)std::toupper((unsigned char)c);
    return s;
}

static std::string lower(std::string s) {
    for (auto& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

static std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> out;
    std::string t;
    while (iss >> t) out.push_back(t);
    return out;
}

static std::string norm_node(std::string s) {
    s = trim(s);
    s = lower(s);
    if (s == "ground") s = "gnd";
    return s;
}

static bool starts_with(const std::string& s, const std::string& p) {
    return s.size() >= p.size() && std::equal(p.begin(), p.end(), s.begin());
}

static double parse_spice_number(const std::string& tok) {
    std::string s = tok;
    s.erase(std::remove(s.begin(), s.end(), '_'), s.end());
    if (s.empty()) return 0.0;

    std::string sufx;
    std::string base = s;

    auto u = upper(s);
    if (u.size() >= 3 && u.substr(u.size() - 3) == "MEG") {
        sufx = "MEG";
        base = s.substr(0, s.size() - 3);
    }
    else {
        char last = u.back();
        if (std::isalpha((unsigned char)last)) {
            sufx = std::string(1, last);
            base = s.substr(0, s.size() - 1);
        }
    }

    double v = 0.0;
    try { v = std::stod(base); }
    catch (...) { v = 0.0; }

    double m = 1.0;
    if (sufx == "T")    m = 1e12;
    else if (sufx == "G")    m = 1e9;
    else if (sufx == "MEG") m = 1e6;
    else if (sufx == "K")    m = 1e3;
    else if (sufx == "M")    m = 1e-3;
    else if (sufx == "U")    m = 1e-6;
    else if (sufx == "N")    m = 1e-9;
    else if (sufx == "P")    m = 1e-12;
    else if (sufx == "F")    m = 1e-15;

    return v * m;
}


Netlist build_netlist_from_text(const std::string& text) {
    Netlist nl;
    std::istringstream iss(text);
    std::string raw;


    while (std::getline(iss, raw)) {
        std::string line = trim(raw);
        if (line.empty()) continue;
        if (starts_with(line, "*") || starts_with(line, ";") || starts_with(line, "//")) continue;
        if (line[0] == '.') continue;

        auto toks = split_ws(line);
        if (toks.empty()) continue;

        std::string name = toks[0];
        char c0 = (char)std::toupper((unsigned char)name[0]);
        auto node2 = [&](int i) { return toks.size() > (size_t)i ? norm_node(toks[i]) : ""; };

        try {
            if (c0 == 'R' && toks.size() >= 4) {
                nl.add_resistor(name, parse_spice_number(toks[3]), { node2(1), node2(2) });
            }
            else if (c0 == 'C' && toks.size() >= 4) {
                nl.add_capacitor(name, parse_spice_number(toks[3]), { node2(1), node2(2) });
            }
            else if (c0 == 'L' && toks.size() >= 4) {
                nl.add_inductor(name, parse_spice_number(toks[3]), { node2(1), node2(2) });
            }
            else if (c0 == 'V' && toks.size() >= 4) {
                std::string type = "DC";
                double freq = 0.0;
                double value = 0.0;
                if (toks.size() >= 5 && upper(toks[3]) == "DC") {
                    value = parse_spice_number(toks[4]);
                }
                else {
                    value = parse_spice_number(toks[3]);
                }
                nl.add_voltage_source(name, value, { node2(1), node2(2) }, type, freq);
            }
            else if ((c0 == 'S' || c0 == 'Q' || c0 == 'D') && toks.size() >= 3) {
                bool is_controlled = false;
                if (upper(line).find("CONTROLLED") != std::string::npos) is_controlled = true;

                // Hack simples para detectar MOSFETs pelo nome (Qxxx)
                if (c0 == 'Q') is_controlled = true;

                if (is_controlled) {
                    nl.add_controlled_switch(name, { node2(1), node2(2) }, false);
                    auto* c = nl.find_component(name);
                    if (c) { c->is_controlled = true; c->gate_signal = false; }
                }
                else {
                    nl.add_switch(name, { node2(1), node2(2) }, false);
                }
            }
            else if (c0 == 'T' && toks.size() >= 5) {
                double ratio = parse_spice_number(toks.back());
                nl.add_transformer(name, ratio, { node2(1), node2(2), node2(3), node2(4) });
            }
        }
        catch (...) {
            // Ignora erros de parse para não quebrar tudo
        }
    }
    return nl;
}