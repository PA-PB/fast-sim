// Solver.h
#ifndef SOLVER_H
#define SOLVER_H


#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 26495 26498 26813 6294 6255)
#endif

#include <Eigen/Dense>
#include <Eigen/LU>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "Netlist.h"
#include "Snapshot.h"
#include <map>
#include <vector>
#include <tuple>
#include <functional>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using SwitchState = std::vector<bool>;

using GDecomp = Eigen::FullPivLU<Matrix>;
struct LutEntry {
    GDecomp lu;
    bool valid = false;
};

// --- NOVA ESTRUTURA PARA EXPORTAR MATRIZES ---
struct SolverDebugInfo {
    Matrix J; // Jacobiana Final
    Matrix B; // Inversa da Jacobiana (Broyden/Newton)
    bool captured = false;
};

class Solver {
public:

    static Matrix build_G_for_state(Netlist& netlist, double dt,
        const std::map<std::string, bool>& switch_state_map);

    static std::pair<std::map<SwitchState, LutEntry>, std::vector<Netlist::Component*>>
        precompute_matrices(Netlist& netlist, double dt);

    static Vector transient_Z_matrix(Netlist& netlist, double dt, double current_time = 0);

    static void update_state(Netlist& netlist, const Vector& solution_vector,
        const std::unordered_map<std::string, int>& node_map,
        double dt);

    static Vector get_state_vector(Netlist& netlist);
    static void set_state_vector(Netlist& netlist, const Vector& states);

    static double get_v(const std::string& node_name, const Vector& sol,
        const std::unordered_map<std::string, int>& node_map);

    // --- MÉTODOS DE SHOOTING ATUALIZADOS ---

    // Versão Simples (mantida para compatibilidade)
    static Vector find_steady_state(Netlist& netlist, double dt, double period,
        int max_iters = 20, double tol = 1e-4, bool verbose = false);

    // Versão Principal: Aceita eps e debug_info
    static Vector find_steady_state_with_switches(
        Netlist& netlist,
        double dt,
        double period,
        std::function<void(double t, Netlist& nl)> control_logic,
        int max_iters = 20,
        double tol = 1e-4,
        double eps = 1e-4,              // <--- Controle de Epsilon
        SolverDebugInfo* debug_out = nullptr, // <--- Opcional para exportar matrizes
        bool verbose = false
    );

    static Snapshot find_steady_state_snapshot(Netlist& netlist, double dt, double period,
        std::function<void(double t, Netlist& nl)> control_logic,
        int max_iters = 20, double tol = 1e-4, bool verbose = false) {

        find_steady_state_with_switches(netlist, dt, period, control_logic, max_iters, tol, 1e-4, nullptr, verbose);
        return Snapshot::from_netlist(netlist);
    }

    static void simulate_time(
        Netlist& netlist,
        const Vector& start_state,
        double dt,
        double total_time,
        std::function<void(double t, const Vector& sol)> callback
    );

private:
    static int row(int n) { return n == 0 ? -1 : n - 1; }
};

#endif // SOLVER_H