#include "roundingsat.h"
#include "gaussian/constants.h"
#include <string>

using std::cout;
using std::endl;
using CMSat::Lit;
using CMSat::lbool;

namespace rs {

bool asynch_interrupt;
DLL_PUBLIC Options options;
DLL_PUBLIC Stats stats;

}  // namespace rs

int main(int argc, char** argv) {
    // parse arguments
    assert(argc > 1);
    string filename = argv[1];
    int thresh = std::stoi(std::string(argv[2]));

    // parse PB formula
    RdSat::SATSolver* solver = new RdSat::SATSolver();
    solver->parse_opb(filename, 0);
    vector<uint32_t>& sampling_set = solver->get_sampling_vars();

    // BSAT
    for (int i=0; i<thresh; i++) {
        vector<Lit> assumps;
        lbool ret = solver->solve(&assumps);

        //ban solution
        if (ret == CMSat::l_True) { 
            vector<lbool> model = solver->get_model();
            cout << "SAT" << endl;
            for (auto var: sampling_set)
                cout << Lit(var, solver->get_model()[var] == CMSat::l_False) << ' ';
            cout << '0' << endl;

            vector<Lit> lits;
            for (const uint32_t var: sampling_set) {
                assert(solver->get_model()[var] != CMSat::l_Undef);
                lits.push_back(Lit(var, solver->get_model()[var] == CMSat::l_True));
            }
            solver->add_clause(lits);
        } else {
            cout << "UNSAT" << endl;
            break;
        }
    }

    delete solver;

    return 0;
}
