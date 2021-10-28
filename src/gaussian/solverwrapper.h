#ifndef GAUSSWRAPPER_H
#define GAUSSWRAPPER_H

#include "../Solver.hpp"
#include "solvertypesmini.h.in"
#include "xor.h"
#include <cstring>
#include "vardata.h"
#include "Vec.h"
#include "propby.h"
#include "gaussian.h"
#include "gqueuedata.h"
#include "gausswatched.h"
#include <limits>

namespace CMSat {

class EGaussian;
class TopLevelGaussAbst;

const uint32_t MAX_VAR = rs::CRef_Undef.ofs>>1;

struct ConfWrapper
{
	int verbosity;
	unsigned maxXORMatrix = 400ULL;
	int do_print_times = 1;
	ConfWrapper (int verb) : verbosity(verb) {}
	struct GaussConfWrapper
	{
		GaussConfWrapper() : min_usefulness_cutoff(0.2) {}
		double min_usefulness_cutoff;
	} gaussconf;
	std::string print_times(const double time_used) const;
};

struct AssignWrapper
{
	AssignWrapper(rs::Solver* solver=NULL): rsolver(solver) {}
	rs::Solver* rsolver;
	const lbool& operator[](const uint32_t idx) const;
};

struct Trail {
	Trail () {}
        Trail (Lit _lit, uint32_t _lev) : lit(_lit), lev(_lev) {}
	Lit lit;
	uint32_t lev;
};

struct TrailWrapper
{
	TrailWrapper(rs::Solver* solver): rsolver(solver) {}
	rs::Solver* rsolver;
	size_t size();
	Trail operator[](int idx);
};

struct VarDataWrapper
{
	VarDataWrapper(rs::Solver* solver): rsolver(solver) {}
	rs::Solver* rsolver;
	VarData operator[](int var);
};

struct ClauseCleaner
{
	Solver* solver;
	ClauseCleaner(Solver* solver) : solver(solver) {}
	bool clean_xor_clauses(vector<Xor>& xors);
	bool clean_one_xor(Xor& x);
};

class Solver
{
	public:
		Solver(rs::Solver* solver = NULL);
		~Solver();
		bool ok;
		bool okay();
		uint32_t decisionLevel();
		size_t nVars();
		lbool value(const Lit p);
		lbool value(const uint32_t x);
		ConfWrapper conf;
		vector<rs::ActValV>& var_act_vsids;
		// vector<uint16_t> seen;
		vec<vec<GaussWatched>> gwatches;
		uint32_t gqhead;
	        vector<EGaussian*> gmatrices;
		vector<GaussQData> gqueuedata;
		vec<PropBy> greasons;
		rs::CeSuper propagate();
		void enqueue(const Lit p);
		void enqueue(const Lit p, const uint32_t level, const PropBy from);
		rs::Ce32 get_reason(rs::CRef ref);
		rs::Ce32 PropBy2Contr(PropBy& reason);
		AssignWrapper assigns;
		TrailWrapper trail;
		VarDataWrapper varData;
		bool init_all_matrices();
		void enlarge_minimal_datastructs(long long n);
		void clear_gauss_matrices();
		bool init_gauss();
		vector<Xor> xors;
		bool add_xor_clause(const std::vector<unsigned>& vars, bool rhs);
		bool add_xor_clause_inter(const std::vector<Lit>& lits, bool rhs, const bool attach=true);
		bool parse_xor(std::string line);
		static Lit rsLit2cmsLit(rs::Lit l);
		static rs::Lit cmsLit2rsLit(Lit l);
		void print_matrix_stats();
		template<class T> void clean_xor_vars_no_prop(T& ps, bool& rhs);
		bool fully_enqueue_this(const Lit lit);
		ClauseCleaner* clauseCleaner;
		size_t trail_size() { return trail.size(); }
		void cancel();
		enum class gauss_ret {g_cont, g_nothing, g_false};
		gauss_ret gauss_jordan_elim();
		gauss_ret run_gauss(rs::CeSuper& confl);
		PropBy confl_propby;
		void clean_xor_no_prop(std::vector<Lit>& ps, bool& rhs);
		void check_need_gauss_jordan_disable();
		bool add_clause(const std::vector<Lit>& clause);
		bool add_every_combination_xor(std::vector<Lit> lits);
		TopLevelGaussAbst *topLevelGauss;
		void run_top_level_gauss();
	private:
		rs::Solver* rsolver;
    		bool xor_clauses_updated = false;
		// uint32_t qhead = 0;
};

// template member function has to be defined in the header
template<class T>
void Solver::clean_xor_vars_no_prop(T& ps, bool& rhs)
{
    std::sort(ps.begin(), ps.end());
    uint32_t p;
    uint32_t i, j;
    for (i = j = 0, p = std::numeric_limits<uint32_t>::max(); i != ps.size(); i++) {
        if (ps[i] == p) {
            //added, but easily removed
            j--;
            p = std::numeric_limits<uint32_t>::max();

            //Flip rhs if neccessary
            if (value(ps[i]) != l_Undef) {
                rhs ^= value(ps[i]) == l_True;
            }

        } else if (value(ps[i]) == l_Undef) {
            //Add and remember as last one to have been added
            ps[j++] = p = ps[i];

            assert(varData[p].removed != Removed::elimed);
        } else {
            //modify rhs instead of adding
            rhs ^= value(ps[i]) == l_True;
        }
    }
    ps.resize(ps.size() - (i - j));
}

} // namespace CMSat

#endif	//GAUSSWRAPPER_H 
