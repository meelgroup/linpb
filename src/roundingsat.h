#include<iostream>
#include<vector>
#include "cryptominisat5/solvertypesmini.h"
#include "globals.hpp"
#include "ConstrExp.hpp"

#ifndef __ROUNDINGSAT_H__
#define __ROUNDINGSAT_H__

using std::string;
using std::vector;
using CMSat::Lit;
using CMSat::lbool;

namespace RdSat {
	class SATSolver {
		public:
			SATSolver();
			// ~SATSolver();
			void setup_sampling_set(const vector<uint32_t>& vars);
			int parse_opb(const string& filename, uint32_t exp_hash);
			uint32_t nVars() { return num_vars;}
			uint32_t new_var();
			lbool solve(vector<Lit>*);
			void add_clause(const vector<Lit>&);
			void add_xor_clause(const vector<uint32_t>&, bool);
			vector<lbool>& get_model() { return model; }
			vector<uint32_t>& get_sampling_vars() { return sampling_set; }
			void set_num_hash(uint32_t num);
			uint32_t get_num_hash() { return num_hash; }
			// simplify(vector<Lit>& assumps);
			int verb = 0;
			void set_verbosity(int);
			// print_stats();
			// get_text_version_info();
            		void set_config_gje(uint32_t);
            double SOLVETIME = 0, CATIME = 0, PROPTIME = 0, GAUSSTIME = 0;
            void print_time();
		private:
			rs::Solver* solver;
			// rs::CeArb* objective;
			bool add_clause_unsat = false;
			bool add_xor_unsat = false;
			uint32_t num_vars = 0;
			string opb_content;
			const string tempfile = "tempfile";
			vector<Lit> pre_assumps;
			vector<lbool> model;
			vector<uint32_t> sampling_set;
			vector<vector<uint64_t>> hashes;
			vector<uint64_t> exp_rhs;
			vector<uint32_t> aux_vars;
			uint32_t num_hash = 0;
			bool first_run = true;
			string exec(const char* cmd);
			string write_to_tempfile(vector<Lit>*);
			lbool parse_result(string);
            		uint32_t config_gje = 4;
	};
}

#endif // __ROUNDINGSAT_H__
