#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include "roundingsat.h"
#include <fstream>
#include <sstream>
#include <limits.h>
#include <algorithm>
#include "globals.hpp"
#include "run.hpp"
#include "parsing.hpp"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ofstream;
using CMSat::Lit;
using CMSat::lbool;


namespace RdSat {
SATSolver::SATSolver()
{
	solver = &rs::run::solver;
	solver->init();
	// objective = rs::run::solver.cePools.takeArb();

	rs::options.verbosity.parse("0");
	rs::options.printSol.parse("1");
	rs::options.dbCleanInc.parse("100");
	if (config_gje < 4) {
		rs::options.configGJE.parse(std::to_string(config_gje));
	}
	// cout << "verbosity: " << rs::options.verbosity.get() << endl;
	// cout << "printSol: " << bool(rs::options.printSol) << endl;
	// cout << "dbCleanInc: " << rs::options.dbCleanInc.get() << endl;

	// solver->initLP(objective);	
	// rs::run::solver.initLP(*objective);
}


void SATSolver::setup_sampling_set(const vector<uint32_t>& vars){
	sampling_set = vars;
}


int SATSolver::parse_opb(const string& filename, uint32_t exp_hash){
	ifstream file(filename);
    	if (!file) { cout << "Could not open " << filename << endl; exit(-1); }
	string header, tmp;
	std::getline(file, header);
	std::istringstream ss(header);
	ss >> tmp >> tmp >> num_vars; // skip * and #variable=	

	sampling_set.clear();
	while(file.peek() == '*') {
		std::getline(file, tmp);
		if (tmp.find("ind") !=  string::npos){
			std::stringstream ss_ind(tmp);
			uint32_t var;
			ss_ind.ignore(5); // ignore "* ind"
			while(ss_ind >> var && var){
				assert(var>0 && var<=num_vars);
				sampling_set.push_back(var-1);
			}

		}
	}
	if (sampling_set.empty()){
		for (uint32_t i=0; i<num_vars; i++){
			sampling_set.push_back(i);
		}
	}

	// sum hash
	if (exp_hash) {
		// auxiliary vars for hash
		for (uint32_t i=1; i<sampling_set.size(); i*=2)
			aux_vars.push_back(new_var());
	}

	file.clear();
	file.seekg(0);
	std::getline(file, tmp);	// skip header
	opb_content.assign(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>() );
	
	file.close();

	return 1;
}

string SATSolver::write_to_tempfile(vector<Lit>* p_assumps) {
	std::ostringstream file;

	// header info
	file << "* #variable= " << num_vars << endl;	
	// content
	file << opb_content;
	// assumption
	for (Lit lit: *p_assumps){
		file << "1 ";
		if (lit == CMSat::lit_Undef) {
			file << "lit_Undef";
		} else {
			file << (lit.sign() ? "~x" : "x") << (lit.var() + 1);
		}
		file << " >= 1 ;" << endl;	
	}	
	// sum hashes (another 2-universal hash)
	if (num_hash > 0) {
		for (uint32_t i=0; i<sampling_set.size(); i++) {
			file << hashes[num_hash-1][i] << " x" << sampling_set[i]+1 << ' ';
		}
		// modulo ........
		uint64_t base = 1;
		for (uint32_t i=0; i<num_hash; i++)
			base *= 2;
		assert(sampling_set.size() < ULLONG_MAX/base);
		int idx = 0;
		for (uint64_t i = 1; i < sampling_set.size(); i *= 2, idx++){
			file << "-" << i*base << " x" << aux_vars[idx] << ' ';
		}
		file << "= " << exp_rhs[num_hash-1] << " ; " << endl;
	}


	// for debug, output temporary opb file		
	// ofstream outfile(tempfile);
	// outfile << file.str();
	// outfile.close();
	

	return file.str();
}


uint32_t SATSolver::new_var(){
	++num_vars;
	solver->setNbVars(num_vars);
	return num_vars;
}

lbool SATSolver::solve(vector<Lit>* p_assumps){
	/*
	assumps.push_back(Lit(0, false));
	assumps.push_back(Lit(1, true));

	vector<Lit> clause;
	clause.push_back(Lit(2, false));
	clause.push_back(Lit(3, true));
	add_clause(clause);
	
	vector<uint32_t> vars;
	vars.push_back(0);
	vars.push_back(1);
	vars.push_back(2);
	vars.push_back(3);
	vars.push_back(4);
	vars.push_back(5);
	vars.push_back(6);
	vars.push_back(7);
	vars.push_back(8);
	vars.push_back(10);
	add_xor_clause(vars, false);
	*/

	#ifdef TEST
	vector<uint32_t> vars;
	vars.push_back(0);
	vars.push_back(1);
	vars.push_back(2);
	add_xor_clause(vars, false);
	#endif
	
	// reset solver
	// 1. only reset when assumps change.
	// 2. how to reset solver 
	// (a.) destruct and construct, which needs to modify copy constructor 
	// (b.) reset attributes which needs to handle ca, constraints and any others relevant. 
	// using strategy (b) now.
	// 3. how to record original constraints, added xor and clauses 
	// (a.) internal string which needs parsing 
	// (b.) internal xor and CeSuper. 
	// using strategy (a) now.
	
	// solver.~Solver();
	// *solver = rs::Solver(); 
	
	// debug write intermediate variables
	// write_to_tempfile(p_assumps);
	if (first_run or pre_assumps != *p_assumps) {
		string opb = write_to_tempfile(p_assumps);

		rs::stats = rs::Stats();

		solver->reset();
		solver->init();
		solver->setNbVars(num_vars);
		rs::CeArb objective = rs::run::solver.cePools.takeArb();

		std::istringstream fin(opb);
		// std::ifstream fin(tempfile);
		
		// corner case, finding UNSAT when reading the file
		std::stringstream ss;
		// change stdout to ss, save stdout
		auto old_buf = std::cout.rdbuf(ss.rdbuf());

		rs::parsing::file_read(fin, *solver, objective);	

		// reset cout to stdout
		std::cout.rdbuf(old_buf);
		lbool ret = parse_result(ss.str());
		if (ret == CMSat::l_False)
			return ret;

		solver->initLP(objective);	

		pre_assumps = *p_assumps;
		first_run = false;
		add_clause_unsat = false;
		add_xor_unsat = false;
	}
	else {
		solver->lastSol = {0};		// clear the last solution to avoid confusion with optim problem.

		if (add_clause_unsat || add_xor_unsat) {
			return CMSat::l_False;
		}

		solver->nconfl_to_reduce = 2000;
		solver->nconfl_to_restart = 0;
		rs::stats.NRESTARTS = 0;
		rs::stats.NCONFL = 0;
		rs::stats.NCLEANUP = 0;
        
        rs::stats.SOLVETIME = 0;
        rs::stats.CATIME = 0;
        rs::stats.PROPTIME = 0;
        rs::stats.GAUSSTIME = 0;
		// rs::stats.NLPPIVOTSROOT = 0;
		// rs::stats.NLPPIVOTSINTERNAL = 0;

		/*
		// local search near previous solution
		std::stringstream ss;
		auto old_buf = std::cout.rdbuf(ss.rdbuf());
		
		rs::CeArb objective = rs::run::solver.cePools.takeArb();
		int bound_inprocess = 1;
		rs::run::run(objective, bound_inprocess);
		std::cout.rdbuf(old_buf);
	
		if (bound_inprocess) {
			lbool ret = parse_result(ss.str());
			return ret;
		}

		cout << "reset" << endl;
		*/
		// reset search path
		// solver->v_vsids_inc = 1.0;
		// solver->c_vsids_inc = 1.0;
  		// solver->activity.resize(num_vars + 1, 1 / rs::actLimitV);
		// // solver->activity.clear();
		// solver->order_heap.recalculate();
		// solver->reset_basis();
		// solver->reset_firstRun();
		// solver->reset_qhead();		// recheck starting assignment for a new running, because constraints may change
	}
	// cout << "nconfl_to_reduce: " << solver->nconfl_to_reduce << endl;
	// cout << "nconfl_to_restart: " << solver->nconfl_to_restart << endl;
	// cout << "solver->v_vsids_inc: " << solver->v_vsids_inc << endl;
	// cout << "solver->c_vsids_inc: " << solver->c_vsids_inc << endl;

	std::stringstream ss;

	// change stdout to ss, save stdout
	auto old_buf = std::cout.rdbuf(ss.rdbuf());

	/*
	// previous strategy, not safe, because constraints may be optimized by assumptions.
	// add assumptions
	uint32_t pre_size = solver->getTrail().size();
	for (auto l : *p_assumps) {
		rs::Lit lit = CMSat::Solver::cmsLit2rsLit(l);
		if (rs::isUnknown(solver->getPos(), lit)) {
			solver->uncheckedEnqueue(lit, rs::CRef_Undef);
		}
	}
	solver->lastSol = {0};		// clear the last solution to avoid confusion with optim problem.
	solver->runPropagation(true);
	*/

	try {
		rs::CeArb objective = rs::run::solver.cePools.takeArb();
		// int bound_inprocess = 0;
		// rs::run::run(objective, bound_inprocess);
		rs::run::run(objective);
		
		/*
		// previous strategy
		// remove aasumptions
		int num_undo = solver->getTrail().size()-pre_size;
		for (int i=0; i<num_undo; i++) {
			solver->undoOne();
		}
		*/

		// reset cout to stdout
		std::cout.rdbuf(old_buf);
	} catch(...) {
		std::cout.rdbuf(old_buf);
		cout << "Abnormal Exit!" << endl;
		cout << ss.str() << endl;
	}
    
    if (verb) 
        cout << ss.str() << endl;	
	lbool ret = parse_result(ss.str());

    // time statistics
    // cout << "Solve time: " << rs::stats.SOLVETIME << endl
    //      << "Conflict analysis time: " << rs::stats.CATIME << endl
    //      << "Propagation time: " << rs::stats.PROPTIME << endl
    //      << "Gauss time: " << rs::stats.GAUSSTIME << endl;
    SOLVETIME += rs::stats.SOLVETIME;
    CATIME += rs::stats.CATIME;
    PROPTIME += rs::stats.PROPTIME;
    GAUSSTIME += rs::stats.GAUSSTIME;

	// old API
	// string res = exec("roundingsat tempfile --verbosity=0 --print-sol=1");
	// lbool ret = parse_result(res);

	// cout << ret << endl;
	return ret;
}

lbool SATSolver::parse_result(string res)
{
	string tmp;
	std::stringstream ss(res);
	while (ss >> tmp && tmp != "UNSATISFIABLE" && tmp != "SATISFIABLE");	
	if (tmp == "UNSATISFIABLE") {
		return lbool(CMSat::l_False);
	}
	else if (tmp != "SATISFIABLE") {
		return lbool(CMSat::l_Undef);
	}
	ss.ignore(2); // ignore LF and v
	model.clear();
	while(ss >> tmp) {
		assert(tmp.length() > 1);
		if (tmp[0] == '-')
			model.push_back(lbool(CMSat::l_False));
		else
			model.push_back(lbool(CMSat::l_True));
	}
	for (uint32_t i=model.size(); i<num_vars; i++){
		model.push_back(lbool(CMSat::l_False));
	}
	assert(model.size() == num_vars);
	
	return lbool(CMSat::l_True);
}

void SATSolver::add_clause(const vector<Lit>& clause)
{
	rs::Ce32 constr = solver->cePools.take32();
	for (auto& lit : clause) {
		assert(lit.var() < num_vars);
		constr->addLhs(1, CMSat::Solver::cmsLit2rsLit(lit));
	}
	constr->addRhs(1);
	// in case return UNSAT when adding a new constraint
	if (solver->addConstraint(constr, rs::Origin::FORMULA).second == rs::ID_Unsat) {
		add_clause_unsat = true;
	}
	
	// ##### old API #####
	std::stringstream ss;

	for (Lit lit: clause){
		ss << "1 ";
		if (lit == CMSat::lit_Undef) {
			ss << "lit_Undef ";
		} else {
			ss << (lit.sign() ? "~x" : "x") << (lit.var() + 1) << ' ';
		}
	}	
	ss << ">= 1 ;" << endl;		

	opb_content += ss.str();
}

void SATSolver::add_xor_clause(const vector<uint32_t>& vars, bool rhs)
{
	if (!solver->getGauss()->add_xor_clause(vars, rhs)) 
		add_xor_unsat = true;

	// ##### old API #####
	std::stringstream ss;
	
	// xor 
	if (config_gje > 0) {
		// GJE handle xor
		ss << "* xor ";
		for (auto var: vars){
			ss << "x" << var+1 << " ";
		}
		ss << rhs << endl;	
	} else {
		// xor -> pb
		// normal encoding
		for (auto var: vars){
		    ss << "1 x" << var+1 << " ";
		}
		assert(vars.size() < INT_MAX/2);
		for (uint32_t i=2; i <= vars.size(); i *= 2){
		    ss << "-" << i << " x" << new_var() << ' ';
		}
		ss << "= " << rhs << " ;" << endl;
	}

	/*
	// elaborate xor (hierarchical encoding)
	assert(vars.size() > 1);
	uint32_t ctrl_var = vars.back();
	vector<uint32_t> cur, next;
	for (uint32_t i=0; i<vars.size()-1; i++)
		cur.push_back(vars[i]);
	// std::random_shuffle(cur.begin(), cur.end());
	while(cur.size() > 1) {
		for (uint32_t i=0; i<cur.size(); i+=2) {
			if (i == cur.size() - 1) {
				next.push_back(cur.back());
			}
			else {
				next.push_back(new_var()-1);
				ss << "1 x" << cur[i]+1 << " ";
				ss << "1 x" << cur[i+1]+1 << " ";
				ss << "-1 x" << next.back()+1 << " ";
				ss << "-2 x" << new_var() << " ";
				ss << "= 0 ; " << endl;
			}
		}
		cur = next;
		next.clear();
	}
	assert(cur.size() == 1);
	if (rhs) {
		ss << "1 x" << ctrl_var+1 << " ";
		ss << "1 x" << cur[0]+1 << " ";
		ss << "= 1 ; " << endl;
	}
	else {
		ss << "1 x" << ctrl_var+1 << " ";
		ss << "-1 x" << cur[0]+1 << " ";
		ss << "= 0 ; " << endl;
	}
	*/

	opb_content += ss.str();
}

double rand1()
{
	return (float) rand()/RAND_MAX;
}

// only used by sum hash
void SATSolver::set_num_hash(uint32_t num)
{
	if (hashes.empty()) {
		vector<uint64_t> hash;
		for (uint32_t i=0; i<sampling_set.size(); i++)
			hash.push_back(rand1()<0.5);
		hashes.push_back(hash);
		exp_rhs.push_back(rand1()<0.5);
	}

	if (num > hashes.size()) {
		uint64_t base = 1;
		for (uint32_t i=0; i<hashes.size(); i++)
			base *= 2;
		for (uint32_t pre=hashes.size(); pre<num; pre++, base*=2) {
			vector<uint64_t> hash;
			for (uint32_t i=0; i<sampling_set.size(); i++)
				hash.push_back(hashes[pre-1][i] + base*(rand1()<0.5));
			hashes.push_back(hash);
			exp_rhs.push_back(exp_rhs[pre-1] + base*(rand1()<0.5));
		}	
	}
	/*
	vector<uint64_t>& hash = hashes.back();
	for (auto e: hash)
		cout << e << ' ';
	cout << endl;
	cout << "rhs: " << rhs.back() << endl;
	*/

	num_hash = num;
}

string SATSolver::exec(const char* cmd) {
    	char buffer[128];
	std::string result = "";
    	FILE* pipe = popen(cmd, "r");
	if (!pipe) throw std::runtime_error("popen() failed!");
	try {
		while (fgets(buffer, sizeof buffer, pipe) != NULL) {
			result += buffer;
		}
	} catch (...) {
		pclose(pipe);
		throw;
	}
	pclose(pipe);
	return result;

}

/*
// For test
#ifdef TEST
int main(int argc, char** argv)
{
	if (argc < 2) {
		cout << "Error: no input file!" << endl;
		return -1;
	}

	SATSolver* solver = new SATSolver();

	solver->parse_opb(string(argv[1]), 0);
	// solver->set_num_hash(10);
	solver->solve(new vector<Lit>());

	return 0;
}
#endif
*/

void SATSolver::set_verbosity(int verbosity) {
    verb = verbosity;
	rs::options.verbosity.parse(std::to_string(verb));
}

void SATSolver::set_config_gje(uint32_t configGJE) {
    config_gje = configGJE;
    rs::options.configGJE.parse(std::to_string(config_gje));
}


void SATSolver::print_time() {
    cout << "Solve time: " << SOLVETIME << endl
         << "Conflict analysis time: " << CATIME << endl
         << "Propagation time: " << PROPTIME << endl
         << "Gauss time: " << GAUSSTIME << endl;
    
}

} // namespace RdSat
