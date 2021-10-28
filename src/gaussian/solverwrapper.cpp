#include "solverwrapper.h"
#include "Options.hpp"
#include <sstream>
#include "toplevelgaussabst.h"

#ifdef USE_M4RI
#include "toplevelgauss.h"
#endif

using namespace CMSat;

DLL_PUBLIC std::string ConfWrapper::print_times(const double time_used) const 
{
	if (do_print_times) {
		std::stringstream ss;
		ss
		<< " T: " << std::setprecision(2) << std::fixed << time_used;

		return ss.str();
	}

	return std::string();
}

const lbool& AssignWrapper::operator[](const uint32_t idx) const 
{
	if (rs::isUnknown(rsolver->getPos(), idx+1)) {
		return l_Undef; 
	}
	else if (rs::isTrue(rsolver->getLevel(), idx+1))
		return l_True;
	else
		return l_False;
}

size_t TrailWrapper::size() 
{
	assert(rsolver != NULL);	
	return rsolver->getTrail().size(); 
}


Trail TrailWrapper::operator[](int idx) 
{
	rs::Lit rslit = rsolver->getTrail()[idx];
	assert(!rs::isUnknown(rsolver->getPos(), rslit));
	return Trail(Solver::rsLit2cmsLit(rslit), rsolver->getLevel()[rslit]);
}

VarData VarDataWrapper::operator[](int var) {
	VarData data;
	if(rs::isTrue(rsolver->getLevel(), var+1)){
		data.level = rsolver->getLevel()[var+1];
	}
	else if(rs::isFalse(rsolver->getLevel(), var+1)){
		data.level = rsolver->getLevel()[-(var+1)];
	} else {
		assert(rs::isUnknown(rsolver->getPos(), var+1));
		data.level = 0;
	}

	return data;
}

// need modification here, add wrapper to var_act_vsids, then handle the case when solver == NULL
Solver::Solver(rs::Solver* solver):rsolver(solver), var_act_vsids(solver->activity),
	assigns(solver), trail(solver), varData(solver), conf(rs::options.verbosity.get())
{
	assert(solver != NULL);
	ok = true;
	clauseCleaner = new ClauseCleaner(this);
	xors.clear();

	topLevelGauss = new TopLevelGaussAbst;
	#ifdef USE_M4RI
	delete topLevelGauss;
	topLevelGauss = new TopLevelGauss(this);
	#endif
}

Solver::~Solver()
{
	delete clauseCleaner;
	clear_gauss_matrices();
	delete topLevelGauss;
}

bool Solver::okay()
{
	return ok;
}	

uint32_t Solver::decisionLevel()
{
	return rsolver->decisionLevel();
}

size_t Solver::nVars()
{
	return rsolver->getNbVars();
}

lbool Solver::value(const Lit p)
{
	return assigns[p.var()] ^ p.sign();
}

lbool Solver::value (const uint32_t x)
{
	return assigns[x];
}

rs::CeSuper Solver::propagate()
{
	return rsolver->runPropagation(true);
}

void Solver::enqueue(const Lit p)
{
	rs::Lit l = p.sign() ? -(p.var()+1) : p.var()+1;
	rsolver->uncheckedEnqueue(l, rs::CRef_Undef);
}

void Solver::enqueue(const Lit p, const uint32_t level, const PropBy from)
{
	assert(p.var() < greasons.size());
	greasons[p.var()] = from;
	assert(rs::CRef_Undef.ofs > p.var() + 1 + MAX_VAR);
	rs::CRef reason = {rs::CRef_Undef.ofs-p.var()-1};	// -1 is to avoid CRef_Undef when p.var() == 0
	rsolver->propagate(Solver::cmsLit2rsLit(p), reason);
}

rs::Ce32 Solver::get_reason(rs::CRef ref)
{
	assert(ref.ofs > MAX_VAR);
	int var = rs::CRef_Undef.ofs - ref.ofs - 1;
	return PropBy2Contr(greasons[var]);
}

rs::Ce32 Solver::PropBy2Contr(PropBy& reason)
{
	vector<Lit>* xor_reason = gmatrices[reason.get_matrix_num()] -> get_reason(reason.get_row_num());	

	rs::Ce32 constr = rsolver->cePools.take32();
	for (auto& lit : *xor_reason) {
		constr->addLhs(1, Solver::cmsLit2rsLit(lit));
	}
	constr->addRhs(1);
    constr->orig = rs::Origin::GAUSS;

	return constr;
}

Lit Solver::rsLit2cmsLit(rs::Lit l)
{
	return Lit(rs::toVar(l)-1, l<0);
}

DLL_PUBLIC rs::Lit Solver::cmsLit2rsLit(Lit l)
{
	return l.sign() ? -(l.var()+1) : l.var()+1;
}

bool Solver::init_all_matrices()
{
    assert(ok);
    assert(decisionLevel() == 0);

    assert(gmatrices.size() == gqueuedata.size());
    for (uint32_t i = 0; i < gmatrices.size(); i++) {
        auto& g = gmatrices[i];
        bool created = false;
        //initial arrary. return true is fine,
        //return false means solver already false;
        if (!g->full_init(created)) {
            return false;
        }
        if (!ok) {
            break;
        }
        if (!created) {
            gqueuedata[i].engaus_disable = true;
            delete g;
            if (conf.verbosity > 5) {
                cout << "DELETED matrix" << endl;
            }
            g = NULL;
        }
    }

    uint32_t j = 0;
    bool modified = false;
    for (uint32_t i = 0; i < gqueuedata.size(); i++) {
        if (gmatrices[i] != NULL) {
            gmatrices[j] = gmatrices[i];
            gmatrices[j]->update_matrix_no(j);
            gqueuedata[j] = gqueuedata[i];

            if (modified) {
                for (size_t var = 0; var < nVars(); var++) {
                    for(GaussWatched* k = gwatches[var].begin();
                        k != gwatches[var].end();
                        k++)
                    {
                        if (k->matrix_num == i) {
                            k->matrix_num = j;
                        }
                    }
                }
            }
            j++;
        } else {
            modified = true;
        }
    }
    gqueuedata.resize(j);
    gmatrices.resize(j);

    return okay();
}

void Solver::enlarge_minimal_datastructs(long long n)
{
	gwatches.insert(2*n);
	greasons.insert(n);
}	

void Solver::clear_gauss_matrices()
{
    xor_clauses_updated = true;
    for(uint32_t i = 0; i < gqueuedata.size(); i++) {
        auto gqd = gqueuedata[i];
        if (conf.verbosity >= 2) {
            cout
            << "c [mat" << i << "] num_props       : "
            << print_value_kilo_mega(gqd.num_props) << endl
            << "c [mat" << i << "] num_conflicts   : "
            << print_value_kilo_mega(gqd.num_conflicts)  << endl;
        }
    }

    if (conf.verbosity >= 1) {
        print_matrix_stats();
    }
    for(EGaussian* g: gmatrices) {
        delete g;
    }
    for(auto& w: gwatches) {
        w.clear();
    }
    gmatrices.clear();
    gqueuedata.clear();
}
		
bool Solver::init_gauss()
{
	clear_gauss_matrices();
	gqhead = trail.size();
	
	// init xors

	// init gmatrices and gqueuedata
	gmatrices.push_back(new EGaussian(this, 0, xors));
	gqueuedata.resize(gmatrices.size());
	ok = true;

	// init matrices
	if (!init_all_matrices()) {
		return false;
	}

	return true;
}

bool Solver::add_clause(const std::vector<Lit>& clause)
{
	rs::Ce32 constr = rsolver->cePools.take32();
	for (auto& lit : clause) {
		assert(lit.var() < rsolver->getNbVars());
		constr->addLhs(1, Solver::cmsLit2rsLit(lit));
	}
	constr->addRhs(1);
	// in case return UNSAT when adding a new constraint
	if (rsolver->addConstraint(constr, rs::Origin::FORMULA).second == rs::ID_Unsat) {
		ok = false;
	}

	return ok;
}

DLL_PUBLIC bool Solver::add_xor_clause(const std::vector<unsigned>& vars, bool rhs)
{
	if (!ok) {
		return false;
	}

	vector<Lit> lits(vars.size());
	for(size_t i = 0; i < vars.size(); i++) {
		assert(vars[i] < rsolver->getNbVars());
		lits[i] = Lit(vars[i], false);
	}
	
	add_xor_clause_inter(lits, rhs);
	return ok;
}

bool Solver::add_xor_clause_inter(const std::vector<Lit>& lits, bool rhs, const bool attach)
{
	assert(ok);
	// assert(rsolver->qhead == trail.size());
	assert(decisionLevel() == 0);

	vector<Lit> ps(lits);
	for(Lit& lit: ps) {
		if (lit.sign()) {
			rhs ^= true;
			lit ^= true;
		}
	}
	clean_xor_no_prop(ps, rhs);

	if (ps.empty()) {
		if (rhs) {
			ok = false;
		}
		return ok;
	}

	if (ps.size() > 2) {
		xor_clauses_updated = true;
	}
	ps[0] ^= rhs;

	if (ps.size() > 2) {
		xors.push_back(Xor(ps, rhs, vector<uint32_t>()));
	} else {
		add_every_combination_xor(ps);
	}

	// debug
	if (conf.verbosity > 0)
		cout << xors.back() << endl;

	return ok;
}

bool Solver::add_every_combination_xor(std::vector<Lit> lits)
{
	// only for 1-xor and 2-xor now
	assert(lits.size() > 0 && lits.size() < 3);

	if (lits.size() == 1) {
		lits[0] ^= 1;
		add_clause(lits);
	} else { // size == 2
		lits[0] ^= 1;
		add_clause(lits);
		lits[0] ^= 1; lits[1] ^= 1;
		add_clause(lits);
	}

	return ok;
}

void Solver::clean_xor_no_prop(std::vector<Lit>& ps, bool& rhs)
{
    std::sort(ps.begin(), ps.end());
    Lit p;
    uint32_t i, j;
    for (i = j = 0, p = lit_Undef; i != ps.size(); i++) {
        assert(ps[i].sign() == false);

        if (ps[i].var() == p.var()) {
            //added, but easily removed
            j--;
            p = lit_Undef;

            //Flip rhs if neccessary
            if (value(ps[i]) != l_Undef) {
                rhs ^= value(ps[i]) == l_True;
            }

        } else if (value(ps[i]) == l_Undef) {
            //Add and remember as last one to have been added
            ps[j++] = p = ps[i];

            assert(varData[p.var()].removed != Removed::elimed);
        } else {
            //modify rhs instead of adding
            rhs ^= value(ps[i]) == l_True;
        }
    }
    ps.resize(ps.size() - (i - j));
}

bool Solver::parse_xor(std::string line)
{
	if (line.substr(0, 5) == "* xor") {
		std::stringstream ss(line);
		string var;
		ss >> var >> var; // skip * and xor
		std::vector<unsigned> vars;
		bool rhs;
		while(ss >> var) {
			if (var[0] == 'x') {
				assert(std::to_string(std::stoi(var.substr(1))) == var.substr(1));
				uint32_t rs_var = std::stoi(var.substr(1));
				assert(rs_var <= rsolver->getNbVars() && "Variable Overflow");
				vars.push_back(rs_var - 1);	// rs::var_id - 1 === CMSat::var_id
			}
			else {
				assert(var.length() == 1);
				assert(var == "1" || var == "0");
				rhs = (var == "1");
				return add_xor_clause(vars, rhs);
			}
		}
	}
	
	assert(false && "Parsing XOR Error");
	return false;
}

void Solver::print_matrix_stats()
{
	for(EGaussian* g: gmatrices) {
		if (g) {
			g->print_matrix_stats(conf.verbosity);
		}
	}
}

bool Solver::fully_enqueue_this(const Lit lit)
{
    const lbool val = value(lit);
    if (val == l_Undef) {
        assert(varData[lit.var()].removed == Removed::none);
        enqueue(lit);
        // ok = propagate<true>().isNULL();
	ok = !propagate();

        if (!ok) {
            return false;
        }
    } else if (val == l_False) {
        ok = false;
        return false;
    }
    return true;
}

bool ClauseCleaner::clean_xor_clauses(vector<Xor>& xors)
{
    assert(solver->ok);
    #ifdef VERBOSE_DEBUG
    for(Xor& x : xors) {
        cout << "orig XOR: " << x << endl;
    }
    #endif

    size_t last_trail = std::numeric_limits<size_t>::max();
    while(last_trail != solver->trail_size()) {
        last_trail = solver->trail_size();
        size_t i = 0;
        size_t j = 0;
        for(size_t size = xors.size(); i < size; i++) {
            Xor& x = xors[i];
            //cout << "Checking to keep xor: " << x << endl;
            const bool keep = clean_one_xor(x);
            if (!solver->ok) {
                return false;
            }

            if (keep) {
                xors[j++] = x;
            } else {
		/*
                solver->removed_xorclauses_clash_vars.insert(
                    solver->removed_xorclauses_clash_vars.end()
                    , x.clash_vars.begin()
                    , x.clash_vars.end()
                );*/
                //cout << "NOT keeping XOR" << endl;
            }
        }
        xors.resize(j);

        #ifdef VERBOSE_DEBUG
        for(Xor& x : xors) {
            cout << "cleaned XOR: " << x << endl;
        }
        #endif
    }
    return solver->okay();
}

bool ClauseCleaner::clean_one_xor(Xor& x)
{
    bool rhs = x.rhs;
    size_t i = 0;
    size_t j = 0;
    for(size_t size = x.size(); i < size; i++) {
        uint32_t var = x[i];
        if (solver->value(var) != l_Undef) {
            rhs ^= solver->value(var) == l_True;
        } else {
            x[j++] = var;
        }
    }
    x.resize(j);
    x.rhs = rhs;

    switch(x.size()) {
        case 0:
            solver->ok &= !x.rhs;
            return false;

        case 1: {
            solver->fully_enqueue_this(Lit(x[0], !x.rhs));
            return false;
        }
        case 2: {
            solver->add_xor_clause_inter(vars_to_lits(x), x.rhs, true);
            return false;
        }
        default: {
            return true;
        }
    }
}

void Solver::cancel()
{
        // if (!all_matrices_disabled) {
            for (uint32_t i = 0; i < gmatrices.size(); i++) {
                if (gmatrices[i] && !gqueuedata[i].engaus_disable) {
                    //cout << "->Gauss canceling" << endl;
                    gmatrices[i]->canceling();
                } else {
                    //cout << "->Gauss NULL" << endl;
                }
            }
        // }
	
}

Solver::gauss_ret Solver::gauss_jordan_elim()
{
    #ifdef VERBOSE_DEBUG
    cout << "Gauss searcher::Gauss_elimination called, declevel: " << decisionLevel() << endl;
    #endif

    for(uint32_t i = 0; i < gqueuedata.size(); i++) {
        if (gqueuedata[i].engaus_disable) {
            continue;
        }
        gqueuedata[i].reset();
        gmatrices[i]->update_cols_vals_set();
    }

    bool confl_in_gauss = false;
    while (gqhead <  trail.size()
        && !confl_in_gauss
    ) {
        const Lit p = trail[gqhead].lit;
        uint32_t currLevel = trail[gqhead].lev;
        gqhead++;

        assert(gwatches.size() > p.var());
        vec<GaussWatched>& ws = gwatches[p.var()];
        GaussWatched* i = ws.begin();
        GaussWatched* j = i;
        const GaussWatched* end = ws.end();
        #ifdef VERBOSE_DEBUG
        cout << "New GQHEAD: " << p << endl;
        #endif

        for (; i != end; i++) {
            if (gqueuedata[i->matrix_num].engaus_disable) {
                //remove watch and continue
                continue;
            }

            gqueuedata[i->matrix_num].new_resp_var = std::numeric_limits<uint32_t>::max();
            gqueuedata[i->matrix_num].new_resp_row = std::numeric_limits<uint32_t>::max();
            gqueuedata[i->matrix_num].do_eliminate = false;
            gqueuedata[i->matrix_num].currLevel = currLevel;

            if (gmatrices[i->matrix_num]->find_truths(
                i, j, p.var(), i->row_n, gqueuedata[i->matrix_num])
            ) {
                continue;
            } else {
                confl_in_gauss = true;
                i++;
                break;
            }
        }

        for (; i != end; i++) {
            *j++ = *i;
        }
        ws.shrink(i-j);

        for (size_t g = 0; g < gqueuedata.size(); g++) {
            if (gqueuedata[g].engaus_disable)
                continue;

            if (gqueuedata[g].do_eliminate) {
                gmatrices[g]->eliminate_col(p.var(), gqueuedata[g]);
                confl_in_gauss |= (gqueuedata[g].ret == gauss_res::confl);
            }
        }
	// propagate only a literal per XOR propagation if using shared watches
	if (rs::options.configGJE.get() == 3) {
        	break;
	}
    }

    #ifdef SLOW_DEBUG
    if (!confl_in_gauss) {
        for (size_t g = 0; g < gqueuedata.size(); g++) {
            if (gqueuedata[g].engaus_disable)
                continue;

            assert(solver->gqhead == solver->trail.size());
            gmatrices[g]->check_invariants();
        }
    }
    #endif

    gauss_ret finret = gauss_ret::g_nothing;
    for (GaussQData& gqd: gqueuedata) {
        if (gqd.engaus_disable)
            continue;

        //There was a conflict but this is not that matrix.
        //Just skip.
        if (confl_in_gauss && gqd.ret != gauss_res::confl) {
            continue;
        }

        switch (gqd.ret) {
            case gauss_res::confl :{
                gqd.num_conflicts++;
                // gqhead = rsolver->qhead = trail.size();
                // gqhead = rsolver->qhead;
                // bool ret = handle_conflict(gqd.confl);
                // if (!ret) return gauss_ret::g_false;
		confl_propby = gqd.confl;		// record conflict
                return gauss_ret::g_cont;
            }

            case gauss_res::prop:
                gqd.num_props++;
                finret = gauss_ret::g_cont;

            case gauss_res::none:
                //nothing
                break;

            default:
                assert(false);
                return gauss_ret::g_nothing;
        }
    }
    #ifdef VERBOSE_DEBUG
    cout << "Exiting GJ" << endl;
    #endif
    return finret;
}

Solver::gauss_ret Solver::run_gauss(rs::CeSuper& confl)
{
	confl_propby = PropBy();
	gauss_ret ret = gauss_jordan_elim();
	if (confl_propby.isNULL())
		confl = rs::CeNull();
	else
		confl = PropBy2Contr(confl_propby);

	return ret;
}

void Solver::check_need_gauss_jordan_disable()
{
	for(uint32_t i = 0; i < gqueuedata.size(); i++) {
		auto& gqd = gqueuedata[i];

		gqd.reset();
		gmatrices[i]->update_cols_vals_set();
	}

	// assert(gqhead <= rsolver->qhead);
}

void Solver::run_top_level_gauss()
{
	#ifdef USE_M4RI
	if (topLevelGauss != NULL) {
	    auto tmp_xors = xors;
	    assert(okay());
	    // solver->ok = finder.xor_together_xors(xors);
	    if (ok) {
		vector<Lit> out_changed_occur;
		// finder.remove_xors_without_connecting_vars(xors);
		topLevelGauss->toplevelgauss(tmp_xors, &out_changed_occur);
		//these may have changed, recalculating occur
		// for(Lit lit: out_changed_occur) {
		//     n_occurs[lit.toInt()] = calc_occ_data(lit);
		//     n_occurs[(~lit).toInt()] = calc_occ_data(~lit);
		// }
	    }
	}
	#endif
}
