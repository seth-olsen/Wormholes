#ifndef SOLVERS_H_INCLUDED
#define SOLVERS_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <vector> // for everything
#include <cmath> // for ICs
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "lapacke.h"
#include "jacobian.h"
#include "sim-header.h"
#include "sim-structs.h"
#include "fda-io.h"
#include "fda-fns.h"
#include "ellis-fns.h"
#include "ellis-proc.h"

using namespace std;

int solve_static(FLDS *f, PAR *p);
int solve2_static(FLDS *f, PAR *p);
int solveAll_static(FLDS *f, PAR *p);
int solve_t0(FLDS *f, PAR *p);
int solve_dynamic(FLDS *f, PAR *p);
int solve2_dynamic(FLDS *f, PAR *p);
int solveAll_dynamic(FLDS *f, PAR *p);
int solve_dynamic_psi_hyp(FLDS *f, PAR *p);
int solve2_dynamic_psi_hyp(FLDS *f, PAR *p);
int solveAll_dynamic_psi_hyp(FLDS *f, PAR *p);

int fields_step(FLDS *f, PAR *p, int i)
{
  if (p->write_ires_xp) {
    f->olderXi = f->oldXi;
    f->olderPi = f->oldPi;
  }
  if (p->write_ires_xp2) {
    f->olderXi2 = f->oldXi2;
    f->olderPi2 = f->oldPi2;
  }
  if (p->write_ires_abp || p->write_mtot) {
    f->oldestPs = f->olderPs;
    
    f->olderAl = f->oldAl;
    f->olderBe = f->oldBe;
    f->olderPs = f->oldPs;
  }
  f->oldAl = f->Al;  f->cnAl = f->Al;  
  f->oldBe = f->Be;  f->cnBe = f->Be;
  f->oldPs = f->Ps;  f->cnPs = f->Ps;
  f->oldXi = f->Xi;  f->cnXi = f->Xi;
  f->oldPi = f->Pi;  f->cnPi = f->Pi;
  f->oldXi2 = f->Xi2;  f->cnXi2 = f->Xi2;
  f->oldPi2 = f->Pi2;  f->cnPi2 = f->Pi2;


  int itn = (*(p->solver))(f, p);
  // ****************** ITERATIVE SOLUTION COMPLETE ******************

  // ************************ APPARENT HORIZON SEARCH *************************
  if (itn < 0) {
    cout << endl << i << "\n\nSTEP EXIT CODE " << itn << endl;
    cout << "\n\nSTEP EXIT ITN " << p->exit_itn << endl;
    cout << "\nt = " << p->t << endl;
    int horizon_code = search_for_horizon(f->Al, f->Be, f->Ps, p);
    if (horizon_code) {
      if (horizon_code == (p->npts)) { record_horizon(p, f->Ps, 0, p->exit_itn, i); }
      else { record_horizon(p, f->Ps, horizon_code, p->exit_itn, i); }
      return horizon_code;
    }
    else { horizon_code = search_for_antihorizon(f->Al, f->Be, f->Ps, p); }
    
    if (horizon_code) {
      record_antihorizon(p, f->Ps, horizon_code, p->exit_itn, i);
      return horizon_code;
    }
    else {
      record_exit(p, itn, p->exit_itn, i);
      return itn;
    }
  }
  
  p->t += p->dt;
  return 0;
}

int solve_static(FLDS *f, PAR *p)
{
  int itn = 1;
  update_xp(f, p);
  dbl res = get_res_xp(f, p);
  while (res > p->tol) {
    update_xp(f, p);
    res = get_res_xp(f, p);
    if (++itn == p->maxit) { return -1; }
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);
  
  return itn;
}
int solve2_static(FLDS *f, PAR *p)
{
  update_xp2(f, p);
  int itn = 1;
  dbl res = get_res_xp2(f, p);
  while (res > p->tol) {
    update_xp2(f, p);
    res = get_res_xp2(f, p);
    if (++itn == p->maxit) { return -1; }
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi2, f->oldPi2, f->Xi2, f->Pi2, (p->lastpt)-1, p->dspn);
  
  return itn;
}
int solveAll_static(FLDS *f, PAR *p)
{
  int itn = solve_static(f, p);
  if (itn < 0) {
    cout << "\nGHOST FIELD ERROR" << endl;
    return itn;
  }
  else {
    itn = solve2_static(f, p);
    if (itn < 0) {
      cout << "\nNORMAL FIELD ERROR" << endl;
    }
  }
  return itn;
}

int solve_t0(FLDS *f, PAR *p)
{
  lapack_int N_0 = 3*(p->npts);
  lapack_int kl = 2;
  lapack_int ku = 2;
  lapack_int nrhs = 1;
  lapack_int ldab = 2*kl + ku + 1;
  lapack_int ldb_0 = N_0;
  vector<lapack_int> ipiv_0(N_0);
  VD res_0(ldb_0, 0);
  lapack_int info = 0;
  dbl ell_tol = (p->ell_tol);
  int ell_maxit = (p->maxit) * 10;
  dbl ell_up_weight = (p->ell_up_weight);

  int ell_itn = 0;
  int ell_cycle = 0;
  dbl res = get_res_abp_t0(res_0, f, p);
  while (res > ell_tol) {
    VD jac(ldab*N_0, 0);
    set_jacCM_abp(jac, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N_0, kl, ku, nrhs,
			 &jac[0], ldab, &ipiv_0[0], &res_0[0], ldb_0);
    if (info != 0) { cout << ell_itn << "\nERROR: cannot solve t0\ninfo = " << info << endl; }
    apply_up_abp(res_0, f->Al, f->Be, f->Ps, p->npts, ell_up_weight);
    res = get_res_abp_t0(res_0, f, p);    
    if (++ell_itn > ell_maxit) {
      ++ell_cycle;
      cout << "\nt = 0\nitn = " << ell_itn << "\nres = " << res << endl;
      int horizon_code = search_for_horizon(f->Al, f->Be, f->Ps, p);
      if (horizon_code) {
	cout << "\n\n***HORIZON IN INITIAL DATA***\n" << endl;
	if (horizon_code == (p->npts)) { record_horizon(p, f->Ps, 0, 0, 0); }
	else { record_horizon(p, f->Ps, horizon_code, 0, 0); }
      }
      
      if (res < (p->tol)) { return ell_itn; }
      else if (ell_cycle > 5) { return -ell_cycle; }
      else {
	cout << "\nt=0 enering CYCLE " << ell_cycle << endl;
	ell_itn = 0;
	reset_metric(f, p);
	res = get_res_abp_t0(res_0, f, p);
	if (ell_cycle == 1) { ell_up_weight = 0.5*(p->ell_up_weight); }
	else if (ell_cycle == 2) { ell_up_weight = 2*(p->ell_up_weight); }
	else if (ell_cycle == 3) {
	  ell_up_weight = (p->ell_up_weight);
	  relax_tol(p);
	}
	else if (ell_cycle == 4) { ell_up_weight = 0.5*(p->ell_up_weight); }
	else if (ell_cycle == 5) { ell_up_weight = 2*(p->ell_up_weight); }
	else {
	  cout << "\nUNDETERMINED ERROR IN T = 0 ELLIPTIC CONSTRAINTS" << endl;
	  return -2;
	}
      }
    }
  }
  return ell_itn;
}

int solve_dynamic(FLDS *f, PAR *p)
{
  VD jac;
  lapack_int info = 0;
  int itn = 0, hyp_itn = 0, ell_itn = 0;
  dbl res = (p->tol) + 1;

  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->tol)) {
      update_xp(f, p);
      res = get_res_xp(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        //return -2;
	if (relax_tol(p) != 0) { return -2; }
	// new
      }
    }    
    res = get_res_abp(f, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCM_abp(jac, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_abp(f->res_ell, f->Al, f->Be, f->Ps, p->npts, p->ell_up_weight);
      res = get_res_abp(f, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	/*
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
	*/
	if (relax_tol(p) != 0) { return -4; }
	// new
      }
    }
    set_abp_cn(f->oldAl, f->oldBe, f->oldPs, f->Al, f->Be, f->Ps,
	       f->cnAl, f->cnBe, f->cnPs, p->npts);
    res = get_res_xp(f, p);
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      /*
      if (res < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
	else { res = 0; }
      }
      else { return -5; }
      */
      if (relax_tol(p) != 0) { return -5; }
      // new
    }   
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);
  
  return itn;
}


int solve_dynamic_psi_hyp(FLDS *f, PAR *p)
{
  VD jac;
  lapack_int info = 0;
  int itn = 0, hyp_itn = 0, ell_itn = 0;
  dbl res = (p->tol) + 1;

  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->ell_tol)) {
      update_psi(f, p);      
      res = get_res_psi(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "PSI hyperbolic solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	/*
	if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -1; }
	  else { res = 0; }
	}
	else { return -1; }
	*/
	if (relax_tol(p) != 0) { return -1; }
	// new
      }
    }
    res = (p->tol) + 1;
    hyp_itn = 0;
    while (res > (p->tol)) {
      update_xp(f, p);
      res = get_res_xp(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        // return -2;
	if (relax_tol(p) != 0) { return -2; }
	// new
      }
    }
    res = get_res_ab(f, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCM_ab(jac, f->Pi, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_ab(f->res_ell, f->Al, f->Be, p->npts, p->ell_up_weight);
      res = get_res_ab(f, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	/*
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
	*/
	if (relax_tol(p) != 0) { return -4; }
	// new
      }
    }
    set_ab_cn(f->oldAl, f->oldBe, f->Al, f->Be, f->cnAl, f->cnBe, p->npts);
    res = get_res_psi(f, p);
    if (res < (p->ell_tol)) {
      res = get_res_xp(f, p);
    }
    else { res = (p->tol) + 1; }
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      /*
      if (get_res_psi(f, p) < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
	else { res = get_res_xp(f, p); }
      }
      else { return -5; }
      */
      if (relax_tol(p) != 0) { return -5; }
      // new
    }
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);
  
  return itn;
}

int solveAll_dynamic(FLDS *f, PAR *p)
{
  VD jac;
  lapack_int info = 0;
  int itn = 0, hyp_itn = 0, ell_itn = 0;
  dbl res = (p->tol) + 1;

  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->tol)) {
      update_xp(f, p);
      res = get_res_xp(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "GHOST HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        //return -1;
	if (relax_tol(p) != 0) { return -1; }
	// new
      }
    }
    res = (p->tol) + 1;
    hyp_itn = 0;
    while (res > (p->tol)) {
      update_xp2(f, p);
      res = get_res_xp2(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "NORMAL HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        //return -2;
	if (relax_tol(p) != 0) { return -2; }
	// new
      }
    }
    res = get_res_abp(f, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCM_abp(jac, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_abp(f->res_ell, f->Al, f->Be, f->Ps, p->npts, p->ell_up_weight);
      res = get_res_abp(f, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	/*
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
	*/
	if (relax_tol(p) != 0) { return -4; }
	// new
      }
    }
    set_abp_cn(f->oldAl, f->oldBe, f->oldPs, f->Al, f->Be, f->Ps,
	       f->cnAl, f->cnBe, f->cnPs, p->npts);
    res = get_res_xp(f, p);
    if (res < p->tol) { res = get_res_xp2(f, p); }
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      /*
      if (res < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
	else { res = 0; }
      }
      else { return -5; }
      */
      if (relax_tol(p) != 0) { return -5; }
      // new
    }   
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);
  dissipationNB2_xp(f->oldXi2, f->oldPi2, f->Xi2, f->Pi2, (p->lastpt)-1, p->dspn);
  
  return itn;
}


int solveAll_dynamic_psi_hyp(FLDS *f, PAR *p)
{
  VD jac;
  lapack_int info = 0;
  int itn = 0, hyp_itn = 0, ell_itn = 0;
  dbl res = (p->tol) + 1;

  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->ell_tol)) {
      update_psi(f, p);      
      res = get_res_psi(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "PSI hyperbolic solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	/*
	if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -1; }
	  else { res = 0; }
	}
	else { return -1; }
	*/
	if (relax_tol(p) != 0) { return -1; }
	// new
      }
    }
    res = (p->tol) + 1;
    hyp_itn = 0;
    while (res > (p->tol)) {
      update_xp(f, p);
      res = get_res_xp(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "GHOST HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        //return -2;
	if (relax_tol(p) != 0) { return -2; }
	// new
      }
    }
    res = (p->tol) + 1;
    hyp_itn = 0;
    while (res > (p->tol)) {
      update_xp2(f, p);
      res = get_res_xp2(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "NORMAL HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        //return -20;
	if (relax_tol(p) != 0) { return -20; }
	// new
      }
    }
    res = get_res_ab(f, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCM_ab(jac, f->Pi, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_ab(f->res_ell, f->Al, f->Be, p->npts, p->ell_up_weight);
      res = get_res_ab(f, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	/*
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
	*/
	if (relax_tol(p) != 0) { return -4; }
	// new
      }
    }
    set_ab_cn(f->oldAl, f->oldBe, f->Al, f->Be, f->cnAl, f->cnBe, p->npts);
    res = get_res_psi(f, p);
    if (res < (p->ell_tol)) {
      res = get_res_xp(f, p);
      if (res < (p->tol)) { res = get_res_xp2(f, p); }
    }
    else { res = (p->tol) + 1; }
    
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      /*
      if (get_res_psi(f, p) < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
        else if (get_res_xp(f, p) < (p->tol)) {
	  if (get_res_xp2(f, p) < (p->tol)) { res = 0; }
	  else { return -5; }
	}
	else { return -5; }
      }
      else { return -5; }
      */
      if (relax_tol(p) != 0) { return -5; }
      // new
    }
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);
  dissipationNB2_xp(f->oldXi2, f->oldPi2, f->Xi2, f->Pi2, (p->lastpt)-1, p->dspn);
  
  return itn;
}

#endif
