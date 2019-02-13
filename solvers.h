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
//#include "ellis-clean.h"

using namespace std;

typedef int (*SOLVER)(VD& , VD& , VD& , VD& , VD& ,
		      VD& , VD& , VD& , VD& , VD& ,
		      VD& , VD& , VD& , VD& , VD& ,
		      VD& , VD& , const VD& , MAPID& ,
		      int , int , int ,
		      int , int , int , int , int , vector<int>& , int);

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
  if (p->write_ires_abp) { f->olderPs = f->oldPs; }
  f->oldAl = f->Al;  f->cnAl = f->Al;  
  f->oldBe = f->Be;  f->cnBe = f->Be;
  f->oldPs = f->Ps;  f->cnPs = f->Ps;
  f->oldXi = f->Xi;  f->cnXi = f->Xi;
  f->oldPi = f->Pi;  f->cnPi = f->Pi;


  int itn = (*(p->solver))(f, p);
  // ****************** ITERATIVE SOLUTION COMPLETE ******************

  // ************************ APPARENT HORIZON SEARCH *************************
  if (itn < 0) {
    cout << endl << i << "\n\nSTEP EXIT CODE " << itn << endl;
    cout << "\n\nSTEP EXIT ITN " << p->exit_itn << endl;
    cout << "\nt = " << p->t << endl;
    int horizon_code = search_for_horizon(f->Al, f->Be, f->Ps, p);
    if (horizon_code) {
      record_horizon(p, f->Ps, horizon_code, p->exit_itn, i);
      return horizon_code;
    }
    else { return itn; }
  }
  p->t += p->dt;
  return 0;
}

int solve_static(FLDS *f, PAR *p)
{
  update_xp(f, p);
  int itn = 1;
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
  //dbl ell_tol = 20 * (p->ell_tol);
  int ell_maxit = (p->maxit) * 10;

  int ell_itn = 0;
  dbl res = get_res_abp(res_0, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
  while (res > (p->ell_tol)) {
    VD jac(ldab*N_0, 0);
    set_jacCMabpslow(jac, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N_0, kl, ku, nrhs,
			 &jac[0], ldab, &ipiv_0[0], &res_0[0], ldb_0);
    if (info != 0) { cout << ell_itn << "\nERROR: cannot solve initial elliptics\ninfo = " << info << endl; }
    
    apply_up_abp(res_0, f->Al, f->Be, f->Ps, p->npts, p->ell_up_weight);
    res = get_res_abp(res_0, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
    
    if (++ell_itn > ell_maxit) {
      cout << "\nt = 0\nitn = " << ell_itn << "\nres = " << res << endl;
      int horizon_found = search_for_horizon(f->Al, f->Be, f->Ps, p);
      if (horizon_found) { return -horizon_found; }
      else if (res < p->tol) { return 0; }
      else { return 2*(p->maxit); }
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
        return -2;
      }
    }    
    get_res_abp(f->res_ell, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCMabpslow(jac, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_abp(f->res_ell, f->Al, f->Be, f->Ps, p->npts, p->ell_up_weight);
      res = get_res_abp(f->res_ell, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
      }
    }
    set_abp_cn(f->oldAl, f->oldBe, f->oldPs, f->Al, f->Be, f->Ps,
	       f->cnAl, f->cnBe, f->cnPs, p->npts);
    res = get_res_xp(f, p);
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      if (res < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
	else { res = 0; }
      }
      else { return -5; }
    }   
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);
  
  return itn;
}

int solve2_dynamic(FLDS *f, PAR *p)
{
  VD jac;
  lapack_int info = 0;
  int itn = 0, hyp_itn = 0, ell_itn = 0;
  dbl res = (p->tol) + 1;

  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->tol)) {
      update_xp2(f, p);
      res = get_res_xp2(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        return -2;
      }
    }    
    get_res_abp(f->res_ell, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCMabpslow(jac, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_abp(f->res_ell, f->Al, f->Be, f->Ps, p->npts, p->ell_up_weight);
      res = get_res_abp(f->res_ell, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
      }
    }
    set_abp_cn(f->oldAl, f->oldBe, f->oldPs, f->Al, f->Be, f->Ps,
	       f->cnAl, f->cnBe, f->cnPs, p->npts);
    res = get_res_xp2(f, p);
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      if (res < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
	else { res = 0; }
      }
      else { return -5; }
    }   
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi2, f->oldPi2, f->Xi2, f->Pi2, (p->lastpt)-1, p->dspn);
  
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
        return -2;
      }
    }
    hyp_itn = 0;
    while (res > (p->tol)) {
      update_xp2(f, p);
      res = get_res_xp2(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "NORMAL HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        return -2;
      }
    }
    get_res_abp(f->res_ell, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCMabpslow(jac, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_abp(f->res_ell, f->Al, f->Be, f->Ps, p->npts, p->ell_up_weight);
      res = get_res_abp(f->res_ell, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
      }
    }
    set_abp_cn(f->oldAl, f->oldBe, f->oldPs, f->Al, f->Be, f->Ps,
	       f->cnAl, f->cnBe, f->cnPs, p->npts);
    res = get_res_xp(f, p);
    if (res < p->tol) { res = get_res_xp2(f, p); }
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      if (res < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
	else { res = 0; }
      }
      else { return -5; }
    }   
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);
  dissipationNB2_xp(f->oldXi2, f->oldPi2, f->Xi2, f->Pi2, (p->lastpt)-1, p->dspn);
  
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
      update_psi_clean(f, p);      
      res = get_res_psi_clean(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "PSI hyperbolic solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -1; }
	  else { res = 0; }
	}
	else { return -1; }
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
        return -2;
      }
    }
    res = get_res_ab_clean(f, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCMabslow(jac, f->Pi, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_ab(f->res_ell, f->Al, f->Be, p->npts, p->ell_up_weight);
      res = get_res_ab_clean(f, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
      }
    }
    set_ab_cn(f->oldAl, f->oldBe, f->Al, f->Be, f->cnAl, f->cnBe, p->npts);
    res = get_res_psi_clean(f, p);
    if (res < (p->ell_tol)) {
      res = get_res_xp(f, p);
    }
    else { res = (p->tol) + 1; }
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      if (get_res_psi_clean(f, p) < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
	else { res = get_res_xp(f, p); }
      }
      else { return -5; }
    }
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);
  
  return itn;
}

int solve2_dynamic_psi_hyp(FLDS *f, PAR *p)
{
  VD jac;
  lapack_int info = 0;
  int itn = 0, hyp_itn = 0, ell_itn = 0;
  dbl res = (p->tol) + 1;

  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->ell_tol)) {
      update_psi_clean(f, p);      
      res = get_res_psi_clean(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "PSI hyperbolic solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -1; }
	  else { res = 0; }
	}
	else { return -1; }
      }
    }
    res = (p->tol) + 1;
    hyp_itn = 0;
    while (res > (p->tol)) {
      update_xp2(f, p);
      res = get_res_xp2(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        return -2;
      }
    }
    res = get_res_ab_clean(f, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCMabslow(jac, f->Pi, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_ab(f->res_ell, f->Al, f->Be, p->npts, p->ell_up_weight);
      res = get_res_ab_clean(f, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
      }
    }
    set_ab_cn(f->oldAl, f->oldBe, f->Al, f->Be, f->cnAl, f->cnBe, p->npts);
    res = get_res_psi_clean(f, p);
    if (res < (p->ell_tol)) {
      res = get_res_xp2(f, p);
    }
    else { res = (p->tol) + 1; }
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      if (get_res_psi_clean(f, p) < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
	else { res = get_res_xp2(f, p); }
      }
      else { return -5; }
    }
  }
  // *********************** kreiss-oliger DISSIPATION ************************
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
      update_psi_clean(f, p);      
      res = get_res_psi_clean(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "PSI hyperbolic solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -1; }
	  else { res = 0; }
	}
	else { return -1; }
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
        return -2;
      }
    }
    hyp_itn = 0;
    while (res > (p->tol)) {
      update_xp2(f, p);
      res = get_res_xp2(f, p);
      if (++hyp_itn > p->maxit) {
	cout << endl << "NORMAL HYPERBOLIC solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        return -2;
      }
    }
    res = get_res_ab_clean(f, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCMabslow(jac, f->Pi, f->Pi2, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_ab(f->res_ell, f->Al, f->Be, p->npts, p->ell_up_weight);
      res = get_res_ab_clean(f, p);
      if (++ell_itn > p->maxit) {
	cout << endl << "ELLIPTIC solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
        if (res < (p->tol)) {
	  if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -4; }
	  else { res = 0; }
	}
        else { return -4; }
      }
    }
    set_ab_cn(f->oldAl, f->oldBe, f->Al, f->Be, f->cnAl, f->cnBe, p->npts);
    res = get_res_psi_clean(f, p);
    if (res < (p->ell_tol)) {
      res = get_res_xp(f, p);
      if (res < (p->tol)) { res = get_res_xp2(f, p); }
    }
    else { res = (p->tol) + 1; }
    
    if (++itn > p->maxit) {
      cout << endl << "FULL solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      if (get_res_psi_clean(f, p) < (p->tol)) {
	if (search_for_horizon(f->Al, f->Be, f->Ps, p)) { return -5; }
        else {
	  res = get_res_xp(f, p);
	  if (res < (p->tol)) { res = get_res_xp2(f, p); }
	}
      }
      else { return -5; }
    }
  }
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);
  dissipationNB2_xp(f->oldXi2, f->oldPi2, f->Xi2, f->Pi2, (p->lastpt)-1, p->dspn);
  
  return itn;
}

#endif
