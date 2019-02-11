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
int solve_dynamic(FLDS *f, PAR *p);
int solve_t0(FLDS *f, PAR *p);

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
  
  // *********************** kreiss-oliger DISSIPATION ************************
  dissipationNB2_xp(f->oldXi, f->oldPi, f->Xi, f->Pi, (p->lastpt)-1, p->dspn);

  // ************************ APPARENT HORIZON SEARCH *************************
  if (itn < 0) {
    cout << endl << i << "\n\nSTEP EXIT CODE " << itn << endl;
    cout << "\n\nSTEP EXIT ITN " << p->exit_itn << endl;
    cout << "\nt = " << p->t << endl;
    /*
    int horizon_code = search_for_horizon(f->Al, f->Be, f->Ps, p);
    if (horizon_code) {
      record_horizon(p, f->Ps, horizon_code, p->exit_itn, i);
      return horizon_code;
    }
    else { return -(p->maxit); }
    */
    return itn;
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
  dbl res = get_res_abp(res_0, f->Xi, f->Pi, f->Al, f->Be, f->Ps, p);
  while (res > (p->ell_tol)) {
    VD jac(ldab*N_0, 0);
    set_jacCMabpslow(jac, f->Xi, f->Pi, f->Al, f->Be, f->Ps, p);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N_0, kl, ku, nrhs,
			 &jac[0], ldab, &ipiv_0[0], &res_0[0], ldb_0);
    if (info != 0) { cout << ell_itn << "\nERROR: cannot solve initial elliptics\ninfo = " << info << endl; }
    
    apply_up_abp(res_0, f->Al, f->Be, f->Ps, p->npts, p->ell_up_weight);
    res = get_res_abp(res_0, f->Xi, f->Pi, f->Al, f->Be, f->Ps, p);
    
    if (++ell_itn > ell_maxit) {
      cout << "\nt = 0\nitn = " << ell_itn << "\nres = " << res << endl;
      int k = p->lastpt;
      if (outgoing_null_b(f->Al, f->Be, f->Ps, p, k) <= 0) { return -k; }
      while (--k > 0) {
	if (outgoing_null(f->Al, f->Be, f->Ps, p, k) <= 0) { return -k; }
      }
      if (outgoing_null_f(f->Al, f->Be, f->Ps, p, 0) <= 0) { return -(p->npts); }
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
    get_res_abp(f->res_ell, f->Xi, f->Pi, f->Al, f->Be, f->Ps, p);
    while (res > (p->ell_tol)) {
      jac = f->jac;
      set_jacCMabpslow(jac, f->Xi, f->Pi, f->Al, f->Be, f->Ps, p);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, p->lp_n, p->lp_kl, p->lp_ku, p->lp_nrhs,
			   &jac[0], p->lp_ldab, p->lp_ipiv, &(f->res_ell[0]), p->lp_ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "itn/ell_itn = " << itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      apply_up_abp(f->res_ell, f->Al, f->Be, f->Ps, p->npts, p->ell_up_weight);
      res = get_res_abp(f->res_ell, f->Xi, f->Pi, f->Al, f->Be, f->Ps, p);
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
  return itn;
}
  
/*

int solve_Hsearch(VD& old_xi, VD& old_pi, VD& old_al, VD& old_be, VD& old_ps,
		  VD& f_xi, VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		  VD& cn_xi, VD& cn_pi, VD& cn_al, VD& cn_be, VD& cn_ps,
		  VD& res_hyp, VD& res_ell, const VD& jac_zero, PAR *p, MAPID& r,
		  int lastpt, int maxit, int i,
		  int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
VD jac = jac_zero;
  lapack_int info = 0;
  int itn = 0, hyp_itn = 0, ell_itn = 0;
  dbl res = (p->tol) + 1;

  int kps = 2 * p->npts;
  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->ell_tol)) {
      //hyp_solve_PSI_ONLY(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);      
      res_hyp[kps] = neumann0res(f_ps, p);
      for (int k = 1; k < lastpt; ++k) {
	res_hyp[kps + k] = fda_hyp_resPs(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
      }
      res_hyp[kps + lastpt] = fdaR_hyp_resPs(f_ps, p, lastpt);
      res = max(  *max_element(res_hyp.begin() + kps, res_hyp.end()),
		  -(*min_element(res_hyp.begin() + kps, res_hyp.end()))  );
      if (++hyp_itn > maxit) {
	cout << endl << i << " PSI hyperbolic solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
	if (res < (p->tol)) {
	  //if (search_for_horizon(f_al, f_be, f_ps, p)) { return -1; }
	  //else { res = 0; }
	  ;
	}
	else { return -1; }
      }
    }
    res = (p->tol) + 1;
    hyp_itn = 0;
    while (res > (p->tol)) {
      //hyp_solve_px(old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
      //res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
      if (++hyp_itn > maxit) {
	cout << endl << i << " hyperbolic solver STUCK at t = " << (p->t) << endl;
        cout << "res = " << res << endl;
	p->exit_itn = itn;
        return -2;
      }
    }    
    //get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
    res = norm_inf(res_ell);
    while (res > (p->ell_tol)) {
      jac = jac_zero;
      set_jacCM_ab_slow(jac, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, p->npts, kl, ku, ldab);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs,
			   &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
      if (info != 0) {
	cout << (p->t) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl;
	cout << "t_itn/itn/ell_itn = " << i <<"/"<< itn <<"/" << ell_itn << "\nres = " << res << endl;
	p->exit_itn = itn;
	return -3;
      }
      //apply_up_ab(res_ell, f_al, f_be, p, p->npts);
      //get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
      res = norm_inf(res_ell);
      if (++ell_itn > maxit) {
	cout << endl << i << " elliptic solver STUCK at t = " << (p->t) << endl;
	cout << "res = " << res << endl;
	p->exit_itn = itn;
        if (res < (p->tol)) {
	  //if (search_for_horizon(f_al, f_be, f_ps, p)) { return -4; }
	  //else { res = 0; }
	  ;
	}
        else { return -4; }
      }
    }
      set2_cn(old_al, old_be, f_al, f_be, cn_al, cn_be, p->npts);
    for (int k = 1; k < lastpt; ++k) {
      res_hyp[kps + k] = fda_hyp_resPs(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    }
    res =  max(  *max_element(res_hyp.begin() + kps, res_hyp.end()),
		   -(*min_element(res_hyp.begin() + kps, res_hyp.end()))  );
    if (res < (p->ell_tol)) { 
      res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
    }
    else { res = (p->tol) + 1; }
    if (++itn > maxit) {
      res = max(norm_inf(res_ell), norm_inf(res_hyp));
      cout << endl << i << " solver STUCK at t = " << (p->t) << "\nres = " << res << endl;
      p->exit_itn = itn;
      if (res < (p->tol)) {
	//if (search_for_horizon(f_al, f_be, f_ps, p)) { return -5; }
	//else { res = 0; }
	;
      }
      else { return -5; }
    }   
  }  
  return itn;
  
}
*/

#endif
