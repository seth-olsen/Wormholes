#ifndef EKG_PROC_H_INCLUDED
#define EKG_PROC_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <vector> // for everything
#include <cmath> // for ICs
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "fda-io.h"
#include "fda-fns.h"
#include "sim-header.h"
#include "sim-structs.h"
#include "ellis-fns.h"
#include "jacobian.h"
//#include "ellis-clean.h"
#include "lapacke.h"

using namespace std;

// update functions
void update_xp(FLDS *f, PAR *p);
dbl get_res_xp(FLDS *f, PAR *p);
void update_psi_clean(FLDS *f, PAR *p);
dbl get_res_psi_clean(FLDS *f, PAR *p);
void get_res_abp_clean(FLDS *f, PAR *p);
dbl get_res_abp(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p);
void apply_up_abp(const VD& res_ell, VD& f_al, VD& f_be, VD& f_ps, int npts, dbl eup_weight);
void set_abp_cn(const VD& old_al, const VD& old_be, const VD& old_ps,
		const VD& f_al, const VD& f_be, const VD& f_ps,
		VD& cn_al, VD& cn_be, VD& cn_ps, int npts);
void apply_up_ab(const VD& res_ell, VD& f_al, VD& f_be, VD& f_ps, int npts, dbl eup_weight);
void set_ab_cn(const VD& old_al, const VD& old_be, const VD& f_al, const VD& f_be,
	       VD& cn_al, VD& cn_be, int npts);
// dissipation functions
void dissipationNB_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		      int one_past_last, dbl dspn);
void dissipationB_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		     int one_past_last, dbl dspn);
void dissipationNB2_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		       int one_past_last, dbl dspn);
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

inline void apply_bcR_xp(FLDS *f, PAR *p)
{
  sommerfeld(f->oldXi, f->Xi, p, p->lastpt);
  sommerfeld(f->oldPi, f->Pi, p, p->lastpt);
}

inline void apply_bc0_xp(FLDS *f, PAR *p)
{
  sommerfeld_f(f->oldXi, f->Xi, p, 0);
  sommerfeld_f(f->oldPi, f->Pi, p, 0);
}
inline void apply_bcR_xp2(FLDS *f, PAR *p)
{
  sommerfeld(f->oldXi2, f->Xi2, p, p->lastpt);
  sommerfeld(f->oldPi2, f->Pi2, p, p->lastpt);
}

inline void apply_bc0_xp2(FLDS *f, PAR *p)
{
  sommerfeld_f(f->oldXi2, f->Xi2, p, 0);
  sommerfeld_f(f->oldPi2, f->Pi2, p, 0);
}

void update_xp(FLDS *f, PAR *p)
{
  apply_bc0_xp(f, p);
  f->cnXi[0] = 0.5 * (f->oldXi[0] + f->Xi[0]);
  f->cnPi[0] = 0.5 * (f->oldPi[0] + f->Pi[0]);
  for (int k = 1; k < (p->lastpt); ++k) {
    // update f_xi & cn_xi
    f->Xi[k] = fda_xi(f, p, k);
    f->cnXi[k] = 0.5 * (f->oldXi[k] + f->Xi[k]);
    // update f_pi & cn_pi
    f->Pi[k] = fda_pi(f, p, k);
    f->cnPi[k] = 0.5 * (f->oldPi[k] + f->Pi[k]);
  }
  // r = R BOUNDARY
  apply_bcR_xp(f, p);
  f->cnXi[p->lastpt] = 0.5 * (f->oldXi[p->lastpt] + f->Xi[p->lastpt]);
  f->cnPi[p->lastpt] = 0.5 * (f->oldPi[p->lastpt] + f->Pi[p->lastpt]);
  return;
}

dbl get_res_xp(FLDS *f, PAR *p)
{
  f->resXi[0] = sommerfeldres_f(f->oldXi, f->Xi, p, 0);
  f->resPi[0] = sommerfeldres_f(f->oldPi, f->Pi, p, 0);
  for (int k = 1; k < (p->lastpt); ++k) {
    f->resXi[k] = fda_resXi(f, p, k);
    f->resPi[k] = fda_resPi(f, p, k);
  }
  f->resXi[p->lastpt] = sommerfeldres(f->oldXi, f->Xi, p, p->lastpt);
  f->resPi[p->lastpt] = sommerfeldres(f->oldPi, f->Pi, p, p->lastpt);
  return max(norm_inf(f->resXi), norm_inf(f->resPi));
}
void update_xp2(FLDS *f, PAR *p)
{
  apply_bc0_xp2(f, p);
  f->cnXi2[0] = 0.5 * (f->oldXi2[0] + f->Xi2[0]);
  f->cnPi2[0] = 0.5 * (f->oldPi2[0] + f->Pi2[0]);
  for (int k = 1; k < (p->lastpt); ++k) {
    // update f_xi & cn_xi
    f->Xi2[k] = fda_xi2(f, p, k);
    f->cnXi2[k] = 0.5 * (f->oldXi2[k] + f->Xi2[k]);
    // update f_pi & cn_pi
    f->Pi2[k] = fda_pi2(f, p, k);
    f->cnPi2[k] = 0.5 * (f->oldPi2[k] + f->Pi2[k]);
  }
  // r = R BOUNDARY
  apply_bcR_xp2(f, p);
  f->cnXi2[p->lastpt] = 0.5 * (f->oldXi2[p->lastpt] + f->Xi2[p->lastpt]);
  f->cnPi2[p->lastpt] = 0.5 * (f->oldPi2[p->lastpt] + f->Pi2[p->lastpt]);
  return;
}

dbl get_res_xp2(FLDS *f, PAR *p)
{
  f->resXi2[0] = sommerfeldres_f(f->oldXi2, f->Xi2, p, 0);
  f->resPi2[0] = sommerfeldres_f(f->oldPi2, f->Pi2, p, 0);
  for (int k = 1; k < (p->lastpt); ++k) {
    f->resXi2[k] = fda_resXi2(f, p, k);
    f->resPi2[k] = fda_resPi2(f, p, k);
  }
  f->resXi2[p->lastpt] = sommerfeldres(f->oldXi2, f->Xi2, p, p->lastpt);
  f->resPi2[p->lastpt] = sommerfeldres(f->oldPi2, f->Pi2, p, p->lastpt);
  return max(norm_inf(f->resXi2), norm_inf(f->resPi2));
}

void update_psi_clean(FLDS *f, PAR *p)
{
  f->Ps[0] = fda0_hyp_ps(f->Ps, p);
  f->cnPs[0] = 0.5 * (f->oldPs[0] + f->Ps[0]);
  for (int k = 1; k < (p->lastpt); ++k) {
    f->Ps[k] = fda_hyp_ps(f->oldPs, f->cnAl, f->cnBe, f->cnPs, p, k);
    f->cnPs[k] = 0.5 * (f->oldPs[k] + f->Ps[k]);
  }
  f->Ps[p->lastpt] = fdaR_hyp_ps(f->Ps, p);
  f->cnPs[p->lastpt] = 0.5 * (f->oldPs[p->lastpt] + f->Ps[p->lastpt]);
  return;
}

dbl get_res_psi_clean(FLDS *f, PAR *p)
{
  f->resPs[0] = fda0_resPs(f->Ps, p);
  for (int k = 1; k < (p->lastpt); ++k) {
    f->resPs[k] = fda_hyp_resPs(f->oldPs, f->Ps, f->cnAl, f->cnBe, f->cnPs, p, k);
  }
  f->resPs[p->lastpt] = fdaR_resPs(f->Ps, p);
  return norm_inf(f->resPs);
}

void get_res_abp_clean(FLDS *f, PAR *p)
{
  int j = 0;
  int jbe = (p->npts);
  int jps = 2*jbe;
  f->res_ell[j] = fda0_resAl(f->Al, p);
  f->res_ell[jbe + j] = fda0_resBe(f->Be, p);
  f->res_ell[jps + j] = fda0_resPs(f->Ps, p);
  for (j = 1; j < (p->lastpt); ++j) {
    f->res_ell[j] = fda_resAl(f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p, j);
    f->res_ell[jbe + j] = fda_resBe(f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p, j);
    f->res_ell[jps + j] = fda_resPs(f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p, j);
  }
  j = p->lastpt;
  f->res_ell[j] = fdaR_resAl(f->Al, p);
  f->res_ell[jbe + j] = fdaR_resBe(f->Be, p);
  f->res_ell[jps + j] = fdaR_resPs(f->Ps, p);
  return;
}

dbl get_res_abp(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p)
{
  int j = 0;
  int jbe = (p->npts);
  int jps = 2*jbe;
  res_ell[j] = fda0_resAl(f_al, p);
  res_ell[jbe + j] = fda0_resBe(f_be, p);
  res_ell[jps + j] = fda0_resPs(f_ps, p);
  for (j = 1; j < (p->lastpt); ++j) {
    res_ell[j] = fda_resAl(f_xi, f_pi, f_xi2, f_pi2, f_al, f_be, f_ps, p, j);
    res_ell[jbe + j] = fda_resBe(f_xi, f_pi, f_xi2, f_pi2, f_al, f_be, f_ps, p, j);
    res_ell[jps + j] = fda_resPs(f_xi, f_pi, f_xi2, f_pi2, f_al, f_be, f_ps, p, j);
  }
  j = p->lastpt;
  res_ell[j] = fdaR_resAl(f_al, p);
  res_ell[jbe + j] = fdaR_resBe(f_be, p);
  res_ell[jps + j] = fdaR_resPs(f_ps, p);
  return norm_inf(res_ell);
}

dbl get_res_ab_clean(FLDS *f, PAR *p)
{
  int jbe = (p->npts);
  f->res_ell[0] = fda0_resAl(f->Al, p);
  f->res_ell[jbe] = fda0_resBe(f->Be, p);
  for (int j = 1; j < (p->lastpt); ++j) {
    f->res_ell[j] = fda_resAl(f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p, j);
    f->res_ell[jbe + j] = fda_resBe(f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p, j);
  }
  f->res_ell[p->lastpt] = fdaR_resAl(f->Al, p);
  f->res_ell[jbe + (p->lastpt)] = fdaR_resBe(f->Be, p);
  return norm_inf(f->res_ell);
}

void apply_up_abp(const VD& res_ell, VD& f_al, VD& f_be, VD& f_ps, int npts, dbl eup_weight)
{
  int jps = 2*npts;
  for (int j = 0; j < npts; ++j) {
    f_al[j] -= eup_weight*res_ell[j];
    f_be[j] -= eup_weight*res_ell[npts + j];
    f_ps[j] -= eup_weight*res_ell[jps + j];
  }
  return;
}
void set_abp_cn(const VD& old_al, const VD& old_be, const VD& old_ps,
		const VD& f_al, const VD& f_be, const VD& f_ps,
		VD& cn_al, VD& cn_be, VD& cn_ps, int npts)
{
  for (int k = 0; k < npts; ++k) {
    cn_al[k] = 0.5 * (old_al[k] + f_al[k]);
    cn_be[k] = 0.5 * (old_be[k] + f_be[k]);
    cn_ps[k] = 0.5 * (old_ps[k] + f_ps[k]);
  }
  return;
}

void apply_up_ab(const VD& res_ell, VD& f_al, VD& f_be, int npts, dbl eup_weight)
{
  for (int j = 0; j < npts; ++j) {
    f_al[j] -= eup_weight*res_ell[j];
    f_be[j] -= eup_weight*res_ell[npts + j];
  }
  return;
}
void set_ab_cn(const VD& old_al, const VD& old_be, const VD& f_al, const VD& f_be,
		VD& cn_al, VD& cn_be, int npts)
{
  for (int k = 0; k < npts; ++k) {
    cn_al[k] = 0.5 * (old_al[k] + f_al[k]);
    cn_be[k] = 0.5 * (old_be[k] + f_be[k]);
  }
  return;
}
  
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// at ind next to boundaries can call dissipate on ind+/-1 or ind+/-2
void dissipationNB_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		      int one_past_last, dbl dspn)
{
  f_xi[1] += antidiss1(dspn, old_xi);
  f_pi[1] += symdiss1(dspn, old_pi);
  for (int k = 2; k < one_past_last; ++k) {
    f_xi[k] += dissipate(dspn, old_xi, k);
    f_pi[k] += dissipate(dspn, old_pi, k);
  }
  return;
}
void dissipationNB2_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		       int one_past_last, dbl dspn)
{
  for (int k = 2; k < one_past_last; ++k) {
    f_xi[k] += dissipate(dspn, old_xi, k);
    f_pi[k] += dissipate(dspn, old_pi, k);
  }
  return;
}
//************ CHECKED: Mon. 10/22 ************
void dissipationB_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		     int one_past_last, dbl dspn)
{
  f_pi[0] += symdiss0(dspn, old_pi);
  f_xi[1] += antidiss1(dspn, old_xi);
  f_pi[1] += symdiss1(dspn, old_pi);
  for (int k = 2; k < one_past_last; ++k) {
    f_xi[k] += dissipate(dspn, old_xi, k);
    f_pi[k] += dissipate(dspn, old_pi, k);
  }
  return;
}
  
#endif

