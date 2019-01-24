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

