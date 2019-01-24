#ifndef EKG_FNS_H_INCLUDED
#define EKG_FNS_H_INCLUDED

#include "fda-fns.h"
#include <map>
#include <vector> // for everything
#include <cmath> // for ICs


////////////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************
// **********************   FDA IMPLEMENTATIONS   **************************
// *************************************************************************
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_xi(FLDS *f, PAR *p, int k)
{
  return f->oldXi[k] + (p->lam2val)*d_c(f->cnPi, k);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_pi(FLDS *f, PAR *p, int k)
{
  return f->oldPi[k] + (p->lam2val)*d_c(f->cnXi, k) + (p->dt)*(p->r[-k])*(f->cnXi[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
// *****************************CHECK
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resXi(FLDS *f, PAR *p, int k)
{
  return (p->indt)*(f->Xi[k] - f->oldXi[k]) - ddr_c(f->cnPi,p,k);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resPi(FLDS *f, PAR *p, int k)
{
  return (p->indt)*(f->Pi[k] - f->oldPi[k]) - ddr_c(f->cnXi,p,k) - (p->r[-k])*(f->cnXi[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_hyp_ps(const VD& old_ps, const VD& cn_xi, const VD& cn_pi, const VD& cn_al,
		      const VD& cn_be, const VD& cn_ps, PAR *p, int k)
{
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// **********************************************************
// **********************************************************
//             INITIAL AND BOUNDARY CONDITIONS
// **********************************************************
// **********************************************************
////////////////////////////////////////////////////////////////////////////////////////////////
// for gaussian field or sin(coeff*r)*cos(coeff*t)/(coeff*r)
inline dbl ic_sol(dbl r, dbl amp, dbl dsq, dbl r0)
{ return amp * exp(-(r - r0)*(r - r0)/dsq); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl ic_xi(dbl r, dbl amp, dbl dsq, dbl r0)
{ return -2 * (r - r0) * amp * exp(-(r - r0)*(r - r0)/dsq) / dsq; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl ic_pi(dbl r, dbl amp, dbl dsq, dbl r0)
{ return ic_xi(r, amp, dsq, r0) + ic_sol(r, amp, dsq, r0)/r; }
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// ******** IRES FUNCTIONS *************************

// *****IMPORTANT NOTE:************
// anything with sq(psi) or /sq() may be wrong
// had to replace sqin and did so in a rush and not sure now if the replacements were right


// centered ires_f1 = f1 - ires_c(f1, f2, k, c1, d1, c2, d2)
//                   - oldf1 - ires_c(oldf1, oldf2, k, c1, d1, c2, d2)
inline dbl iresxi_c(const VD& older_xi, const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be,
		    const VD& old_ps, const VD& f_xi, int k, dbl lam, dbl al_ps2)
{ 
  return 0;
}

inline dbl irespi_c(const VD& older_pi, const VD& old_xi, const VD& old_pi, const VD& old_al, const VD& old_be,
		    const VD& old_ps, const VD& f_pi, int k, dbl lam, dbl dr, dbl r, dbl al_ps2, dbl lam6_part)
{
  return 0;
}

inline dbl irespsihyp_c(const VD& older_ps, const VD& old_ps, const VD& f_ps, int k, dbl lam6_part)
{ return 0; }

inline dbl irespsi_c(const VD& xi, const VD& pi,
		     const VD& alpha, const VD& beta,
		     const VD& psi, int k, dbl lam, dbl dr, dbl r)
{
  return (psi[k+1] - psi[k-1]) / (2*dr);
    //d2_c(psi,k) + sq(dr)*psi[k]*M_PI*(sq(xi[k]) + sq(pi[k])) + dr*d_c(psi,k)/r
    //+ pw5(psi[k])*sq(r*d_urinv_c(beta,k,dr,r)) / (48*sq(alpha[k]));
}

inline dbl iresbeta_c(const VD& xi, const VD& pi,
		      const VD& alpha, const VD& beta,
		      const VD& psi, int k, dbl lam, dbl dr, dbl r)
{
  return (beta[k+1] - beta[k-1]) / (2*dr);
}

inline dbl iresalpha_c(const VD& xi, const VD& pi,
		       const VD& alpha, const VD& beta,
		       const VD& psi, int k, dbl lam, dbl dr, dbl r)
{
  return (alpha[k+1] - alpha[k-1]) / (2*dr);
}


inline dbl irespsi_f(const VD& xi, const VD& pi,
		     const VD& alpha, const VD& beta,
		     const VD& psi, int k, dbl lam, dbl dr, dbl r)
{
  return (-3*psi[k] + 4*psi[k+1] - psi[k+2]) / (2*dr);
}

inline dbl iresbeta_f(const VD& xi, const VD& pi,
		      const VD& alpha, const VD& beta,
		      const VD& psi, int k, dbl lam, dbl dr, dbl r)
{
  return (-3*beta[k] + 4*beta[k+1] - beta[k+2]) / (2*dr);
}

inline dbl iresalpha_f(const VD& xi, const VD& pi,
		       const VD& alpha, const VD& beta,
		       const VD& psi, int k, dbl lam, dbl dr, dbl r)
{
  return (-3*alpha[k] + 4*alpha[k+1] - alpha[k+2]) / (2*dr);
}

// *********************************************
// **************  DIAGNOSTICS  ****************
// *********************************************
inline dbl mass_aspect(const VD& alpha, const VD& beta, const VD& psi, PAR *p, int k)
{
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl mass_aspectR(const VD& alpha, const VD& beta, const VD& psi, PAR *p, int k)
{
  return 0;
}
inline void get_maspect(VD& maspect, const VD& f_al, const VD& f_be, const VD& f_ps,
			PAR *p, int i_last, int i_save) {
  int s = i_save;
  for (int k = 1; k < i_last; ++k) {
    maspect[k] = mass_aspect(f_al, f_be, f_ps, p, s);
    s += i_save;
  }
  maspect[i_last] = mass_aspectR(f_al, f_be, f_ps, p, s);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_null(const VD& alpha, const VD& beta,
			 const VD& psi, PAR *p, int k)
{
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_null_f(const VD& alpha, const VD& beta,
			   const VD& psi, PAR *p, int k)
{
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_null_b(const VD& alpha, const VD& beta,
			   const VD& psi, PAR *p, int k)
{
  return 0;
}
inline void get_outnull(VD& outnull, const VD& f_al, const VD& f_be, const VD& f_ps,
			PAR *p, int i_last, int i_save) {
  outnull[0] = outgoing_null_f(f_al, f_be, f_ps, p, 0);
  int s = i_save;
  for (int k = 1; k < i_last; ++k) {
    outnull[k] = outgoing_null(f_al, f_be, f_ps, p, s);
    s += i_save;
  }
  outnull[i_last] = outgoing_null_b(f_al, f_be, f_ps, p, s);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl sRicci(const VD& f_xi, const VD& f_pi, const VD& f_ps, int k)
{
  return (sq(f_xi[k]) - sq(f_pi[k])) / pw4(f_ps[k]);
}
void get_ricci(VD& ricci, const VD& f_xi, const VD& f_pi, const VD& f_ps, const vector< pair<int,int> >& indices)
{
  for (auto k : indices) { ricci[k.first] = sRicci(f_xi, f_pi, f_ps, k.second); }
}
////////////////////////////////////////////////////////////////////////////////////
// HORIZON SEARCH
int search_for_horizon(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p)
{
  if (outgoing_null_b(f_al, f_be, f_ps, p, p->lastpt) <= 0) { return p->lastpt; }
  int k = p->lastpt;
  while (--k > 0) { if (outgoing_null(f_al, f_be, f_ps, p, k) <= 0) { return k; } }
  if (outgoing_null_f(f_al, f_be, f_ps, p, 0) <= 0) { return p->npts; }
  return 0;
}

#endif
