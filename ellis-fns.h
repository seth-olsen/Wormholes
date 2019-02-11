#ifndef EKG_FNS_H_INCLUDED
#define EKG_FNS_H_INCLUDED

#include "fda-fns.h"
#include <map>
#include <vector> // for everything
#include <cmath> // for ICs

inline dbl r2(PAR *p, int k)
{
  return (sq(p->r[k]) + (p->lsq));
}
////////////////////////////////////////////////////////////////////////////////////////////////
// *************************************************************************
// **********************   FDA IMPLEMENTATIONS   **************************
// *************************************************************************
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_xi(FLDS *f, PAR *p, int k)
{
  return f->oldXi[k] + (p->lam2val) *
    ((f->cnAl[k+1])*(f->cnPi[k+1])/sq(f->cnPs[k+1]) + (f->cnBe[k+1])*(f->cnXi[k+1])
     - (f->cnAl[k-1])*(f->cnPi[k-1])/sq(f->cnPs[k-1]) - (f->cnBe[k-1])*(f->cnXi[k-1]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
dbl fda_pi(FLDS *f, PAR *p, int k)
{
  dbl lam6part = 1 + (p->lam6val)*(d_c(f->cnBe,k) + (f->cnBe[k])*(6*dln_c(f->cnPs,k) + 4*(p->dr)*(p->r[-k])));
  return ((f->oldPi[k])*(2 - lam6part) + ((p->lam2val) / (r2(p,k)*pw4(f->cnPs[k]))) *
	  ((r2(p,k+1)*pw4(f->cnPs[k+1]))*((f->cnAl[k+1])*(f->cnXi[k+1])/sq(f->cnPs[k+1]) + (f->cnBe[k+1])*(f->cnPi[k+1])) -
	   (r2(p,k-1)*pw4(f->cnPs[k-1]))*((f->cnAl[k-1])*(f->cnXi[k-1])/sq(f->cnPs[k-1]) + (f->cnBe[k-1])*(f->cnPi[k-1]))))
    / lam6part;
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resXi(FLDS *f, PAR *p, int k)
{
  return (p->indt)*(f->Xi[k] - f->oldXi[k]) - (p->in2dr) *
    ((f->cnAl[k+1])*(f->cnPi[k+1])/sq(f->cnPs[k+1]) + (f->cnBe[k+1])*(f->cnXi[k+1])
     - (f->cnAl[k-1])*(f->cnPi[k-1])/sq(f->cnPs[k-1]) - (f->cnBe[k-1])*(f->cnXi[k-1]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resPi(FLDS *f, PAR *p, int k)
{
  return (p->indt)*(f->Pi[k] - f->oldPi[k]) - ((p->in2dr) / (r2(p,k)*pw4(f->cnPs[k]))) *
    ((r2(p,k+1)*pw4(f->cnPs[k+1]))*((f->cnAl[k+1])*(f->cnXi[k+1])/sq(f->cnPs[k+1]) + (f->cnBe[k+1])*(f->cnPi[k+1])) -
     (r2(p,k-1)*pw4(f->cnPs[k-1]))*((f->cnAl[k-1])*(f->cnXi[k-1])/sq(f->cnPs[k-1]) + (f->cnBe[k-1])*(f->cnPi[k-1])))
    + (f->cnPi[k])*(p->in3dr)*(d_c(f->cnBe,k) + (f->cnBe[k])*(6*dln_c(f->cnPs,k) + 4*(p->dr)*(p->r[-k])));
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_xi2(FLDS *f, PAR *p, int k)
{
  return f->oldXi2[k] + (p->lam2val) *
    ((f->cnAl[k+1])*(f->cnPi2[k+1])/sq(f->cnPs[k+1]) + (f->cnBe[k+1])*(f->cnXi2[k+1])
     - (f->cnAl[k-1])*(f->cnPi2[k-1])/sq(f->cnPs[k-1]) - (f->cnBe[k-1])*(f->cnXi2[k-1]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
dbl fda_pi2(FLDS *f, PAR *p, int k)
{
  dbl lam6part = 1 + (p->lam6val)*(d_c(f->cnBe,k) + (f->cnBe[k])*(6*dln_c(f->cnPs,k) + 4*(p->dr)*(p->r[-k])));
  return ((f->oldPi2[k])*(2 - lam6part) + ((p->lam2val) / (r2(p,k)*pw4(f->cnPs[k]))) *
	  ((r2(p,k+1)*pw4(f->cnPs[k+1]))*((f->cnAl[k+1])*(f->cnXi2[k+1])/sq(f->cnPs[k+1]) + (f->cnBe[k+1])*(f->cnPi2[k+1])) -
	   (r2(p,k-1)*pw4(f->cnPs[k-1]))*((f->cnAl[k-1])*(f->cnXi2[k-1])/sq(f->cnPs[k-1]) + (f->cnBe[k-1])*(f->cnPi2[k-1]))))
    / lam6part;
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resXi2(FLDS *f, PAR *p, int k)
{
  return (p->indt)*(f->Xi2[k] - f->oldXi2[k]) - (p->in2dr) *
    ((f->cnAl[k+1])*(f->cnPi2[k+1])/sq(f->cnPs[k+1]) + (f->cnBe[k+1])*(f->cnXi2[k+1])
     - (f->cnAl[k-1])*(f->cnPi2[k-1])/sq(f->cnPs[k-1]) - (f->cnBe[k-1])*(f->cnXi2[k-1]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resPi2(FLDS *f, PAR *p, int k)
{
  return (p->indt)*(f->Pi2[k] - f->oldPi2[k]) - ((p->in2dr) / (r2(p,k)*pw4(f->cnPs[k]))) *
    ((r2(p,k+1)*pw4(f->cnPs[k+1]))*((f->cnAl[k+1])*(f->cnXi2[k+1])/sq(f->cnPs[k+1]) + (f->cnBe[k+1])*(f->cnPi2[k+1])) -
     (r2(p,k-1)*pw4(f->cnPs[k-1]))*((f->cnAl[k-1])*(f->cnXi2[k-1])/sq(f->cnPs[k-1]) + (f->cnBe[k-1])*(f->cnPi2[k-1])))
    + (f->cnPi2[k])*(p->in3dr)*(d_c(f->cnBe,k) + (f->cnBe[k])*(6*dln_c(f->cnPs,k) + 4*(p->dr)*(p->r[-k])));
}
////////////////////////////////////////////////////////////////////////////////////////////////
dbl fda_hyp_ps(const VD& old_ps, const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, int k)
{
  dbl dt12_part = (p->lam6val) * (0.25*d_c(cn_be,k) + cn_be[k]*(p->dr)*(p->r[-k]));
  return ( old_ps[k]*(1 + dt12_part) + (p->lam2val)*cn_be[k]*d_c(cn_ps,k) ) / (1 - dt12_part);
}
inline dbl fdaR_hyp_ps(const VD& f_ps, PAR *p)
{
  return (p->cpsi_rhs)*((p->inrmax) - (p->jacRRm1)*f_ps[(p->lastpt)-1] - (p->jacRRm2)*f_ps[(p->lastpt)-2]);
}
inline dbl fda0_hyp_ps(const VD& f_ps, PAR *p)
{
  return (p->cpsi_rhs)*((p->inrmax) - (p->jacRRm1)*f_ps[1] - (p->jacRRm2)*f_ps[2]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_hyp_resPs(const VD& old_ps, const VD& f_ps, const VD& cn_al,
			 const VD& cn_be, const VD& cn_ps, PAR *p, int k)
{
  return (p->indt)*(f_ps[k] - old_ps[k]) - cn_be[k]*ddr_c(cn_ps,p,k)
    - cn_ps[k]*(p->one_third)*(0.5*ddr_c(cn_be,p,k) + cn_be[k]*(p->r[-k]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resPs(const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		     const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k)
{
  return ddr2_c(f_ps,p,k) + M_PI*(sq(f_xi2[k]) + sq(f_pi2[k]) - sq(f_xi[k]) - sq(f_pi[k]))
    + f_ps[k]*(0.25*(p->lsq)/sq((p->lsq) + sq(p->r[k])) + 2*(p->r[-k])*ddr_c(f_ps,p,k))
    + (p->twelfth)*pw5(f_ps[k])*sq(ddr_c(f_be,p,k) - (p->r[-k])*f_be[k]) / sq(f_al[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resBe(const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		     const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k)
{
  return ddr2_c(f_be,p,k) - f_be[k]*(p->lsq)/sq((p->lsq) + sq(p->r[k]))
    + (p->twelve_pi)*f_al[k]*(f_xi2[k]*f_pi2[k] - f_xi[k]*f_pi[k]) / sq(f_ps[k])
    + (ddr_c(f_be,p,k) - (p->r[-k])*f_be[k])*(2*(p->r[-k]) + 6*ddrln_c(f_ps,p,k) - ddrln_c(f_al,p,k));
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resAl(const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		     const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k)
{
  return ddr2_c(f_al,p,k) + (p->eight_pi)*f_al[k]*(sq(f_pi[k]) - sq(f_pi2[k]))
    + 2*ddr_c(f_al,p,k)*((p->r[-k]) + ddrln_c(f_ps,p,k))
    - (p->two_thirds)*pw4(f_ps[k])*sq(ddr_c(f_be,p,k) - (p->r[-k])*f_be[k]) / f_al[k];
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fdaR_resPs(const VD& f_ps, PAR *p)
{ return ddr_b(f_ps,p,(p->lastpt)) + (p->inrmax)*(f_ps[(p->lastpt)] - 1); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fdaR_resBe(const VD& f_be, PAR *p)
{ return ddr_b(f_be,p,(p->lastpt)) + (p->inrmax)*f_be[(p->lastpt)]; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fdaR_resAl(const VD& f_al, PAR *p)
{ return ddr_b(f_al,p,(p->lastpt)) + (p->inrmax)*(f_al[(p->lastpt)] - 1); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda0_resPs(const VD& f_ps, PAR *p)
{ return ddr_f(f_ps,p,0) - (p->inrmax)*(f_ps[0] - 1); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda0_resBe(const VD& f_be, PAR *p)
{ return ddr_f(f_be,p,0) - (p->inrmax)*f_be[0]; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda0_resAl(const VD& f_al, PAR *p)
{ return ddr_f(f_al,p,0) - (p->inrmax)*(f_al[0] - 1); }
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
