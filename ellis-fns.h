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
    + f_ps[k]*(0.25*(p->lsq)/sq(r2(p,k)) + 2*(p->r[-k])*ddr_c(f_ps,p,k))
    + (p->twelfth)*pw5(f_ps[k])*sq(ddr_c(f_be,p,k) - (p->r[-k])*f_be[k]) / sq(f_al[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl fda_resBe(const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		     const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k)
{
  return ddr2_c(f_be,p,k) - f_be[k]*(p->lsq)/sq(r2(p,k))
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

inline dbl iresxi_c(FLDS *f, PAR *p, int k)
{ 
  return (p->indt)*(1.5*(f->Xi[k]) - 2*(f->oldXi[k]) + 0.5*(f->olderXi[k])) - (p->in2dr) *
    ((f->Be[k])*d_c(f->Xi,k) + (f->Xi[k])*d_c(f->Be,k) + (f->Al[k])*d_c(f->Pi,k)/sq(f->Ps[k])
     + ((f->Pi[k])/pw3(f->Ps[k]))*((f->Ps[k])*d_c(f->Al,k) - 2*(f->Al[k])*d_c(f->Ps,k)));
}

inline dbl irespi_c(FLDS *f, PAR *p, int k)
{
  return (p->indt)*(1.5*(f->Pi[k]) - 2*(f->oldPi[k]) + 0.5*(f->olderPi[k]) +
		    4*(f->Pi[k])*(1.5*(f->Ps[k]) - 2*(f->oldPs[k]) + 0.5*(f->olderPs[k]))/(f->Ps[k]))
    - (p->in2dr)*((f->Be[k])*d_c(f->Pi,k) + (f->Pi[k])*d_c(f->Be,k) + (f->Al[k])*d_c(f->Xi,k)/sq(f->Ps[k])
		  + ((f->Xi[k])/pw3(f->Ps[k]))*((f->Ps[k])*d_c(f->Al,k) - 2*(f->Al[k])*d_c(f->Ps,k)))
    - ((f->Al[k])*(f->Xi[k])/sq(f->Ps[k]) + (f->Be[k])*(f->Pi[k]))*(4*ddrln_c(f->Ps,p,k) + 2*(p->r[-k]));
}

void get_ires_xp(WRS *wr, FLDS *f, PAR *p)
{
  (wr->p_iresXi).wr_field[0] = sommerfeldres_f(f->oldXi, f->Xi, p, 0);
  (wr->p_iresPi).wr_field[0] = sommerfeldres_f(f->oldPi, f->Pi, p, 0);
  for (int k = 1; k < p->lastwr; ++k) {
    (wr->p_iresXi).wr_field[k] = iresxi_c(f, p, (p->inds[k]).second);
    (wr->p_iresPi).wr_field[k] = irespi_c(f, p, (p->inds[k]).second);
  }
  (wr->p_iresXi).wr_field[p->lastwr] = sommerfeldres(f->oldXi, f->Xi, p, p->lastpt);
  (wr->p_iresPi).wr_field[p->lastwr] = sommerfeldres(f->oldPi, f->Pi, p, p->lastpt);
  return;
}

inline dbl iresxi2_c(FLDS *f, PAR *p, int k)
{ 
  return (p->indt)*(1.5*(f->Xi2[k]) - 2*(f->oldXi2[k]) + 0.5*(f->olderXi2[k])) - (p->in2dr) *
    ((f->Be[k])*d_c(f->Xi2,k) + (f->Xi2[k])*d_c(f->Be,k) + (f->Al[k])*d_c(f->Pi2,k)/sq(f->Ps[k])
     + ((f->Pi2[k])/pw3(f->Ps[k]))*((f->Ps[k])*d_c(f->Al,k) - 2*(f->Al[k])*d_c(f->Ps,k)));
}

inline dbl irespi2_c(FLDS *f, PAR *p, int k)
{
  return (p->indt)*(1.5*(f->Pi2[k]) - 2*(f->oldPi2[k]) + 0.5*(f->olderPi2[k]) +
		    4*(f->Pi2[k])*(1.5*(f->Ps[k]) - 2*(f->oldPs[k]) + 0.5*(f->olderPs[k]))/(f->Ps[k]))
    - (p->in2dr)*((f->Be[k])*d_c(f->Pi2,k) + (f->Pi2[k])*d_c(f->Be,k) + (f->Al[k])*d_c(f->Xi2,k)/sq(f->Ps[k])
		  + ((f->Xi2[k])/pw3(f->Ps[k]))*((f->Ps[k])*d_c(f->Al,k) - 2*(f->Al[k])*d_c(f->Ps,k)))
    - ((f->Al[k])*(f->Xi2[k])/sq(f->Ps[k]) + (f->Be[k])*(f->Pi2[k]))*(4*ddrln_c(f->Ps,p,k) + 2*(p->r[-k]));
}

void get_ires_xp2(WRS *wr, FLDS *f, PAR *p)
{
  (wr->p_iresXi2).wr_field[0] = sommerfeldres_f(f->oldXi2, f->Xi2, p, 0);
  (wr->p_iresPi2).wr_field[0] = sommerfeldres_f(f->oldPi2, f->Pi2, p, 0);
  for (int k = 1; k < p->lastwr; ++k) {
    (wr->p_iresXi2).wr_field[k] = iresxi2_c(f, p, (p->inds[k]).second);
    (wr->p_iresPi2).wr_field[k] = irespi2_c(f, p, (p->inds[k]).second);
  }
  (wr->p_iresXi2).wr_field[p->lastwr] = sommerfeldres(f->oldXi2, f->Xi2, p, p->lastpt);
  (wr->p_iresPi2).wr_field[p->lastwr] = sommerfeldres(f->oldPi2, f->Pi2, p, p->lastpt);
  return;
}

inline dbl irespsihyp_c(FLDS *f, PAR *p, int k)
{
  return fda_resPs(f->Xi, f->Pi, f->Xi2, f->Pi2, f->Al, f->Be, f->Ps, p, k);
}

inline dbl irespsi_c(FLDS *f, PAR *p, int k)
{
  return (p->indt)*(1.5*(f->Ps[k]) - 2*(f->oldPs[k]) + 0.5*(f->olderPs[k])) - f->Be[k]*ddr_c(f->Ps,p,k)
    - f->Ps[k]*(p->one_third)*(0.5*ddr_c(f->Be,p,k) + f->Be[k]*(p->r[-k]));
}

inline dbl iresbeta_c(FLDS *f, PAR *p, int k)
{
  return  ddr2_c(f->Be,p,k) - f->Be[k]*(p->lsq)/sq(r2(p,k))
    + (p->twelve_pi)*(f->Al[k])*(f->Xi2[k]*f->Pi2[k] - f->Xi[k]*f->Pi[k]) / sq(f->Ps[k])
    + sqrt(r2(p,k))*((f->Be[k+1])/sqrt(r2(p,k+1)) - (f->Be[k-1])/sqrt(r2(p,k+1)))*
    (2*(p->r[-k]) + 6*ddrlog_c(f->Ps,p,k) - ddrlog_c(f->Al,p,k));
}

inline dbl iresalpha_c(FLDS *f, PAR *p, int k)
{
  return ddr2_c(f->Al,p,k) + (p->eight_pi)*(f->Al[k])*(sq(f->Pi[k]) - sq(f->Pi2[k]))
    + (p->indr)*d_c(f->Al,k)*((p->r[-k]) + ddrlog_c(f->Ps,p,k))
    - ((p->in2dr)*(p->in3dr)*r2(p,k)*pw4(f->Ps[k])/(f->Al[k])) *
    sq((f->Be[k+1])/sqrt(r2(p,k+1)) - (f->Be[k-1])/sqrt(r2(p,k+1)));
}


inline dbl irespsi_f(FLDS *f, PAR *p, int k)
{
  return (-3*f->Ps[k] + 4*f->Ps[k+1] - f->Ps[k+2]) * (p->in2dr);
}

inline dbl iresbeta_f(FLDS *f, PAR *p, int k)
{
  return (-3*f->Be[k] + 4*f->Be[k+1] - f->Be[k+2]) * (p->in2dr);
}

inline dbl iresalpha_f(FLDS *f, PAR *p, int k)
{
  return (-3*f->Al[k] + 4*f->Al[k+1] - f->Al[k+2]) * (p->in2dr);
}

void get_ires_abp(WRS *wr, FLDS *f, PAR *p)
{
  (wr->p_iresAl).wr_field[0] = fda0_resAl(f->Al, p);
  (wr->p_iresBe).wr_field[0] = fda0_resBe(f->Be, p);
  (wr->p_iresPs).wr_field[0] = fda0_resPs(f->Ps, p);
  if (p->psi_hyp) {
    for (int k = 1; k < p->lastwr; ++k) {
      (wr->p_iresAl).wr_field[k] = iresalpha_c(f, p, (p->inds[k]).second);
      (wr->p_iresBe).wr_field[k] = iresbeta_c(f, p, (p->inds[k]).second);
      (wr->p_iresPs).wr_field[k] = irespsihyp_c(f, p, (p->inds[k]).second);
    }
  }
  else {
    for (int k = 1; k < p->lastwr; ++k) {
      (wr->p_iresAl).wr_field[k] = iresalpha_c(f, p, (p->inds[k]).second);
      (wr->p_iresBe).wr_field[k] = iresbeta_c(f, p, (p->inds[k]).second);
      (wr->p_iresPs).wr_field[k] = irespsi_c(f, p, (p->inds[k]).second);
    }
  }
  (wr->p_iresAl).wr_field[p->lastwr] = fdaR_resAl(f->Al, p);
  (wr->p_iresBe).wr_field[p->lastwr] = fdaR_resBe(f->Be, p);
  (wr->p_iresPs).wr_field[p->lastwr] = fdaR_resPs(f->Ps, p);
  return;
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
