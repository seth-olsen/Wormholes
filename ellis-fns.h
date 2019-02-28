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
  return ddr2_c(f_ps,p,k) + 2*(p->r[-k])*ddr_c(f_ps,p,k) + f_ps[k] *
    (0.25*(p->lsq)/sq(r2(p,k)) + M_PI*(sq(f_xi2[k]) + sq(f_pi2[k]) - sq(f_xi[k]) - sq(f_pi[k])))
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
  return  ddr2_c(f->Be,p,k) - (f->Be[k])*(p->lsq)/sq(r2(p,k))
    + (p->twelve_pi)*(f->Al[k])*((f->Xi2[k])*(f->Pi2[k]) - (f->Xi[k])*(f->Pi[k])) / sq(f->Ps[k])
    + sqrt(r2(p,k))*((f->Be[k+1])/sqrt(r2(p,k+1)) - (f->Be[k-1])/sqrt(r2(p,k-1)))*(p->in2dr)*
    (2*(p->r[-k]) + 6*ddrlog_c(f->Ps,p,k) - ddrlog_c(f->Al,p,k));
}

inline dbl iresalpha_c(FLDS *f, PAR *p, int k)
{
  return ddr2_c(f->Al,p,k) + (p->eight_pi)*(f->Al[k])*(sq(f->Pi[k]) - sq(f->Pi2[k]))
    + (p->indr)*d_c(f->Al,k)*((p->r[-k]) + ddrlog_c(f->Ps,p,k))
    - ((p->in2dr)*(p->in3dr)*r2(p,k)*pw4(f->Ps[k])/(f->Al[k])) *
    sq((f->Be[k+1])/sqrt(r2(p,k+1)) - (f->Be[k-1])/sqrt(r2(p,k-1)));
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
  return 0.5*sqrt(r2(p,k))*sq(psi[k]) *
    ( 1 + r2(p,k)*( sq(sq(psi[k])*(ddr_c(beta,p,k) - (p->r[-k])*beta[k])/(3*alpha[k]))
		    - sq((p->r[-k]) + (p->indr)*dln_c(psi,k)) ) );
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl mass_aspectR(const VD& alpha, const VD& beta, const VD& psi, PAR *p, int k)
{
  return 0.5*sqrt(r2(p,k))*sq(psi[k]) *
    ( 1 + r2(p,k)*( sq(sq(psi[k])*(ddr_b(beta,p,k) - (p->r[-k])*beta[k])/(3*alpha[k]))
		   - sq((p->r[-k]) + (p->indr)*dln_b(psi,k)) ) );
}
inline dbl mass_aspect0(const VD& alpha, const VD& beta, const VD& psi, PAR *p)
{
  return 0.5*sqrt(r2(p,0))*sq(psi[0]) *
    ( 1 + r2(p,0)*( sq(sq(psi[0])*(ddr_f(beta,p,0) + (p->r[-(p->lastpt)])*beta[0])/(3*alpha[0]))
		    - sq((p->indr)*dln_f(psi,0) - (p->r[-(p->lastpt)])) ) );
}
inline void get_maspect(VD& maspect, const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p) {
  maspect[0] = mass_aspect0(f_al, f_be, f_ps, p);
  for (int k = 1; k < p->lastwr; ++k) {
    maspect[k] = mass_aspect(f_al, f_be, f_ps, p, (p->inds[k]).second);
  }
  maspect[p->lastwr] = mass_aspectR(f_al, f_be, f_ps, p, p->lastpt);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl half_outgoing_null(const VD& alpha, const VD& beta,
			      const VD& psi, PAR *p, int k)
{
  return (ddr_c(beta,p,k) - (p->r[-k])*beta[k])/(3*alpha[k])
    + ((p->r[-k]) + (p->indr)*dln_c(psi,k))/sq(psi[k]); 
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl half_outgoing_null_rev(const VD& alpha, const VD& beta,
				  const VD& psi, PAR *p, int k)
{
  return (ddr_c(beta,p,k) - (p->r[-k])*beta[k])/(3*alpha[k])
    - ((p->r[-k]) + (p->indr)*dln_c(psi,k))/sq(psi[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl half_outgoing_null0(const VD& alpha, const VD& beta,
			       const VD& psi, PAR *p)
{
  return (ddr_f(beta,p,0) + beta[0]*(p->r[-(p->lastpt)]))/(3*alpha[0])
    - ((p->indr)*dln_f(psi,0) - (p->r[-(p->lastpt)]))/sq(psi[0]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl half_outgoing_nullR(const VD& alpha, const VD& beta,
			       const VD& psi, PAR *p, int k)
{
  return (ddr_b(beta,p,k) - (p->r[-k])*beta[k])/(3*alpha[k])
    + ((p->r[-k]) + (p->indr)*dln_b(psi,k))/sq(psi[k]);
    
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_null(const VD& alpha, const VD& beta,
			 const VD& psi, PAR *p, int k)
{
  return (1 + (p->indr)*dln_c(psi,k)/(p->r[-k]))/sq(psi[k])
    + (ddr_c(beta,p,k)/(p->r[-k]) - beta[k])/(3*alpha[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_null_rev(const VD& alpha, const VD& beta,
			 const VD& psi, PAR *p, int k)
{
  return (ddr_c(beta,p,k)/(p->r[-k]) - beta[k])/(3*alpha[k])
    - (1 + (p->indr)*dln_c(psi,k)/(p->r[-k]))/sq(psi[k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_null0(const VD& alpha, const VD& beta,
			  const VD& psi, PAR *p)
{
  return ((p->indr)*dln_f(psi,0)/(p->r[-(p->lastpt)]) - 1)/sq(psi[0])
    - (ddr_f(beta,p,0)/(p->r[-(p->lastpt)]) + beta[0])/(3*alpha[0]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl outgoing_nullR(const VD& alpha, const VD& beta,
			  const VD& psi, PAR *p, int k)
{
  return (1 + (p->indr)*dln_b(psi,k)/(p->r[-k]))/sq(psi[k])
    + (ddr_b(beta,p,k)/(p->r[-k]) - beta[k])/(3*alpha[k]);
    
}
inline void get_outnull(VD& outnull, const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p) {
  outnull[0] = outgoing_null0(f_al, f_be, f_ps, p);
  for (int k = 1; k < p->zerowr; ++k) {
    outnull[k] = outgoing_null_rev(f_al, f_be, f_ps, p, (p->inds[k]).second);
  }
  for (int k = (p->zerowr) + 1; k < p->lastwr; ++k) {
    outnull[k] = outgoing_null(f_al, f_be, f_ps, p, (p->inds[k]).second);
  }
  outnull[p->lastwr] = outgoing_nullR(f_al, f_be, f_ps, p, p->lastpt);
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl sRicci(const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2, const VD& f_ps, int k)
{
  return (sq(f_xi[k]) - sq(f_pi[k]) - sq(f_xi2[k]) + sq(f_pi2[k])) / pw4(f_ps[k]);
}
void get_ricci(VD& ricci, const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
	       const VD& f_ps, const vector< pair<int,int> >& indices)
{
  for (auto k : indices) { ricci[k.first] = sRicci(f_xi, f_pi, f_xi2, f_pi2, f_ps, k.second); }
}
////////////////////////////////////////////////////////////////////////////////////
// HORIZON SEARCH
int search_for_horizon(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p)
{
  int k = p->lastpt;
  if (half_outgoing_nullR(f_al, f_be, f_ps, p, p->lastpt) <= 0) { return k; }
  while (--k > p->zeropt) {
    if (half_outgoing_null(f_al, f_be, f_ps, p, k) <= 0) { return k; }
  }
  while (--k > 0) {
    if (half_outgoing_null_rev(f_al, f_be, f_ps, p, k) <= 0) { return k; }
  }
  if (half_outgoing_null0(f_al, f_be, f_ps, p) <= 0) { return p->npts; }
  return 0;
}


// EINSTEIN EQUATIONS

dbl iresEEtt(dbl dt2ps, dbl dx2al, dbl dx2be, dbl dx2ps,
	     dbl dxdtbe, dbl dxdtps,
	     dbl dtal, dbl dtbe, dbl dtps,
	     dbl dxal, dbl dxbe, dbl dxps,
	     dbl xi, dbl pi, dbl xi2, dbl pi2,
	     dbl al, dbl be, dbl ps, dbl x, dbl xsq, dbl lsq, dbl four_pi)
{
  return (-4*pow(dtps,2))/pow(ps,2) - (4*dt2ps)/ps + (4*dtps*(al*(dtal + 
be*dxal) - pow(be,2)*pow(ps,3)*(-2*dtps + 2*be*dxps +			
dxbe*ps)))/(pow(al,2)*ps) + ((pow(al,2) -				
pow(be,2)*pow(ps,4))*(pow(al,2)*(-2*dxal*dxps + dx2al*ps) +		
dtal*pow(ps,4)*(2*dtps - dxbe*ps) -					
2*pow(be,2)*pow(ps,3)*(al*pow(dxps,2) + (al*dx2ps - dxal*dxps)*ps) + 
be*pow(ps,3)*(4*al*dtps*dxps - 2*(dtps*dxal + dtal*dxps + 
al*(-2*dxdtps + 3*dxbe*dxps))*ps + (-(al*dx2be) + 
dxal*dxbe)*pow(ps,2)) + al*pow(ps,3)*(-2*pow(dtps,2) + 2*(-dt2ps + 
2*dtps*dxbe + dtbe*dxps)*ps + (-pow(dxbe,2) + 
dxdtbe)*pow(ps,2))))/(pow(al,3)*pow(ps,5)) + 2*(dtbe - (be*(dtal + 
be*dxal))/al + be*(-dxbe + (4*dtps)/ps) + (al*dxal)/pow(ps,4) - 
(2*pow(be,2)*dxps)/ps + (pow(be,3)*pow(ps,3)*(-2*dtps + 2*be*dxps + 
dxbe*ps))/pow(al,2))*((2*dxps)/ps + x/(lsq + xsq)) + ((-pow(al,2) + 
pow(be,2)*pow(ps,4))*(pow(al,2)*(lsq + xsq)*(2*dxal*dxps*(lsq + xsq) 
+ ps*(2*x*dxal + dx2al*(lsq + xsq))) + pow(al,3)*(lsq*ps + 4*(lsq + 
xsq)*(2*x*dxps + dx2ps*(lsq + xsq))) + (-dtal + 
be*dxal)*pow(ps,4)*(lsq + xsq)*((-6*dtps + dxbe*ps)*(lsq + xsq) + 
2*be*(x*ps + 3*dxps*(lsq + xsq))) - 
al*pow(ps,3)*(pow(be,2)*(18*pow(dxps,2)*pow(lsq + xsq,2) + 
pow(ps,2)*(2*lsq + xsq) + 2*ps*(lsq + xsq)*(8*x*dxps + 3*dx2ps*(lsq + 
xsq))) + (lsq + xsq)*(18*pow(dtps,2)*(lsq + xsq) - 2*(-3*dt2ps + 
4*dtps*dxbe + 3*dtbe*dxps)*ps*(lsq + xsq) + pow(ps,2)*(-2*x*dtbe + 
pow(dxbe,2)*(lsq + xsq) - dxdtbe*(lsq + xsq))) + be*(lsq + 
xsq)*(-36*dtps*dxps*(lsq + xsq) + pow(ps,2)*(4*x*dxbe + dx2be*(lsq + 
xsq)) + 2*ps*(-8*x*dtps - 6*dxdtps*(lsq + xsq) + 7*dxbe*dxps*(lsq + 
xsq))))))/(pow(al,3)*pow(ps,5)*pow(lsq + xsq,2))
    + four_pi*(sq(be*xi + al*pi/sq(ps)) + sq(be*pi + al*xi/sq(ps)));
}

dbl iresEEtr(dbl dt2ps, dbl dx2al, dbl dx2be, dbl dx2ps,
	     dbl dxdtbe, dbl dxdtps,
	     dbl dtal, dbl dtbe, dbl dtps,
	     dbl dxal, dbl dxbe, dbl dxps,
	     dbl xi, dbl pi, dbl xi2, dbl pi2,
	     dbl al, dbl be, dbl ps, dbl x, dbl xsq, dbl lsq, dbl four_pi)
{
  return (4*pow(al,2)*dtps*dxal*ps*pow(lsq + xsq,2) + 2*be*(-dtal + 
be*dxal)*pow(ps,5)*(lsq + xsq)*(-2*dtps*(lsq + xsq) + be*(x*ps + 
2*dxps*(lsq + xsq))) + pow(al,3)*(4*(dtps*dxps - dxdtps*ps)*pow(lsq + 
xsq,2) + be*ps*(lsq*ps + 4*(lsq + xsq)*(2*x*dxps + dx2ps*(lsq + 
xsq)))) - al*be*pow(ps,4)*(-2*(lsq + xsq)*(x*dtbe*pow(ps,2) - 
4*pow(dtps,2)*(lsq + xsq) + 2*(-dt2ps + dtbe*dxps)*ps*(lsq + xsq)) + 
pow(be,2)*(8*pow(dxps,2)*pow(lsq + xsq,2) + pow(ps,2)*(2*lsq + xsq) + 
4*ps*(lsq + xsq)*(3*x*dxps + dx2ps*(lsq + xsq))) + 2*be*(lsq + 
xsq)*(x*dxbe*pow(ps,2) - 8*dtps*dxps*(lsq + xsq) + 2*ps*(-3*x*dtps - 
2*dxdtps*(lsq + xsq) + dxbe*dxps*(lsq + 
xsq)))))/(pow(al,3)*pow(ps,2)*pow(lsq + xsq,2))
    + four_pi*(be*(sq(xi) + sq(pi)) + 2*al*xi*pi/sq(ps));
}

dbl iresEErr(dbl dt2ps, dbl dx2al, dbl dx2be, dbl dx2ps,
	     dbl dxdtbe, dbl dxdtps,
	     dbl dtal, dbl dtbe, dbl dtps,
	     dbl dxal, dbl dxbe, dbl dxps,
	     dbl xi, dbl pi, dbl xi2, dbl pi2,
	     dbl al, dbl be, dbl ps, dbl x, dbl xsq, dbl lsq, dbl four_pi)
{
  return (2*pow(al,2)*dxal*ps*(lsq + xsq)*(x*ps + 2*dxps*(lsq + xsq)) + 
pow(al,3)*(-(lsq*pow(ps,2)) + 4*x*dxps*ps*(lsq + xsq) + 
4*pow(dxps,2)*pow(lsq + xsq,2)) + 2*(-dtal + be*dxal)*pow(ps,5)*(lsq 
+ xsq)*(-2*dtps*(lsq + xsq) + be*(x*ps + 2*dxps*(lsq + xsq))) - 
al*pow(ps,4)*(-2*(lsq + xsq)*(x*dtbe*pow(ps,2) - 4*pow(dtps,2)*(lsq + 
xsq) + 2*(-dt2ps + dtbe*dxps)*ps*(lsq + xsq)) + 
pow(be,2)*(8*pow(dxps,2)*pow(lsq + xsq,2) + pow(ps,2)*(2*lsq + xsq) + 
4*ps*(lsq + xsq)*(3*x*dxps + dx2ps*(lsq + xsq))) + 2*be*(lsq + 
xsq)*(x*dxbe*pow(ps,2) - 8*dtps*dxps*(lsq + xsq) + 2*ps*(-3*x*dtps - 
2*dxdtps*(lsq + xsq) + dxbe*dxps*(lsq + 
xsq)))))/(pow(al,3)*pow(ps,2)*pow(lsq + xsq,2))
    + four_pi*(sq(xi) + sq(pi));
}

dbl iresEEthth(dbl dt2ps, dbl dx2al, dbl dx2be, dbl dx2ps,
	       dbl dxdtbe, dbl dxdtps,
	       dbl dtal, dbl dtbe, dbl dtps,
	       dbl dxal, dbl dxbe, dbl dxps,
	       dbl xi, dbl pi, dbl xi2, dbl pi2,
	       dbl al, dbl be, dbl ps, dbl x, dbl xsq, dbl lsq, dbl four_pi)
{
  return lsq/(lsq + xsq) - (2*pow(dxps,2)*(lsq + xsq))/pow(ps,2) + (x*dxal + 
dx2al*(lsq + xsq))/al + (2*(x*dxps + dx2ps*(lsq + xsq)))/ps + ((-dtal 
+ be*dxal)*pow(ps,3)*((-4*dtps + dxbe*ps)*(lsq + xsq) + be*(x*ps + 
4*dxps*(lsq + xsq))))/pow(al,3) - 
(pow(ps,2)*(pow(be,2)*(lsq*pow(ps,2) + 8*pow(dxps,2)*pow(lsq + xsq,2) 
+ 2*ps*(lsq + xsq)*(3*x*dxps + 2*dx2ps*(lsq + xsq))) + (lsq + 
xsq)*(8*pow(dtps,2)*(lsq + xsq) - 2*(-2*dt2ps + 3*dtps*dxbe + 
2*dtbe*dxps)*ps*(lsq + xsq) + pow(ps,2)*(-(x*dtbe) + pow(dxbe,2)*(lsq 
+ xsq) - dxdtbe*(lsq + xsq))) + be*(lsq + xsq)*(-16*dtps*dxps*(lsq + 
xsq) + pow(ps,2)*(2*x*dxbe + dx2be*(lsq + xsq)) + 2*ps*(-3*x*dtps - 
4*dxdtps*(lsq + xsq) + 5*dxbe*dxps*(lsq + xsq)))))/(pow(al,2)*(lsq + 
xsq))
    + four_pi*(lsq + xsq)*(sq(pi) - sq(xi));
}

inline void get_EE_adm(WRS *wr, FLDS *f, PAR *p)
{
  dbl dt2ps; dbl dx2al; dbl dx2be; dbl dx2ps; dbl dxdtbe; dbl dxdtps;
  dbl dtal; dbl dtbe; dbl dtps; dbl dxal; dbl dxbe; dbl dxps;
  dbl xi; dbl pi; /*dbl xi2; dbl pi2;*/ dbl al; dbl be; dbl ps;
  dbl x; dbl xsq; dbl lsq = p->lsq;
  int k;
  for (int j = 1; j < p->lastwr; ++j) {
    k = (p->inds[j]).second;
    x = p->r[k];
    xsq = sq(x);
    dt2ps = sq(p->indt)*(2*(f->Ps[k]) - 5*(f->oldPs[k]) + 4*(f->olderPs[k]) - f->oldestPs[k]);
    dx2al = ddr2_c(f->Al, p, k); dx2be = ddr2_c(f->Be, p, k); dx2ps = ddr2_c(f->Ps, p, k);
    dxdtbe = (p->in2dr)*(p->indt)*(1.5*(f->Be[k+1] - f->Be[k-1])
				   - 2*(f->oldBe[k+1] - f->oldBe[k-1])
				   + 0.5*(f->olderBe[k+1] - f->olderBe[k-1]));
    dxdtps = (p->in2dr)*(p->indt)*(1.5*(f->Ps[k+1] - f->Ps[k-1])
				   - 2*(f->oldPs[k+1] - f->oldPs[k-1])
				   + 0.5*(f->olderPs[k+1] - f->olderPs[k-1]));
    dtal = (p->indt)*(1.5*(f->Al[k]) - 2*(f->oldAl[k]) + 0.5*(f->olderAl[k]));
    dtbe = (p->indt)*(1.5*(f->Be[k]) - 2*(f->oldBe[k]) + 0.5*(f->olderBe[k]));
    dtps = (p->indt)*(1.5*(f->Ps[k]) - 2*(f->oldPs[k]) + 0.5*(f->olderPs[k]));
    dxal = ddr_c(f->Al, p, k); dxbe = ddr_c(f->Be, p, k); dxps = ddr_c(f->Ps, p, k);
    xi = f->Xi[k]; pi = f->Pi[k]; //xi2 = f->Xi2[k]; pi2 = f->Pi2[k];
    al = f->Al[k]; be = f->Be[k]; ps = f->Ps[k];

    (wr->p_EEtt).wr_field[j] = (-4*pow(dtps,2))/pow(ps,2) - (4*dt2ps)/ps + (4*dtps*(al*(dtal + 
be*dxal) - pow(be,2)*pow(ps,3)*(-2*dtps + 2*be*dxps +			
dxbe*ps)))/(pow(al,2)*ps) + ((pow(al,2) -				
pow(be,2)*pow(ps,4))*(pow(al,2)*(-2*dxal*dxps + dx2al*ps) +		
dtal*pow(ps,4)*(2*dtps - dxbe*ps) -					
2*pow(be,2)*pow(ps,3)*(al*pow(dxps,2) + (al*dx2ps - dxal*dxps)*ps) + 
be*pow(ps,3)*(4*al*dtps*dxps - 2*(dtps*dxal + dtal*dxps + 
al*(-2*dxdtps + 3*dxbe*dxps))*ps + (-(al*dx2be) + 
dxal*dxbe)*pow(ps,2)) + al*pow(ps,3)*(-2*pow(dtps,2) + 2*(-dt2ps + 
2*dtps*dxbe + dtbe*dxps)*ps + (-pow(dxbe,2) + 
dxdtbe)*pow(ps,2))))/(pow(al,3)*pow(ps,5)) + 2*(dtbe - (be*(dtal + 
be*dxal))/al + be*(-dxbe + (4*dtps)/ps) + (al*dxal)/pow(ps,4) - 
(2*pow(be,2)*dxps)/ps + (pow(be,3)*pow(ps,3)*(-2*dtps + 2*be*dxps + 
dxbe*ps))/pow(al,2))*((2*dxps)/ps + x/(lsq + xsq)) + ((-pow(al,2) + 
pow(be,2)*pow(ps,4))*(pow(al,2)*(lsq + xsq)*(2*dxal*dxps*(lsq + xsq) 
+ ps*(2*x*dxal + dx2al*(lsq + xsq))) + pow(al,3)*(lsq*ps + 4*(lsq + 
xsq)*(2*x*dxps + dx2ps*(lsq + xsq))) + (-dtal + 
be*dxal)*pow(ps,4)*(lsq + xsq)*((-6*dtps + dxbe*ps)*(lsq + xsq) + 
2*be*(x*ps + 3*dxps*(lsq + xsq))) - 
al*pow(ps,3)*(pow(be,2)*(18*pow(dxps,2)*pow(lsq + xsq,2) + 
pow(ps,2)*(2*lsq + xsq) + 2*ps*(lsq + xsq)*(8*x*dxps + 3*dx2ps*(lsq + 
xsq))) + (lsq + xsq)*(18*pow(dtps,2)*(lsq + xsq) - 2*(-3*dt2ps + 
4*dtps*dxbe + 3*dtbe*dxps)*ps*(lsq + xsq) + pow(ps,2)*(-2*x*dtbe + 
pow(dxbe,2)*(lsq + xsq) - dxdtbe*(lsq + xsq))) + be*(lsq + 
xsq)*(-36*dtps*dxps*(lsq + xsq) + pow(ps,2)*(4*x*dxbe + dx2be*(lsq + 
xsq)) + 2*ps*(-8*x*dtps - 6*dxdtps*(lsq + xsq) + 7*dxbe*dxps*(lsq + 
xsq))))))/(pow(al,3)*pow(ps,5)*pow(lsq + xsq,2))
    + (p->four_pi)*(sq(be*xi + al*pi/sq(ps)) + sq(be*pi + al*xi/sq(ps)));
    
    (wr->p_EEtx).wr_field[j] = (4*pow(al,2)*dtps*dxal*ps*pow(lsq + xsq,2) + 2*be*(-dtal + 
be*dxal)*pow(ps,5)*(lsq + xsq)*(-2*dtps*(lsq + xsq) + be*(x*ps + 
2*dxps*(lsq + xsq))) + pow(al,3)*(4*(dtps*dxps - dxdtps*ps)*pow(lsq + 
xsq,2) + be*ps*(lsq*ps + 4*(lsq + xsq)*(2*x*dxps + dx2ps*(lsq + 
xsq)))) - al*be*pow(ps,4)*(-2*(lsq + xsq)*(x*dtbe*pow(ps,2) - 
4*pow(dtps,2)*(lsq + xsq) + 2*(-dt2ps + dtbe*dxps)*ps*(lsq + xsq)) + 
pow(be,2)*(8*pow(dxps,2)*pow(lsq + xsq,2) + pow(ps,2)*(2*lsq + xsq) + 
4*ps*(lsq + xsq)*(3*x*dxps + dx2ps*(lsq + xsq))) + 2*be*(lsq + 
xsq)*(x*dxbe*pow(ps,2) - 8*dtps*dxps*(lsq + xsq) + 2*ps*(-3*x*dtps - 
2*dxdtps*(lsq + xsq) + dxbe*dxps*(lsq + 
xsq)))))/(pow(al,3)*pow(ps,2)*pow(lsq + xsq,2))
    + (p->four_pi)*(be*(sq(xi) + sq(pi)) + 2*al*xi*pi/sq(ps));
    
    (wr->p_EExx).wr_field[j] = (2*pow(al,2)*dxal*ps*(lsq + xsq)*(x*ps + 2*dxps*(lsq + xsq)) + 
pow(al,3)*(-(lsq*pow(ps,2)) + 4*x*dxps*ps*(lsq + xsq) + 
4*pow(dxps,2)*pow(lsq + xsq,2)) + 2*(-dtal + be*dxal)*pow(ps,5)*(lsq 
+ xsq)*(-2*dtps*(lsq + xsq) + be*(x*ps + 2*dxps*(lsq + xsq))) - 
al*pow(ps,4)*(-2*(lsq + xsq)*(x*dtbe*pow(ps,2) - 4*pow(dtps,2)*(lsq + 
xsq) + 2*(-dt2ps + dtbe*dxps)*ps*(lsq + xsq)) + 
pow(be,2)*(8*pow(dxps,2)*pow(lsq + xsq,2) + pow(ps,2)*(2*lsq + xsq) + 
4*ps*(lsq + xsq)*(3*x*dxps + dx2ps*(lsq + xsq))) + 2*be*(lsq + 
xsq)*(x*dxbe*pow(ps,2) - 8*dtps*dxps*(lsq + xsq) + 2*ps*(-3*x*dtps - 
2*dxdtps*(lsq + xsq) + dxbe*dxps*(lsq + 
xsq)))))/(pow(al,3)*pow(ps,2)*pow(lsq + xsq,2))
    + (p->four_pi)*(sq(xi) + sq(pi));
    
    (wr->p_EEhh).wr_field[j] = (lsq/(lsq + xsq) - (2*pow(dxps,2)*(lsq + xsq))/pow(ps,2) + (x*dxal + 
dx2al*(lsq + xsq))/al + (2*(x*dxps + dx2ps*(lsq + xsq)))/ps + ((-dtal 
+ be*dxal)*pow(ps,3)*((-4*dtps + dxbe*ps)*(lsq + xsq) + be*(x*ps + 
4*dxps*(lsq + xsq))))/pow(al,3) - 
(pow(ps,2)*(pow(be,2)*(lsq*pow(ps,2) + 8*pow(dxps,2)*pow(lsq + xsq,2) 
+ 2*ps*(lsq + xsq)*(3*x*dxps + 2*dx2ps*(lsq + xsq))) + (lsq + 
xsq)*(8*pow(dtps,2)*(lsq + xsq) - 2*(-2*dt2ps + 3*dtps*dxbe + 
2*dtbe*dxps)*ps*(lsq + xsq) + pow(ps,2)*(-(x*dtbe) + pow(dxbe,2)*(lsq 
+ xsq) - dxdtbe*(lsq + xsq))) + be*(lsq + xsq)*(-16*dtps*dxps*(lsq + 
xsq) + pow(ps,2)*(2*x*dxbe + dx2be*(lsq + xsq)) + 2*ps*(-3*x*dtps - 
4*dxdtps*(lsq + xsq) + 5*dxbe*dxps*(lsq + xsq)))))/(pow(al,2)*(lsq + 
							       xsq)))/(lsq + xsq)
    + (p->four_pi)*(sq(pi) - sq(xi));

    /*
    wr->EEtt[j] = iresEEtt(dt2ps, dx2al, dx2be, dx2ps, dxdtbe, dxdtps,
			   dtal, dtbe, dtps, dxal, dxbe, dxps,
			   xi, pi, xi2, pi2, al, be, ps, x, xsq, p->lsq, p->four_pi);
    wr->EEtx[j] = iresEEtr(dt2ps, dx2al, dx2be, dx2ps, dxdtbe, dxdtps,
			   dtal, dtbe, dtps, dxal, dxbe, dxps,
			   xi, pi, xi2, pi2, al, be, ps, x, xsq, p->lsq, p->four_pi);
    wr->EExx[j] = iresEErr(dt2ps, dx2al, dx2be, dx2ps, dxdtbe, dxdtps,
			   dtal, dtbe, dtps, dxal, dxbe, dxps,
			   xi, pi, xi2, pi2, al, be, ps, x, xsq, p->lsq, p->four_pi);
    wr->EEhh[j] = iresEEthth(dt2ps, dx2al, dx2be, dx2ps, dxdtbe, dxdtps,
			     dtal, dtbe, dtps, dxal, dxbe, dxps,
			     xi, pi, xi2, pi2, al, be, ps, x, xsq, p->lsq, p->four_pi);
    */

    (wr->p_hamiltonian).wr_field[j] =
      ((p->eight_pi)*ps*(pow(pi,2) + pow(xi,2)) + (-2*pow(al,2)*(lsq*ps + 4*(lsq +
xsq)*(2*x*dxps + dx2ps*(lsq + xsq))) + 2*pow(ps,3)*(-2*dtps*(lsq +
xsq) + be*(x*ps + 2*dxps*(lsq + xsq)))*(2*(-3*dtps + dxbe*ps)*(lsq +
xsq) + be*(x*ps + 6*dxps*(lsq + xsq))))/(pow(al,2)*pow(lsq +
xsq,2)))/pow(ps,5);

    (wr->p_momentum).wr_field[j] =
      (-(p->eight_pi)*pi*xi + (2*dxal*ps*(lsq + xsq)*(-2*dtps*(lsq + xsq) +
be*(x*ps + 2*dxps*(lsq + xsq))) - 2*al*(2*(dtps*dxps -
dxdtps*ps)*pow(lsq + xsq,2) + be*(lsq*pow(ps,2) -
2*pow(dxps,2)*pow(lsq + xsq,2) + 2*ps*(lsq + xsq)*(x*dxps +
dx2ps*(lsq + xsq)))))/(pow(al,2)*pow(lsq + xsq,2)))/pow(ps,2);

    (wr->p_kext).wr_field[j] =
      ((-6*dtps + dxbe*ps)*(lsq + xsq) + 2*be*(x*ps + 3*dxps*(lsq + xsq)))
      /(al*ps*(lsq + xsq));

    (wr->p_dtkext).wr_field[j] =
      -((pow(al,2)*(lsq + xsq)*(2*dxal*dxps*(lsq + xsq) + ps*(2*x*dxal +
dx2al*(lsq + xsq))) + al*ps*((p->eight_pi)*pow(al,2)*pow(pi,2)*pow(lsq +
xsq,2) - pow(ps,2)*pow(-2*dtps + 2*be*dxps + dxbe*ps,2)*pow(lsq +
xsq,2) - 2*pow(ps,2)*pow(-2*dtps*(lsq + xsq) + be*(x*ps + 2*dxps*(lsq
+ xsq)),2)) + be*pow(ps,3)*(2*x*al*ps*((-6*dtps + dxbe*ps)*(lsq +
xsq) + 2*be*(x*ps + 3*dxps*(lsq + xsq))) + al*dxps*(lsq +
xsq)*((-6*dtps + dxbe*ps)*(lsq + xsq) + 2*be*(x*ps + 3*dxps*(lsq +
xsq))) + dxal*ps*(lsq + xsq)*((-6*dtps + dxbe*ps)*(lsq + xsq) +
2*be*(x*ps + 3*dxps*(lsq + xsq))) - al*ps*(lsq + xsq)*(2*x*(-6*dtps +
dxbe*ps) + (-6*dxdtps + dxbe*dxps + dx2be*ps)*(lsq + xsq) +
2*be*(7*x*dxps + ps + 3*dx2ps*(lsq + xsq)) + 2*dxbe*(x*ps +
3*dxps*(lsq + xsq)))))/(pow(al,2)*pow(ps,5)*pow(lsq + xsq,2)));
  }
  return;
}

#endif
