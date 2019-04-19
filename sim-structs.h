#ifndef SIM_STRUCTS_H_INCLUDED
#define SIM_STRUCTS_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include "sim-header.h"

/*
typedef struct bbhutil_params BBHP;
typedef struct sim_fields FLDS;
typedef struct sim_writers WRS;
typedef struct sim_params PAR;
typedef struct sim_site_values SSV;
*/

struct bbhutil_params {
  str filename = "";
  VD *full_field = NULL;
  VD wr_field;
  char *file = NULL;
  int *shape = NULL;
  int rank = 1;
  dbl *coords = NULL;
  dbl *data = NULL;
  bool write_this = false;
  D_FN compute = NULL;
} ;


struct sim_fields {
  VD Al;
  VD oldAl;
  VD cnAl;
  //VD resAl;
  VD olderAl;
  VD Be;
  VD oldBe;
  VD cnBe;
  //VD resBe;
  VD olderBe; 
  VD Ps;
  VD oldPs;
  VD cnPs;
  VD resPs;
  VD olderPs;
  VD oldestPs;
  VD Xi;
  VD oldXi;
  VD cnXi;
  VD resXi;
  VD olderXi;
  VD Pi;
  VD oldPi;
  VD cnPi;
  VD resPi;
  VD olderPi;
  VD Xi2;
  VD oldXi2;
  VD cnXi2;
  VD resXi2;
  VD olderXi2;
  VD Pi2;
  VD oldPi2;
  VD cnPi2;
  VD resPi2;
  VD olderPi2;
  VD res_ell;
  VD jac;
  VD areal;
} ;


struct sim_writers {
  BBHP p_Al;
  BBHP p_resAl;
  BBHP p_iresAl;
  BBHP p_Be;
  BBHP p_resBe;
  BBHP p_iresBe;
  BBHP p_Ps;
  BBHP p_resPs;
  BBHP p_iresPs;
  BBHP p_Xi;
  BBHP p_resXi;
  BBHP p_iresXi;
  BBHP p_Pi;
  BBHP p_resPi;
  BBHP p_iresPi;
  BBHP p_Xi2;
  BBHP p_resXi2;
  BBHP p_iresXi2;
  BBHP p_Pi2;
  BBHP p_resPi2;
  BBHP p_iresPi2;
  BBHP p_maspect;
  BBHP p_outnull;
  BBHP p_revnull;
  BBHP p_areal;
  BBHP p_ricci;
  BBHP p_EEtt;
  BBHP p_EEtx;
  BBHP p_EExx;
  BBHP p_EEhh;
  BBHP p_momentum;
  BBHP p_hamiltonian;
  BBHP p_kext;
  BBHP p_dtkext;
  SSVSET get_site_vals = NULL;
  SSVSET get_bound_vals = NULL;
  void (*set_old_fields)(FLDS *) = NULL;
  BBHP p_amin;
  BBHP p_xamin;
} ;


struct sim_params {
  str outfile = "ellis";
  int lastpt = 1000;
  int save_pt = 1;
  int nsteps = 2000;
  int save_step = 1;
  int resn_factor = 1;
  int check_step = 1;
  int maxit = 100;
  int norm_type = 0; // type of l-norm to compute
  int signal_code = 0; // code for multipurpose signaling
  dbl lam = 0.25; // dt/dr
  dbl rmax = 100;
  dbl rmin = -rmax;
  dbl lsq = 1;
  dbl dspn = 0.5; // dissipation coefficient
  dbl tol = 0.0000000001; // iterative method tolerance
  dbl ell_tol = 0.01*tol;
  dbl ell_up_weight = 0.5;
  dbl ic_Dsq = 25.0; // gaussian width
  dbl ic_r0 = 20.0; // gaussian center
  dbl ic_Amp = 0; // gaussian amplitude
  dbl ic2_Dsq = 25.0; // gaussian width
  dbl ic2_r0 = 20.0; // gaussian center
  dbl ic2_Amp = 0.001; // gaussian amplitude
  bool psi_hyp = true; // update psi with hyperbolic evolution eqn after IC?
  bool somm_cond = true; // sommerfeld condition at outer bound?
  bool dspn_bound = false; // dissipate boundary points?
  bool dspn_psi = false; // dissipate psi (only activated if psi_hyp=true)?
  bool static_metric = false; // ignore dynamics of metric?
  bool clean_hyp = true; // zero initial data for pi?
  bool clean_ell = false; // use clean hyperbolic update functions (slower)?
  bool write_res = false; // write residuals?
  bool write_ricci = false; // write ricci?
  bool write_itn = false; // write iteration counts?
  bool write_mtot = true; // total diagnostics?
  bool write_maspect = false; // write mass aspect?
  bool write_outnull = false; // write outgoing null expansion?
  bool write_xp = true; // write ghost xi and pi?
  bool write_xp2 = false; // write normal xi and pi?
  bool write_abp = true; // write metric fields (alpha, beta, psi)?
  bool write_ires_xp = false; // write ires for xi and pi?
  bool write_ires_xp2 = false;
  bool write_ires_abp = false; // write ires for metric variables?
  bool horizon_search = false; // search for apparent horizon after each step?
  bool sym_pert = false; // symmetric perturbation?
  // for convergence purposes
  str hold_const = "lambda"; // "lambda", "dt", or "dr"
  bool same_times = true;
  bool same_grids = true;

  int n_ell = 3;

  int lastwr = lastpt/save_pt;
  int wr_shape = lastwr + 1;
  int zerowr = lastwr/2;
  dbl wr_dr = (rmax - rmin) / ((dbl) lastwr);
  dbl coord_lims[2] = {rmin, rmax};
  vector< pair<int,int> > inds;
  int check_diagnostics = save_step * check_step;
  
  MAPID r {{0, rmin}};
  dbl t = 0;
  int exit_itn = maxit;

  // DERIVED PARAMS
  int npts = lastpt + 1;
  int zeropt = lastpt/2;
  dbl norm_factor = 1 / ((dbl) npts);
  dbl dr = (rmax - rmin) / ((dbl) lastpt);
  dbl dt = lam * dr;

  // LAPACK PARAMS
  int lp_n = n_ell * npts;
  int lp_kl = 2;
  int lp_ku = 2;
  int lp_nrhs = 1;
  int lp_ldab = 2*lp_kl + lp_ku + 1;
  int lp_ldb = lp_n;
  int *lp_ipiv = NULL;
  vector<int> ipiv; 
  
  // FREQUENTLY USED
  dbl one_third = 1.0 / 3.0;
  dbl two_thirds = 2 * one_third;
  //dbl four_thirds = 2 * two_thirds;
  dbl twelfth = 0.25 * one_third;
  dbl five_twelfths = 5 * twelfth;
  dbl four_pi = 4 * M_PI;
  dbl eight_pi = 2 * four_pi;
  dbl twelve_pi = 3 * four_pi;
  dbl lam2val = 0.5 * lam;
  dbl lam6val = lam2val * one_third;
  dbl drsq = dr * dr;
  dbl indr = 1 / dr;
  dbl in2dr = 0.5 * indr;
  dbl in3dr = one_third * indr;
  dbl indrsq = indr * indr;
  dbl neg2indrsq = -2 * indrsq;
  dbl indt = 1 / dt;
  dbl inrmax = 1 / rmax;
  // SPECIFIC TERMS
  dbl csomm_old = 1 - 0.75*lam - 0.5*dt*inrmax;
  dbl csomm_rhs = 1 / (2 - csomm_old);
  dbl jacRR = 3*in2dr + inrmax;
  dbl jacRRm1 = -4 * in2dr;
  dbl jacRRm2 = in2dr;
  dbl cpsi_rhs = 1 / jacRR;

  int (*solver)(FLDS *, PAR *) = NULL;
  
} ;

  

struct sim_site_values {
  dbl lsq = 1;
  dbl x; dbl xsq; dbl rsq;
  dbl dt2ps;
  dbl dx2al; dbl dx2be; dbl dx2ps;
  dbl dxdtbe; dbl dxdtps;
  dbl dtxi; dbl dtpi; dbl dtxi2; dbl dtpi2; 
  dbl dtal; dbl dtbe; dbl dtps;
  dbl dxxi; dbl dxpi; dbl dxxi2; dbl dxpi2;
  dbl dxal; dbl dxbe; dbl dxps;
  dbl xi; dbl pi; dbl xi2; dbl pi2;
  dbl al; dbl be; dbl ps;
} ;



#endif
