#ifndef JACOBIAN_H_INCLUDED
#define JACOBIAN_H_INCLUDED

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
#include "sim-header.h"
#include "sim-structs.h"
#include "fda-io.h"
#include "fda-fns.h"
#include "ellis-fns.h"
#include "ellis-proc.h"


using namespace std;

//  for LAPACKE_dgbsv(): jac[ (kl + ku + 1) + (ldab - 1)*j + i ]  =  jac[ i, j ]
void set_jacCM_abp(VD& jac, const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		   const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, MAPID& r,
		   int npts, int kl, int ku, int ldab);
void set_jacCM_ab(VD& jac, const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		  const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, MAPID& r,
		  int npts, int kl, int ku, int ldab);
inline int jac_ind(int i, int j) { return (4 + i + 6*j); } // (kl + ku + i + (2*kl + ku)*j); } 
// ***********************  JACOBIAN FUNCTIONS  ***********************
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl jac_aa(const VD& f_pi, const VD& f_pi2, const VD& f_al, const VD& f_be, const VD& f_ps,
		  PAR *p, int k)
{
  return (p->neg2indrsq) + (p->eight_pi)*(sq(f_pi[k]) - sq(f_pi2[k])) +
    (p->two_thirds)*pw4(f_ps[k])*sq(ddr_c(f_be,p,k) - f_be[k]*(p->r[-k])) / sq(f_al[k]);
}

inline dbl jac_aa_pm(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k, int p_m)
{
  return (p->indrsq) + p_m*(p->indr)*(ddrln_c(f_ps,p,k) + (p->r[-k]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl jac_bb(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k)
{
  return (p->neg2indrsq) - 1/r2(p,k)
    - (p->r[-k])*(p->r[-k] + 6*ddrln_c(f_ps,p,k) - ddrln_c(f_al,p,k));
}

inline dbl jac_bb_pm(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k, int p_m)
{
  return (p->indrsq) + p_m*(p->in2dr)*(2*(p->r[-k]) + 6*ddrln_c(f_ps,p,k) - ddrln_c(f_al,p,k));
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl jac_pp(const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		  const VD& f_al, const VD& f_be, const VD& f_ps,
		  PAR *p, int k)
{
  return (p->neg2indrsq) + M_PI*(sq(f_xi2[k]) + sq(f_pi2[k]) - sq(f_xi[k]) - sq(f_pi[k]))
    + 0.25*(p->lsq)/sq(r2(p,k)) +
    (p->five_twelfths)*pw4(f_ps[k])*sq(ddr_c(f_be,p,k) - f_be[k]*(p->r[-k])) / sq(f_al[k]);
}

inline dbl jac_pp_pm(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k, int p_m)
{
  return (p->indrsq) + p_m*(p->indr)*(p->r[-k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// **********************************************************
// **********************************************************
//                    POPULATING JACOBIAN
// **********************************************************
// **********************************************************
////////////////////////////////////////////////////////////////////////////////////////////////


void set_jacCM_abp(VD& jac, const VD& f_xi, const VD& f_pi, const VD& f_xi2, const VD& f_pi2,
		   const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p)
{
  int j = 0;
  int jbe = (p->npts);
  int jps = jbe + jbe;
  // ROW 0, COL 0
  jac[jac_ind(j,j)] = -(p->jacRR);
  jac[jac_ind(jbe,jbe)] = -(p->jacRR);
  jac[jac_ind(jps,jps)] = -(p->jacRR);
  // ROW 0, COL 1
  jac[jac_ind(j,j + 1)] = -(p->jacRRm1);
  jac[jac_ind(jbe,jbe + 1)] = -(p->jacRRm1);
  jac[jac_ind(jps,jps + 1)] = -(p->jacRRm1);
  // ROW 0, COL 2
  jac[jac_ind(j,j + 2)] = -(p->jacRRm2);
  jac[jac_ind(jbe,jbe + 2)] = -(p->jacRRm2);
  jac[jac_ind(jps,jps + 2)] = -(p->jacRRm2);  
  
  for (j = 1; j < (p->lastpt); ++j) {
    ++jbe;
    ++jps;
    // ROW j, COL j-1
    jac[jac_ind(j,j - 1)] = jac_aa_pm(f_al, f_be, f_ps, p, j, -1);
    jac[jac_ind(jbe,jbe - 1)] = jac_bb_pm(f_al, f_be, f_ps, p, j, -1);
    jac[jac_ind(jps,jps - 1)] = jac_pp_pm(f_al, f_be, f_ps, p, j, -1);
    // ROW j, COL j
    jac[jac_ind(j,j)] = jac_aa(f_pi, f_pi2, f_al, f_be, f_ps, p, j);
    jac[jac_ind(jbe,jbe)] = jac_bb(f_al, f_be, f_ps, p, j);
    jac[jac_ind(jps,jps)] = jac_pp(f_xi, f_pi, f_xi2, f_pi2, f_al, f_be, f_ps, p, j);
    // ROW j+1, COL j
    jac[jac_ind(j,j + 1)] = jac_aa_pm(f_al, f_be, f_ps, p, j, 1);
    jac[jac_ind(jbe,jbe + 1)] = jac_bb_pm(f_al, f_be, f_ps, p, j, 1);
    jac[jac_ind(jps,jps + 1)] = jac_pp_pm(f_al, f_be, f_ps, p, j, 1);
  }
  j = (p->lastpt);
  ++jbe;
  ++jps;
  // ROW lastpt, COL lastpt-2
  jac[jac_ind(j,j - 2)] = (p->jacRRm2);
  jac[jac_ind(jbe,jbe - 2)] = (p->jacRRm2);
  jac[jac_ind(jps,jps - 2)] = (p->jacRRm2);
  // ROW lastpt, COL lastpt-1
  jac[jac_ind(j,j - 1)] = (p->jacRRm1);
  jac[jac_ind(jbe,jbe - 1)] = (p->jacRRm1);
  jac[jac_ind(jps,jps - 1)] = (p->jacRRm1);
  // ROW lastpt, COL lastpt
  jac[jac_ind(j,j)] = (p->jacRR);
  jac[jac_ind(jbe,jbe)] = (p->jacRR);
  jac[jac_ind(jps,jps)] = (p->jacRR);
  return;
}

void set_jacCM_ab(VD& jac, const VD& f_pi, const VD& f_pi2,
		  const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p)
{
  int j = 0;
  int jbe = (p->npts);
  // ROW 0, COL 0
  jac[jac_ind(j,j)] = -(p->jacRR);
  jac[jac_ind(jbe,jbe)] = -(p->jacRR);
  // ROW 0, COL 1
  jac[jac_ind(j,j + 1)] = -(p->jacRRm1);
  jac[jac_ind(jbe,jbe + 1)] = -(p->jacRRm1);
  // ROW 0, COL 2
  jac[jac_ind(j,j + 2)] = -(p->jacRRm2);
  jac[jac_ind(jbe,jbe + 2)] = -(p->jacRRm2);
  for (j = 1; j < (p->lastpt); ++j) {
    ++jbe;
    // ROW j, COL j-1
    jac[jac_ind(j,j - 1)] = jac_aa_pm(f_al, f_be, f_ps, p, j, -1);
    jac[jac_ind(jbe,jbe - 1)] = jac_bb_pm(f_al, f_be, f_ps, p, j, -1);
    // ROW j, COL j
    jac[jac_ind(j,j)] = jac_aa(f_pi, f_pi2, f_al, f_be, f_ps, p, j);
    jac[jac_ind(jbe,jbe)] = jac_bb(f_al, f_be, f_ps, p, j);
    // ROW j+1, COL j
    jac[jac_ind(j,j + 1)] = jac_aa_pm(f_al, f_be, f_ps, p, j, 1);
    jac[jac_ind(jbe,jbe + 1)] = jac_bb_pm(f_al, f_be, f_ps, p, j, 1);
  }
  j = (p->lastpt);
  ++jbe;
  // ROW lastpt, COL lastpt-2
  jac[jac_ind(j,j - 2)] = (p->jacRRm2);
  jac[jac_ind(jbe,jbe - 2)] = (p->jacRRm2);
  // ROW lastpt, COL lastpt-1
  jac[jac_ind(j,j - 1)] = (p->jacRRm1);
  jac[jac_ind(jbe,jbe - 1)] = (p->jacRRm1);
  // ROW lastpt, COL lastpt
  jac[jac_ind(j,j)] = (p->jacRR);
  jac[jac_ind(jbe,jbe)] = (p->jacRR);
  return;
}

#endif
