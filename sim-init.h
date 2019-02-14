#ifndef SIM_INIT_H_INCLUDED
#define SIM_INIT_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <cmath> // for ICs
#include <vector> // for everything
#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "lapacke.h"
#include "fda-io.h"
#include "fda-fns.h"
#include "ellis-fns.h"
#include "jacobian.h"
//#include "ellis-clean.h"
#include "ellis-proc.h"
#include "solvers.h"
#include "sim-header.h"
#include "sim-structs.h"

void bbhp_init(BBHP *bp, PAR *p, str fieldname, VD *p_field, const VD& zeros)
{
  bp->full_field = p_field;
  bp->wr_field = zeros;
  bp->filename = fieldname + "-" + p->outfile + ".sdf";
  bp->file = &(bp->filename[0]);
  bp->shape = &(p->wr_shape);
  bp->rank = 1;
  bp->coords = &(p->coord_lims[0]);  
  bp->data = &(bp->wr_field[0]);
  bp->write_this = true;
  return;
}

vector<BBHP *> writers_init(WRS *wr, FLDS *f, PAR *p)
{
  VD zeros(p->wr_shape, 0);
  vector<BBHP *> out_vec;
  if (p->write_abp) {
    bbhp_init(&(wr->p_Al), p, "Al", &(f->Al), zeros);
    bbhp_init(&(wr->p_Be), p, "Be", &(f->Be), zeros);
    bbhp_init(&(wr->p_Ps), p, "Ps", &(f->Ps), zeros);
    out_vec.push_back(&(wr->p_Al));
    out_vec.push_back(&(wr->p_Be));
    out_vec.push_back(&(wr->p_Ps));
    if (p->write_res) {
      bbhp_init(&(wr->p_resAl), p, "ResAl", &(f->resAl), zeros);
      bbhp_init(&(wr->p_resBe), p, "ResBe", &(f->resBe), zeros);
      bbhp_init(&(wr->p_resPs), p, "ResPs", &(f->resPs), zeros);
      out_vec.push_back(&(wr->p_resAl));
      out_vec.push_back(&(wr->p_resBe));
      out_vec.push_back(&(wr->p_resPs));
    }
  }
  if (p->write_xp) {
    bbhp_init(&(wr->p_Xi), p, "Xi", &(f->Xi), zeros);
    bbhp_init(&(wr->p_Pi), p, "Pi", &(f->Pi), zeros);
    out_vec.push_back(&(wr->p_Xi));
    out_vec.push_back(&(wr->p_Pi));
    if (p->write_res) {
      bbhp_init(&(wr->p_resXi), p, "ResXi", &(f->resXi), zeros);
      bbhp_init(&(wr->p_resPi), p, "ResPi", &(f->resPi), zeros);
      out_vec.push_back(&(wr->p_resXi));
      out_vec.push_back(&(wr->p_resPi));
    }
  }
  if (p->write_xp2) {
    bbhp_init(&(wr->p_Xi2), p, "Xi2", &(f->Xi2), zeros);
    bbhp_init(&(wr->p_Pi2), p, "Pi2", &(f->Pi2), zeros);
    out_vec.push_back(&(wr->p_Xi2));
    out_vec.push_back(&(wr->p_Pi2));
    if (p->write_res) {
      bbhp_init(&(wr->p_resXi2), p, "ResXi2", &(f->resXi2), zeros);
      bbhp_init(&(wr->p_resPi2), p, "ResPi2", &(f->resPi2), zeros);
      out_vec.push_back(&(wr->p_resXi2));
      out_vec.push_back(&(wr->p_resPi2));
    }
  }
  if (p->write_ires_abp) {
    bbhp_init(&(wr->p_iresAl), p, "iresAl", NULL, zeros);
    bbhp_init(&(wr->p_iresBe), p, "iresBe", NULL, zeros);
    bbhp_init(&(wr->p_iresPs), p, "iresPs", NULL, zeros);
  }
  if (p->write_ires_xp) {
    bbhp_init(&(wr->p_iresXi), p, "iresXi", NULL, zeros);
    bbhp_init(&(wr->p_iresPi), p, "iresPi", NULL, zeros);
  }
  if (p->write_ires_xp2) {
    bbhp_init(&(wr->p_iresXi2), p, "iresXi2", NULL, zeros);
    bbhp_init(&(wr->p_iresPi2), p, "iresPi2", NULL, zeros);
  }
  if (p->write_maspect) { bbhp_init(&(wr->p_maspect), p, "maspect", NULL, zeros); }
  if (p->write_outnull) { bbhp_init(&(wr->p_outnull), p, "outnull", NULL, zeros); }
  if (p->write_ricci) { bbhp_init(&(wr->p_ricci), p, "ricci", NULL, zeros); }
  return out_vec;
}

int fields_init(FLDS *f, PAR *p)
{
  VD zeros(p->npts, 0);
  f->Xi = zeros;
  f->Pi = zeros;
  f->Xi2 = zeros;
  f->Pi2 = zeros;
  f->Al = zeros;
  f->Be = zeros;
  f->Ps = zeros;
  // PUT INITIAL CONDITIONS
  for (int k = 0; k < (p->zeropt) + 1; ++k) {
    f->Al[k] = 1;
    //f->Be[k] = 0;
    f->Ps[k] = 1;
  }
  for (int k = (p->zeropt) + 1; k < (p->npts); ++k) {
    f->Al[k] = 1;
    //f->Be[k] = 0;
    f->Ps[k] = 1;
    f->Xi[k] = ic_xi(p->r[k], p->ic_Amp, p->ic_Dsq, p->ic_r0);
    f->Pi[k] = ic_pi(p->r[k], p->ic_Amp, p->ic_Dsq, p->ic_r0);
    f->Xi2[k] = ic_xi(p->r[k], p->ic2_Amp, p->ic2_Dsq, p->ic2_r0);
    f->Pi2[k] = ic_pi(p->r[k], p->ic2_Amp, p->ic2_Dsq, p->ic2_r0);
  }

  f->oldXi = f->Xi;
  f->cnXi = f->Xi;
  f->resXi = zeros;
  f->oldPi = f->Pi;
  f->cnPi = f->Pi;
  f->resPi = zeros;
  if (p->write_ires_xp) {
    f->olderXi = zeros;
    f->olderPi = zeros;
  }
  f->oldXi2 = f->Xi2;
  f->cnXi2 = f->Xi2;
  f->resXi2 = zeros;
  f->oldPi2 = f->Pi2;
  f->cnPi2 = f->Pi2;
  f->resPi2 = zeros;
  if (p->write_ires_xp2) {
    f->olderXi2 = zeros;
    f->olderPi2 = zeros;
  }

  if (!(p->static_metric)) {
    int t0_itn = solve_t0(f, p);
    if (t0_itn < 0) {
      if (t0_itn > -(p->npts)) {
	record_horizon(p, f->Ps, -t0_itn, 0, 0);
      }
      else if (t0_itn == -(p->npts)) {
	record_horizon(p, f->Ps, 0, 0, 0);
      }
      else { cout << "\nUNDETERMINED ERROR IN T = 0 ELLIPTIC CONSTRAINTS" << endl; }
      return t0_itn;
    }
    else if (t0_itn == 2*(p->maxit)) {
      cout << "\nUNDETERMINED ERROR IN T = 0 ELLIPTIC CONSTRAINTS, MAXIT REACHED" << endl;
      return t0_itn;
    }
  }
  f->oldAl = f->Al;
  f->cnAl = f->Al;
  f->oldBe = f->Be;
  f->cnBe = f->Be;
  
  f->oldPs = f->Ps;
  f->cnPs = f->Ps;
  if (p->psi_hyp) { f->resPs = zeros; }
  if (p->write_ires_abp) {
    //f->olderAl = zeros;
    //f->olderBe = zeros;
    f->olderPs = zeros;
  }
  VD ell_res_zeros((p->lp_ldb), 0);
  VD jac_zeros(((p->lp_ldab)*(p->lp_n)), 0);
  f->res_ell = ell_res_zeros;
  f->jac = jac_zeros;
  
  return 0;
}

int params_init(PAR *p, int argc, char **argv)
{
  map<str, str *> p_str {{"-outfile",&(p->outfile)}};
  map<str, int *> p_int {{"-lastpt",&(p->lastpt)}, {"-save_pt",&(p->save_pt)},
      {"-nsteps",&(p->nsteps)}, {"-save_step",&(p->save_step)},
      {"-norm_type",&(p->norm_type)}, {"-maxit",&(p->maxit)},
      {"-check_step",&(p->check_step)}, {"-resn_factor",&(p->resn_factor)}};
  map<str, dbl *> p_dbl {{"-lam",&(p->lam)}, {"-lsq",&(p->lsq)}, {"-rmin",&(p->rmin)}, {"-rmax",&(p->rmax)},
      {"-dspn",&(p->dspn)}, {"-tol",&(p->tol)}, {"-ell_tol",&(p->ell_tol)}, {"-ell_up_weight",&(p->ell_up_weight)},
      {"-ic_Dsq",&(p->ic_Dsq)}, {"-ic_r0",&(p->ic_r0)}, {"-ic_Amp",&(p->ic_Amp)},
      {"-ic2_Dsq",&(p->ic2_Dsq)}, {"-ic2_r0",&(p->ic2_r0)}, {"-ic2_Amp",&(p->ic2_Amp)}};
  map<str, bool *> p_bool { {"-psi_hyp",&(p->psi_hyp)}, {"-static_metric",&(p->static_metric)},
      {"-somm_cond",&(p->somm_cond)}, {"-dspn_bound",&(p->dspn_bound)}, {"-dr3_up",&(p->dr3_up)}, {"-dspn_psi",&(p->dspn_psi)},
      {"-write_res",&(p->write_res)},{"-write_ricci",&(p->write_ricci)}, {"-write_itn",&(p->write_itn)},
      {"-write_mtot",&(p->write_mtot)},{"-write_maspect",&(p->write_maspect)}, {"-write_outnull",&(p->write_outnull)},
      {"-write_xp",&(p->write_xp)}, {"-write_xp2",&(p->write_xp2)}, {"-write_abp",&(p->write_abp)},
      {"-write_ires_xp",&(p->write_ires_xp)}, {"-write_ires_xp2",&(p->write_ires_xp2)}, {"-write_ires_abp",&(p->write_ires_abp)},
      {"-clean_hyp",&(p->clean_hyp)}, {"-clean_ell",&(p->clean_ell)}, {"-horizon_search",&(p->horizon_search)}};
  map<str, str> params;
  if (argc > 1) { param_collect(argv, argc, params); }
  else { file_param_collect("ellis-parameters.txt", params); }
  param_set(params, p_str, p_int, p_dbl, p_bool);
  
  // *****************************************CHECKS**********************
  // check that grid size (lastpt = npts-1) is divisible by save_pt 
  if (((p->lastpt) % (2*(p->save_pt))) != 0) {
    cout << "ERROR: save_pt = " << p->save_pt << " entered for grid size " << p->lastpt << endl;
    p->save_pt -= ((p->lastpt) % (2*(p->save_pt)));
    cout << "--> corrected: save_pt = " << p->save_pt << endl;
  }
  // check that p->norm_type is valid
  if (p->norm_type < 0 || p->norm_type > 2) {
    cout << "ERROR: norm_type=" << p->save_pt << " not valid -> using inf-norm" << endl;
    p->norm_type = 0;
  }
  if (p->horizon_search) { cout << "\nsearching for horizon..." << endl; }
  //************************************************************************
  //********************SETTING PARAMS FROM OLD PROGRAM*********************
  if (p->ic2_Amp == 0) {
    if (p->psi_hyp) {
      p->n_ell = 2;
      p->solver = solve_dynamic_psi_hyp;
    }
    else if (p->static_metric) { p->solver = solve_static; }
    else { p->solver = solve_dynamic; }
  }
  else if (p->ic_Amp == 0) {
    if (p->psi_hyp) {
      p->n_ell = 2;
      p->solver = solve2_dynamic_psi_hyp;
    }
    else if (p->static_metric) { p->solver = solve2_static; }
    else { p->solver = solve2_dynamic; }
  }
  else {
    if (p->psi_hyp) {
      p->n_ell = 2;
      p->solver = solveAll_dynamic_psi_hyp;
    }
    else if (p->static_metric) { p->solver = solveAll_static; }
    else { p->solver = solveAll_dynamic; }
  }
  // bbhutil parameters for writing data to sdf
  p->rmin = -(p->rmax);
  p->lastwr = p->lastpt / p->save_pt;
  p->wr_shape = p->lastwr + 1;
  p->zerowr = p->lastwr / 2;
  p->coord_lims[0] = p->rmin; p->coord_lims[1] = p->rmax;
  p->wr_dr = (p->rmax - p->rmin) / ((dbl) p->lastwr);
  p->lastpt = p->lastpt * p->resn_factor;
  p->save_pt = p->save_pt * p->resn_factor;
  p->nsteps = p->nsteps * p->resn_factor;
  p->save_step = p->save_step * p->resn_factor;
  p->outfile = to_string(p->resn_factor) + "-" + p->outfile;
  // derived parameters
  p->npts = p->lastpt + 1;
  p->zeropt = p->lastpt / 2;
  p->norm_factor = 1 / ((dbl) p->npts);
  if (p->norm_type == 1) { p->norm_factor = sqrt(p->norm_factor); }
  p->dr = (p->rmax - p->rmin) / ((dbl) p->lastpt);
  p->dt = p->lam * p->dr;
  p->check_diagnostics = p->save_step * p->check_step;
  p->r[0] = (p->rmin);
  for (int k = 1; k < (p->npts); ++k) {
    p->r[k] = (p->rmin) + k*(p->dr);
    p->r[-k] = (p->r[k]) / (sq(p->r[k]) + (p->lsq));
  }
  // IMPORTANT: to get x/r^2 at r[0] = xmin use -r[-lastpt] for r[-0]
  if (p->r[p->zeropt] != 0 || p->r[p->lastpt] != p->rmax) {
    cout << "\n***COORDINATE INIT ERROR***\n" << endl;
  }
  p->t = 0;
  for (int k = 0; k < p->wr_shape; ++k) {
    (p->inds).push_back({k, (p->save_pt)*k});
  }
  if ((p->inds[p->lastwr]).second != p->lastpt || (p->inds[p->zerowr]).second != p->zeropt) {
    cout << "\n***INDEX INIT ERROR***\n" << endl;
  }
  
  // lapack object declaration
  p->lp_n = (p->n_ell) * (p->npts);
  p->lp_kl = 2;
  p->lp_ku = 2;
  p->lp_nrhs = 1;
  p->lp_ldab = 2*(p->lp_kl) + (p->lp_ku) + 1;
  p->lp_ldb = p->lp_n;
  vector<lapack_int> ipiv_zeros(p->lp_n, 0);
  p->ipiv = ipiv_zeros;
  p->lp_ipiv = &(p->ipiv[0]);

  p->lam2val = 0.5 * (p->lam);
  p->lam6val = (p->lam2val) * (p->one_third);
  p->drsq = (p->dr) * (p->dr);
  p->indr = 1 / (p->dr);
  p->in2dr = 0.5 * (p->indr);
  p->in3dr = (p->one_third) * (p->indr);
  p->indrsq = sq(p->indr);
  p->neg2indrsq = -2 * (p->indrsq);
  p->indt = 1 / (p->dt);
  p->inrmax = 1 / (p->rmax);
  // SPECIFIC TERMS
  p->csomm = 0.75*(p->lam) + 0.5*(p->dt)*(p->inrmax);
  p->csomm_rhs = 1 / (1 + (p->csomm));
  p->csomm_old = 1 - (p->csomm);
  p->jacRR = 3*(p->in2dr) + (p->inrmax);
  p->jacRRm1 = -4 * (p->in2dr);
  p->jacRRm2 = (p->in2dr);
  p->cpsi_rhs = 1 / (p->jacRR);

  // PARAMETER DATA OUTPUT
  str param_data = get_param_string(p_str, p_int, p_dbl, p_bool);
  ofstream specs;
  str specs_name = p->outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << param_data;
  specs.close();
  cout << param_data << endl;
  return 0;
}


#endif

