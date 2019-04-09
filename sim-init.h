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
#include "ellis-proc.h"
#include "solvers.h"
#include "sim-header.h"
#include "sim-structs.h"

void bbhp_init(BBHP *bp, PAR *p, str fieldname, VD *p_field, const VD& zeros, D_FN comp)
{
  bp->full_field = p_field;
  if (fieldname == "") { bp->filename = fieldname; }
  else {
    bp->filename = fieldname + "-" + p->outfile + ".sdf";
    bp->file = &(bp->filename[0]);
  } 
  bp->shape = &(p->wr_shape);
  bp->rank = 1;
  bp->coords = &(p->coord_lims[0]);  
  bp->write_this = true;
  if (((bp->full_field) != NULL) && ((p->save_pt) == 1)) {
    bp->data = &((*(bp->full_field))[0]);
  }
  else {
    bp->wr_field = zeros;
    bp->data = &(bp->wr_field[0]);
  }
  bp->compute = comp;
  return;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

vector<BBHP *> writers_init(WRS *wr, FLDS *f, PAR *p)
{
  VD zeros(p->wr_shape, 0);
  vector<BBHP *> out_vec;
  if (p->write_abp) {
    bbhp_init(&(wr->p_Al), p, "Al", &(f->Al), zeros, NULL);
    bbhp_init(&(wr->p_Be), p, "Be", &(f->Be), zeros, NULL);
    bbhp_init(&(wr->p_Ps), p, "Ps", &(f->Ps), zeros, NULL);
    out_vec.push_back(&(wr->p_Al));
    out_vec.push_back(&(wr->p_Be));
    out_vec.push_back(&(wr->p_Ps));
    if (p->write_res) {
      bbhp_init(&(wr->p_resAl), p, "ResAl", NULL, zeros, NULL);
      bbhp_init(&(wr->p_resBe), p, "ResBe", NULL, zeros, NULL);
      bbhp_init(&(wr->p_resPs), p, "ResPs", NULL, zeros, NULL);
    }
  }
  if (p->write_xp) {
    bbhp_init(&(wr->p_Xi), p, "Xi", &(f->Xi), zeros, NULL);
    bbhp_init(&(wr->p_Pi), p, "Pi", &(f->Pi), zeros, NULL);
    out_vec.push_back(&(wr->p_Xi));
    out_vec.push_back(&(wr->p_Pi));
    if (p->write_res) {
      bbhp_init(&(wr->p_resXi), p, "ResXi", &(f->resXi), zeros, NULL);
      bbhp_init(&(wr->p_resPi), p, "ResPi", &(f->resPi), zeros, NULL);
      out_vec.push_back(&(wr->p_resXi));
      out_vec.push_back(&(wr->p_resPi));
    }
  }
  if (p->write_xp2) {
    bbhp_init(&(wr->p_Xi2), p, "Xi2", &(f->Xi2), zeros, NULL);
    bbhp_init(&(wr->p_Pi2), p, "Pi2", &(f->Pi2), zeros, NULL);
    out_vec.push_back(&(wr->p_Xi2));
    out_vec.push_back(&(wr->p_Pi2));
    if (p->write_res) {
      bbhp_init(&(wr->p_resXi2), p, "ResXi2", &(f->resXi2), zeros, NULL);
      bbhp_init(&(wr->p_resPi2), p, "ResPi2", &(f->resPi2), zeros, NULL);
      out_vec.push_back(&(wr->p_resXi2));
      out_vec.push_back(&(wr->p_resPi2));
    }
  }
  if (p->write_ires_abp) {
    bbhp_init(&(wr->p_iresAl), p, "iresAl", NULL, zeros, NULL);
    bbhp_init(&(wr->p_iresBe), p, "iresBe", NULL, zeros, NULL);
    bbhp_init(&(wr->p_iresPs), p, "iresPs", NULL, zeros, NULL);
  }
  if (p->write_ires_xp) {
    bbhp_init(&(wr->p_iresXi), p, "iresXi", NULL, zeros, NULL);
    bbhp_init(&(wr->p_iresPi), p, "iresPi", NULL, zeros, NULL);
  }
  if (p->write_ires_xp2) {
    bbhp_init(&(wr->p_iresXi2), p, "iresXi2", NULL, zeros, NULL);
    bbhp_init(&(wr->p_iresPi2), p, "iresPi2", NULL, zeros, NULL);
  }
  if (p->write_maspect) { bbhp_init(&(wr->p_maspect), p, "maspect", NULL, zeros, NULL); }
  if (p->write_outnull) { bbhp_init(&(wr->p_outnull), p, "outnull", NULL, zeros, NULL); }
  if (p->write_ricci) { bbhp_init(&(wr->p_ricci), p, "ricci", NULL, zeros, NULL); }
  if (p->write_mtot) {
    bbhp_init(&(wr->p_EEtt), p, "EEtt", NULL, zeros, NULL);
    bbhp_init(&(wr->p_EEtx), p, "EEtx", NULL, zeros, NULL);
    bbhp_init(&(wr->p_EExx), p, "EExx", NULL, zeros, NULL);
    bbhp_init(&(wr->p_EEhh), p, "EEhh", NULL, zeros, NULL);
    bbhp_init(&(wr->p_hamiltonian), p, "cHamiltonian", NULL, zeros, NULL);
    bbhp_init(&(wr->p_momentum), p, "cMomentum", NULL, zeros, NULL);
    bbhp_init(&(wr->p_kext), p, "cKext", NULL, zeros, NULL);
    bbhp_init(&(wr->p_dtkext), p, "cDtKext", NULL, zeros, NULL);
  }
  return out_vec;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

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
  dbl amp0 = sqrt((p->lsq)/(p->four_pi));
  int Xeq0 = (p->zeropt);
  for (int k = 0; k < (p->npts); ++k) {
    f->Al[k] = 1;
    //f->Be[k] = 0;
    f->Ps[k] = 1;
    f->Xi[k] = amp0 / r2(p,k);
  }
  
  if (abs(p->ic_Amp) > (p->tol)) {
    if ((p->ic_r0) < (p->ic_Dsq)) {
      for (int k = 0; k < Xeq0; ++k) {
	f->Xi[k] += ic_xi(p->r[k], p->ic_Amp, p->ic_Dsq, p->ic_r0);
	if (!(p->clean_hyp)) { f->Pi[k] = ic_pi(p->r[k], p->ic_Amp, p->ic_Dsq, p->ic_r0); }
      }
    }
    f->Xi[Xeq0] += ic_xi(p->r[Xeq0], p->ic_Amp, p->ic_Dsq, p->ic_r0);
    if ((p->ic_r0) > -(p->ic_Dsq)) {
      for (int k = Xeq0 + 1; k < (p->npts); ++k) {
	f->Xi[k] += ic_xi(p->r[k], p->ic_Amp, p->ic_Dsq, p->ic_r0);
	if (!(p->clean_hyp)) { f->Pi[k] = ic_pi(p->r[k], p->ic_Amp, p->ic_Dsq, p->ic_r0); }
      }
    }
    if (p->sym_pert) {
      for (int k = 0; k < Xeq0; ++k) {
	f->Xi[k] = f->Xi[(p->lastpt) - k];
	if (!(p->clean_hyp)) { f->Pi[k] = f->Pi[(p->lastpt) - k]; }
      }
    }
    if (!(p->clean_hyp)) { f->Pi[Xeq0] = (4*(f->Pi[Xeq0+1]) - (f->Pi[Xeq0+2]) +
					  4*(f->Pi[Xeq0-1]) - (f->Pi[Xeq0-2])) / 6.0; }
  }
  
  if (abs(p->ic2_Amp) > (p->tol)) {
    if ((p->ic2_r0) < (p->ic2_Dsq)) {
      for (int k = 0; k < Xeq0; ++k) {
	f->Xi2[k] += ic_xi(p->r[k], p->ic2_Amp, p->ic2_Dsq, p->ic2_r0);
	if (!(p->clean_hyp)) { f->Pi2[k] = ic_pi(p->r[k], p->ic2_Amp, p->ic2_Dsq, p->ic2_r0); }
      }
    }
    f->Xi2[Xeq0] += ic_xi(p->r[Xeq0], p->ic2_Amp, p->ic2_Dsq, p->ic2_r0);
    if ((p->ic2_r0) > -(p->ic2_Dsq)) {
      for (int k = Xeq0 + 1; k < (p->npts); ++k) {
	f->Xi2[k] += ic_xi(p->r[k], p->ic2_Amp, p->ic2_Dsq, p->ic2_r0);
	if (!(p->clean_hyp)) { f->Pi2[k] = ic_pi(p->r[k], p->ic2_Amp, p->ic2_Dsq, p->ic2_r0); }
      }
    }

    if (p->sym_pert) {
      for (int k = 0; k < Xeq0; ++k) {
	f->Xi2[k] = f->Xi2[(p->lastpt) - k];
	if (!(p->clean_hyp)) { f->Pi2[k] = f->Pi2[(p->lastpt) - k]; }
      }
    }
    if (!(p->clean_hyp)) { f->Pi2[Xeq0] = (4*(f->Pi2[Xeq0+1]) - (f->Pi2[Xeq0+2]) +
					   4*(f->Pi2[Xeq0-1]) - (f->Pi2[Xeq0-2])) / 6.0; }
  }

  f->oldXi = f->Xi;
  f->cnXi = f->Xi;
  f->resXi = zeros;
  f->oldPi = f->Pi;
  f->cnPi = f->Pi;
  f->resPi = zeros;
  if (p->write_ires_xp) {
    f->olderXi = f->Xi;
    f->olderPi = f->Pi;
  }
  f->oldXi2 = f->Xi2;
  f->cnXi2 = f->Xi2;
  f->resXi2 = zeros;
  f->oldPi2 = f->Pi2;
  f->cnPi2 = f->Pi2;
  f->resPi2 = zeros;
  if (p->write_ires_xp2) {
    f->olderXi2 = f->Xi2;
    f->olderPi2 = f->Pi2;
  }

  if (!(p->static_metric)) {
    int t0_itn = solve_t0(f, p);
    if (t0_itn < 0) { return t0_itn; }
  }
  f->oldAl = f->Al;
  f->cnAl = f->Al;
  f->oldBe = f->Be;
  f->cnBe = f->Be;
  
  f->oldPs = f->Ps;
  f->cnPs = f->Ps;
  if (p->psi_hyp) { f->resPs = zeros; }
  if (p->write_ires_abp || p->write_mtot) {
    f->olderAl = f->Al;
    f->olderBe = f->Be;
    f->olderPs = f->Ps;
    f->oldestPs = f->Ps;
  }
  VD ell_res_zeros((p->lp_ldb), 0);
  VD jac_zeros(((p->lp_ldab)*(p->lp_n)), 0);
  f->res_ell = ell_res_zeros;
  f->jac = jac_zeros;
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

int params_init(PAR *p, int argc, char **argv)
{
  map<str, str *> p_str { {"-outfile", &(p->outfile)}, {"-hold_const", &(p->hold_const)} };
  map<str, int *> p_int = get_p_int(p);
  map<str, dbl *> p_dbl = get_p_dbl(p);
  map<str, bool *> p_bool = get_p_bool(p);
  map<str, str> params;
  bool sim_cmd_line = false;
  if (argc > 2) {
    param_collect(argv, argc, params);
    if (((str) argv[0]) == "./ellis") {
      sim_cmd_line = true;
      cout << "\ntaking sim parameters from command line\n" << endl;
    }
  }
  else if (argc == 2) { file_param_collect(argv[1], params); }
  else { file_param_collect("ellis-parameters.txt", params); }
  param_set(params, p_str, p_int, p_dbl, p_bool);
  // check that grid size (lastpt = npts-1) is divisible by 2*save_pt 
  if (((p->lastpt) % (2*(p->save_pt))) != 0) {
    cout << "ERROR: save_pt = " << p->save_pt << " entered for grid size " << p->lastpt << endl;
    p->save_pt -= ((p->lastpt) % (2*(p->save_pt)));
    cout << "--> corrected: save_pt = " << p->save_pt << endl;
  }
  // set SOLVER
  if ((p->ic2_Amp) < (p->tol)) {
    if (p->static_metric) { p->solver = solve_static; }
    else if (p->psi_hyp) {
      p->n_ell = 2;
      p->solver = solve_dynamic_psi_hyp;
    }
    else { p->solver = solve_dynamic; }
  }
  else {
    if (p->static_metric) { p->solver = solveAll_static; }
    else if (p->psi_hyp) {
      p->n_ell = 2;
      p->solver = solveAll_dynamic_psi_hyp;
    }
    else { p->solver = solveAll_dynamic; }
  }
  
  // derived parameters
  //////////////////////////////////////////////////////////////////////////////////////
  // IMPORTANT NOTE: will take (lastpt & nsteps)*=resn if in sim and use cmd line args,
  //                 otherwise lastpt & nsteps will be taken from param file directly
  //////////////////////////////////////////////////////////////////////////////////////
  if (sim_cmd_line) { 
    p->lastpt = (p->lastpt) * (p->resn_factor);
    p->nsteps = (p->nsteps) * (p->resn_factor);
    if (p->same_grids) { p->save_pt = (p->save_pt) * (p->resn_factor); }
    if (p->same_times) { p->save_step = (p->save_step) * (p->resn_factor); }
  }
  p->npts = (p->lastpt) + 1;
  p->zeropt = (p->lastpt) / 2;
  p->norm_factor = 1 / ((dbl) (p->npts));
  if (p->norm_type == 1) { p->norm_factor = sqrt(p->norm_factor); }
  // ***set rmin because can't set negative numbers with user input***
  p->rmin = -(p->rmax);
  p->dr = ((p->rmax) - (p->rmin)) / ((dbl) (p->lastpt));
  p->dt = (p->lam) * (p->dr);
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
  //////////////////////////////////////////////////////////////////////////////////////
  // IMPORTANT NOTE: will take (save_pt & save_step)*=resn if in sim and use cmd line args,
  //                 otherwise save_pt & save_step will be taken from param file directly
  //////////////////////////////////////////////////////////////////////////////////////
  // bbhutil parameters for writing data to sdf
  p->coord_lims[0] = (p->rmin); p->coord_lims[1] = (p->rmax);
  p->lastwr = (p->lastpt) / (p->save_pt);
  p->wr_shape = (p->lastwr) + 1;
  p->zerowr = (p->lastwr) / 2;
  p->wr_dr = ((p->rmax) - (p->rmin)) / ((dbl) (p->lastwr));
  p->check_diagnostics = (p->save_step) * (p->check_step);
  if (p->same_grids) {
    for (int k = 0; k < p->wr_shape; ++k) {
      (p->inds).push_back({k, (p->save_pt)*k});
    }
    if (((p->inds[p->lastwr]).second != (p->lastpt)) || ((p->inds[p->zerowr]).second != (p->zeropt))) {
      cout << "\n***INDEX INIT ERROR***\n" << endl;
    }
  }
  else if ((p->save_pt) != 1) { cout << "\nWARNING: SAME_GRIDS=FALSE and SAVE_PT != 1\n" << endl; }
  else if (((p->lastwr) != (p->lastpt)) || ((p->wr_shape) != (p->npts)) || ((p->zerowr) != (p->zeropt))) {
    cout << "\nERROR: number of grid points not set correctly\n" << endl;
    return -2;
  }
  else { p->wr_dr = (p->dr); }
  
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
  // frequently used
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
  // specific terms
  p->csomm_old = 1 - 0.75*(p->lam) - 0.5*(p->dt)*(p->inrmax);
  p->csomm_rhs = 1 / (2 - (p->csomm_old));
  p->jacRR = 3*(p->in2dr) + (p->inrmax);
  p->jacRRm1 = -4 * (p->in2dr);
  p->jacRRm2 = (p->in2dr);
  p->cpsi_rhs = 1 / (p->jacRR);
  // parameter data output
  if ((sim_cmd_line) || (argc == 1)) {
    ofstream specs;
    str specs_name = to_string(p->resn_factor) + "-" + p->outfile + ".txt";
    specs.open(specs_name, ofstream::out);
    specs << get_paramFile_string(p_str, p_int, p_dbl, p_bool);
    specs.close();
  }
  p->outfile = to_string(p->resn_factor) + "-" + p->outfile;
  cout << get_param_string(p_str, p_int, p_dbl, p_bool) << endl;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

vector<BBHP *> analysis_init(WRS *wr, FLDS *f, PAR *p, int argc, char **argv) {
  bool keep_metric_dt = false;
  bool keep_scalar_dt = false;
  vector<BBHP *> writer_vec;
  int new_argc = 2;
  char *new_argv[] = {argv[0], argv[1]};
  int err_code = params_init(p, new_argc, new_argv);
  if (err_code != 0) {
    cout << "\nPARAM INIT error code = " << err_code << endl;
    return writer_vec;
  }
  if ( ((p->lastwr) != (p->lastpt)) || ((p->wr_shape) != (p->npts))
       || ((p->zerowr) != (p->zeropt)) || ((p->wr_dr) != (p->dr))
       || ((p->save_pt) != 1) || ((p->save_step) != 1) ) {
    cout << "\nERROR: FIELDS NOT WRITTEN AT SIM RESOLUTION\n" << endl;
    return writer_vec;
  }
  VD zeros((p->npts), 0);
  
  f->Xi = zeros;
  f->Pi = zeros;
  f->Xi2 = zeros;
  f->Pi2 = zeros;
  f->Al = zeros;
  f->Be = zeros;
  f->Ps = zeros;

  bbhp_init(&(wr->p_Xi), p, "Xi", &(f->Xi), zeros, NULL);
  bbhp_init(&(wr->p_Pi), p, "Pi", &(f->Pi), zeros, NULL);
  if (p->write_xp2) {
    bbhp_init(&(wr->p_Xi2), p, "Xi2", &(f->Xi2), zeros, NULL);
    bbhp_init(&(wr->p_Pi2), p, "Pi2", &(f->Pi2), zeros, NULL);
  }
  bbhp_init(&(wr->p_Al), p, "Al", &(f->Al), zeros, NULL);
  bbhp_init(&(wr->p_Be), p, "Be", &(f->Be), zeros, NULL);
  bbhp_init(&(wr->p_Ps), p, "Ps", &(f->Ps), zeros, NULL);

  for (int arg_ind = 2; arg_ind < argc; ++arg_ind) {
    str arg = argv[arg_ind];
    if ((arg == "ResAl") || (arg == "ALL") || (arg == "ral")) {
      bbhp_init(&(wr->p_resAl), p, "ResAl", NULL, zeros, compute_ResAl);
      writer_vec.push_back(&(wr->p_resAl));
      keep_metric_dt = true;
    }
    if ((arg == "ResBe") || (arg == "ALL") || (arg == "rbe")) {
      bbhp_init(&(wr->p_resBe), p, "ResBe", NULL, zeros, compute_ResBe);
      writer_vec.push_back(&(wr->p_resBe));
      keep_metric_dt = true;
    }
    if ((arg == "ResPs") || (arg == "ALL") || (arg == "rps")) {
      bbhp_init(&(wr->p_resPs), p, "ResPs", NULL, zeros, compute_ResPs);
      writer_vec.push_back(&(wr->p_resPs));
      keep_metric_dt = true;
    }
    if ((arg == "ResXi") || (arg == "ALL") || (arg == "rxi")) {
      bbhp_init(&(wr->p_resXi), p, "ResXi", NULL, zeros, compute_ResXi);
      writer_vec.push_back(&(wr->p_resXi));
      keep_scalar_dt = true;
    }
    if ((arg == "ResPi") || (arg == "ALL") || (arg == "rpi")) {
      bbhp_init(&(wr->p_resPi), p, "ResPi", NULL, zeros, compute_ResPi);
      writer_vec.push_back(&(wr->p_resPi));
      keep_scalar_dt = true;
    }
    if ((arg == "ResXi2") || (arg == "ALL") || (arg == "rxi2")) {
      bbhp_init(&(wr->p_resXi2), p, "ResXi2", NULL, zeros, compute_ResXi2);
      writer_vec.push_back(&(wr->p_resXi2));
      keep_scalar_dt = true;
    }
    if ((arg == "ResPi2") || (arg == "ALL") || (arg == "rpi2")) {
      bbhp_init(&(wr->p_resPi2), p, "ResPi2", NULL, zeros, compute_ResPi2);
      writer_vec.push_back(&(wr->p_resPi2));
      keep_scalar_dt = true;
    }
    if ((arg == "iresAl") || (arg == "ALL") || (arg == "ial")) {
      bbhp_init(&(wr->p_iresAl), p, "iresAl", NULL, zeros, compute_iresAl);
      writer_vec.push_back(&(wr->p_iresAl));
      keep_metric_dt = true;
    }
    if ((arg == "iresBe") || (arg == "ALL") || (arg == "ibe")) {
      bbhp_init(&(wr->p_iresBe), p, "iresBe", NULL, zeros, compute_iresBe);
      writer_vec.push_back(&(wr->p_iresBe));
      keep_metric_dt = true;
    }
    if ((arg == "iresPs") || (arg == "ALL") || (arg == "ips")) {
      bbhp_init(&(wr->p_iresPs), p, "iresPs", NULL, zeros, compute_iresPs);
      writer_vec.push_back(&(wr->p_iresPs));
      keep_metric_dt = true;
    }
    if ((arg == "iresXi") || (arg == "ALL") || (arg == "ixi")) {
      bbhp_init(&(wr->p_iresXi), p, "iresXi", NULL, zeros, compute_iresXi);
      writer_vec.push_back(&(wr->p_iresXi));
      keep_scalar_dt = true;
    }
    if ((arg == "iresPi") || (arg == "ALL") || (arg == "ipi")) {
      bbhp_init(&(wr->p_iresPi), p, "iresPi", NULL, zeros, compute_iresPi);
      writer_vec.push_back(&(wr->p_iresPi));
      keep_scalar_dt = true;
    }
    if ((arg == "iresXi2") || (arg == "ALL") || (arg == "ixi2")) {
      bbhp_init(&(wr->p_iresXi2), p, "iresXi2", NULL, zeros, compute_iresXi2);
      writer_vec.push_back(&(wr->p_iresXi2));
      keep_scalar_dt = true;
    }
    if ((arg == "iresPi2") || (arg == "ALL") || (arg == "ipi2")) {
      bbhp_init(&(wr->p_iresPi2), p, "iresPi2", NULL, zeros, compute_iresPi2);
      writer_vec.push_back(&(wr->p_iresPi2));
      keep_scalar_dt = true;
    }
    if ((arg == "outnull") || (arg == "ALL") || (arg == "null")) {
      bbhp_init(&(wr->p_outnull), p, "outnull", NULL, zeros, compute_outnull);
      writer_vec.push_back(&(wr->p_outnull));
      keep_metric_dt = true;
    }
    if ((arg == "revnull") || (arg == "ALL") || (arg == "rnull")) {
      bbhp_init(&(wr->p_revnull), p, "revnull", NULL, zeros, compute_revnull);
      writer_vec.push_back(&(wr->p_revnull));
      keep_metric_dt = true;
    }
    if ((arg == "maspect") || (arg == "ALL") || (arg == "mass")) {
      bbhp_init(&(wr->p_maspect), p, "maspect", NULL, zeros, compute_maspect);
      writer_vec.push_back(&(wr->p_maspect));
      keep_metric_dt = true;
    }
    if ((arg == "ricci") || (arg == "ALL") || (arg == "ric")) {
      bbhp_init(&(wr->p_ricci), p, "ricci", NULL, zeros, compute_ricci);
      writer_vec.push_back(&(wr->p_ricci));
      keep_metric_dt = true;
    }
    if ((arg == "EEtt") || (arg == "ALL") || (arg == "tt")) {
      bbhp_init(&(wr->p_EEtt), p, "EEtt", NULL, zeros, compute_EEtt);
      writer_vec.push_back(&(wr->p_EEtt));
      keep_metric_dt = true;
    }
    if ((arg == "EEtx") || (arg == "ALL") || (arg == "tx") || (arg == "xt")) {
      bbhp_init(&(wr->p_EEtx), p, "EEtx", NULL, zeros, compute_EEtx);
      writer_vec.push_back(&(wr->p_EEtx));
      keep_metric_dt = true;
    }
    if ((arg == "EExx") || (arg == "ALL") || (arg == "xx")) {
      bbhp_init(&(wr->p_EExx), p, "EExx", NULL, zeros, compute_EExx);
      writer_vec.push_back(&(wr->p_EExx));
      keep_metric_dt = true;
    }
    if ((arg == "EEhh") || (arg == "ALL") || (arg == "hh")) {
      bbhp_init(&(wr->p_EEhh), p, "EEhh", NULL, zeros, compute_EEhh);
      writer_vec.push_back(&(wr->p_EEhh));
      keep_metric_dt = true;
    }
    if ((arg == "cHamiltonian") || (arg == "ALL") || (arg == "hamiltonian") || (arg == "ham")) {
      bbhp_init(&(wr->p_hamiltonian), p, "cHamiltonian", NULL, zeros, compute_hamiltonian);
      writer_vec.push_back(&(wr->p_hamiltonian));
      keep_metric_dt = true;
    }
    if ((arg == "cMomentum") || (arg == "ALL") || (arg == "momentum") || (arg == "mom")) {
      bbhp_init(&(wr->p_momentum), p, "cMomentum", NULL, zeros, compute_momentum);
      writer_vec.push_back(&(wr->p_momentum));
      keep_metric_dt = true;
    }
    if ((arg == "cKext") || (arg == "ALL") || (arg == "kext")) {
      bbhp_init(&(wr->p_kext), p, "cKext", NULL, zeros, compute_kext);
      writer_vec.push_back(&(wr->p_kext));
      keep_metric_dt = true;
    }
    if ((arg == "cDtKext") || (arg == "ALL") || (arg == "dtkext")) {
      bbhp_init(&(wr->p_dtkext), p, "cDtKext", NULL, zeros, compute_dtkext);
      writer_vec.push_back(&(wr->p_dtkext));
      keep_metric_dt = true;
    }
  }
  if (keep_metric_dt) {
    f->oldAl = f->Al;
    f->oldBe = f->Be;
    f->oldPs = f->Ps;
    f->olderAl = f->Al;
    f->olderBe = f->Be;
    f->olderPs = f->Ps;
    f->oldestPs = f->Ps;
  }
  if (keep_scalar_dt) {
    f->oldXi = f->Xi;  
    f->oldPi = f->Pi;
    f->olderXi = f->Xi;
    f->olderPi = f->Pi;
    f->oldXi2 = f->Xi2;
    f->oldPi2 = f->Pi2;
    f->olderXi2 = f->Xi2;
    f->olderPi2 = f->Pi2;
  }

  if (keep_metric_dt && keep_scalar_dt) {
    wr->get_site_vals = get_all_vals;
    wr->get_bound_vals = get_all_boundvals;
    wr->set_old_fields = set_old_all;
  } // BOTH
  else if (keep_metric_dt && !(keep_scalar_dt)) {
    wr->get_site_vals = get_metric_vals;
    wr->get_bound_vals = get_metric_boundvals;
    wr->set_old_fields = set_old_metric;
  } // JUST METRIC
  else if (!(keep_metric_dt) && keep_scalar_dt) {
    wr->get_site_vals = get_scalar_vals;
    wr->get_bound_vals = get_scalar_boundvals;
    wr->set_old_fields = set_old_scalar;
  } // JUST SCALAR
  else {
    cout << "\nERROR DETERMINING WHICH DIAGNOSTICS TO WRITE" << endl;
    vector<BBHP *> error_msg;
    return error_msg;
  }
  return writer_vec;
}

#endif

