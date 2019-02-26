#ifndef FDA_IO_H_INCLUDED
#define FDA_IO_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <vector> // for everything
#include <cmath> // for ICs
#include <string> // for parameter input
#include <sstream> // for printing doubles to full precision
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "sim-structs.h"
#include "sim-header.h"
#include "fda-fns.h"
#include "ellis-fns.h"
#include "ellis-proc.h"

inline str bool_to_str(bool bool_in)
{
  return ( (bool_in) ? "TRUE" : "FALSE" );
}

inline void write_sdf(const BBHP *bp, dbl t)
{
  gft_out_bbox(bp->file, t, bp->shape, bp->rank, bp->coords, bp->data);
}

inline void prepare_write(const VD& fld, VD& wr_fld, PAR *p)
{
  for (auto ind : (p->inds)) {
    wr_fld[ind.first] = fld[ind.second];
  }
}

inline void write_bbhp(BBHP *bp, PAR *p)
{
  prepare_write(*(bp->full_field), bp->wr_field, p);
  gft_out_bbox(bp->file, p->t, bp->shape, bp->rank, bp->coords, bp->data);
}

inline void write_bbhp_vec(vector<BBHP *>& bp_vec, PAR *p)
{
  for (BBHP *bp : bp_vec) { write_bbhp(bp, p); }
}

VD make_vector(int len, dbl val)
{
  VD new_vec(len, val);
  return new_vec;
}

map<str, int *> get_p_int(PAR *p)
{
  map<str, int *> p_int {{"-lastpt",&(p->lastpt)}, {"-save_pt",&(p->save_pt)},
   {"-nsteps",&(p->nsteps)}, {"-save_step",&(p->save_step)},
   {"-norm_type",&(p->norm_type)}, {"-maxit",&(p->maxit)},
   {"-check_step",&(p->check_step)}, {"-resn_factor",&(p->resn_factor)}};
  return p_int;
}
map<str, dbl *> get_p_dbl(PAR *p)
{
  map<str, dbl *> p_dbl {{"-lam",&(p->lam)}, {"-ell_up_weight",&(p->ell_up_weight)},
   {"-rmin",&(p->rmin)}, {"-rmax",&(p->rmax)}, {"-dspn",&(p->dspn)},
   {"-tol",&(p->tol)}, {"-ell_tol",&(p->ell_tol)}, {"-lsq",&(p->lsq)},
   {"-ic_Dsq",&(p->ic_Dsq)}, {"-ic_r0",&(p->ic_r0)}, {"-ic_Amp",&(p->ic_Amp)},
   {"-ic2_Dsq",&(p->ic2_Dsq)}, {"-ic2_r0",&(p->ic2_r0)}, {"-ic2_Amp",&(p->ic2_Amp)}};
  return p_dbl;
}

map<str, bool *> get_p_bool(PAR *p)
{
  map<str, bool *> p_bool { {"-psi_hyp",&(p->psi_hyp)}, {"-dspn_psi",&(p->dspn_psi)},
   {"-dspn_bound",&(p->dspn_bound)}, {"-static_metric",&(p->static_metric)},
   {"-somm_cond",&(p->somm_cond)}, {"-horizon_search",&(p->horizon_search)},
   {"-write_res",&(p->write_res)}, {"-write_ricci",&(p->write_ricci)},
   {"-write_itn",&(p->write_itn)}, {"-write_mtot",&(p->write_mtot)},
   {"-write_maspect",&(p->write_maspect)}, {"-write_outnull",&(p->write_outnull)},
   {"-write_xp",&(p->write_xp)}, {"-write_ires_xp",&(p->write_ires_xp)},
   {"-write_abp",&(p->write_abp)}, {"-write_ires_abp",&(p->write_ires_abp)},
   {"-write_xp2",&(p->write_xp2)}, {"-write_ires_xp2",&(p->write_ires_xp2)},
   {"-clean_hyp",&(p->clean_hyp)}, {"-clean_ell",&(p->clean_ell)}};
  return p_bool;
}

void param_collect(char **source, int num, map<str, str>& dest) {
  for (int arg = 1; arg < num; ++arg) {
    if (source[arg][0] == '-') {
      dest[source[arg]] = source[arg+1];
    }
  }
}

void file_param_collect(str filename, map<str, str>& dest) {
  vector<str> source;
  str entry;
  ifstream p_file(filename);
  if (p_file.is_open()) {
    while (getline(p_file, entry)) { source.push_back(entry); }
    for (unsigned int arg = 0; arg < source.size(); ++arg) {
      if (source[arg][0] == '-') { dest[source[arg]] = source[arg+1]; }
    }
  }
  else { cout << "\n***UNABLE TO OPEN PARAMETER FILE***\n" << endl; }
}

void param_set(map<str, str>& p_all, map<str, str *>& p_str,
	       map<str, int *>& p_int, map<str, dbl *>& p_dbl,
	       map<str, bool *>& p_bool) {
  for (pair<str, str> p : p_all) {
    if (p_str.count(p.first)) { *p_str[p.first] = p.second; }
    else if (p_int.count(p.first)) { *p_int[p.first] = atoi(&p.second[0]); }
    else if (p_dbl.count(p.first)) { *p_dbl[p.first] = atof(&p.second[0]); }
    else if (p_bool.count(p.first)) { *p_bool[p.first] = (bool) atoi(&p.second[0]); }
  }
}

// read fields using bbhutil
void read_step(const vector<char *>& files, int times[], const vector<dbl *>& fields, int nfields) {
  for (int k = 0; k < nfields; ++k) {
    gft_read_brief(files[k], times[k], fields[k]);
  }
  return;
}

void write_diagnostics(WRS *wr, FLDS *f, PAR *p)
{
  if ((p->write_res) && (p->write_abp)) {
    int npts = (p->npts);
    VD res_0(3*npts, 0);
    dbl res = get_res_abp_t0(res_0, f, p);
    if (res > p->ell_tol) { cout << p->t << " res = " << res << endl; }
    for (int k = 0; k < p->wr_shape; ++k) {
      (wr->p_resAl).wr_field[k] = res_0[(p->inds[k]).second];
      (wr->p_resBe).wr_field[k] = res_0[npts + (p->inds[k]).second];
      (wr->p_resPs).wr_field[k] = res_0[2*npts + (p->inds[k]).second];
    }
    write_sdf(&(wr->p_resAl), p->t);
    write_sdf(&(wr->p_resBe), p->t);
    write_sdf(&(wr->p_resPs), p->t);
  }
  // write ires
  if (p->write_ires_xp) {
    get_ires_xp(wr, f, p);
    write_sdf(&(wr->p_iresXi), p->t);
    write_sdf(&(wr->p_iresPi), p->t);
  }
  if (p->write_ires_xp2) {
    get_ires_xp2(wr, f, p);
    write_sdf(&(wr->p_iresXi2), p->t);
    write_sdf(&(wr->p_iresPi2), p->t);
  }
  if (p->write_ires_abp) {
    get_ires_abp(wr, f, p);
    write_sdf(&(wr->p_iresAl), p->t);
    write_sdf(&(wr->p_iresBe), p->t);
    write_sdf(&(wr->p_iresPs), p->t);
  }
  // write ricci
  if (p->write_ricci) {
    get_ricci((wr->p_ricci).wr_field, f->Xi, f->Pi, f->Xi2, f->Pi2, f->Ps, p->inds);
    write_sdf(&(wr->p_ricci), p->t);
  }
  // write outnull
  if (p->write_outnull) {
    get_outnull((wr->p_outnull).wr_field, f->Al, f->Be, f->Ps, p);
    write_sdf(&(wr->p_outnull), p->t);
  }
  // write maspect
  if (p->write_maspect) {
    get_maspect((wr->p_maspect).wr_field, f->Al, f->Be, f->Ps, p);
    write_sdf(&(wr->p_maspect), p->t);
  }
  // ******************************************************************
  // EINSTEIN EQUATION RESIDUALS
  if (p->write_mtot) {
    get_EE_adm(wr, f, p);
    write_sdf(&(wr->p_EEtt), p->t); // einstein tt residual
    write_sdf(&(wr->p_EEtx), p->t); // einstein tr residual
    write_sdf(&(wr->p_EExx), p->t); // einstein rr residual
    write_sdf(&(wr->p_EEhh), p->t);  // einstein thth residual
    write_sdf(&(wr->p_hamiltonian), p->t); // hamiltonian constraint residual
    write_sdf(&(wr->p_momentum), p->t); // momentum constraint residual
    write_sdf(&(wr->p_kext), p->t); // maximal slicing residual
    write_sdf(&(wr->p_dtkext), p->t);  // slicing time derivative residual
  }
  // ******************************************************************
  return;
}

/*
// write fields using bbhutil
void wr_step(int nfields, const vector<char *>& files,
	     dbl time, int *shape, int rank, dbl *coordinates,
	     const vector<dbl *>& fields) {
  for (int k = 0; k < nfields; ++k) {
    gft_out_bbox(files[k], time, shape, rank, coordinates, fields[k]);
  }
  return;
}

// GET COARSENED ARRAY FOR WRITING
void get_wr_f(const VD& f, VD& wr, int wr_shape, int save_pt)
{
  int k, s = 0;
  for (k = 0; k < wr_shape; ++k) {
    wr[k] = f[s];
    s += save_pt;
  }
  return;
}
*/

str get_param_string(map<str, str *>& p_str, map<str, int *>& p_int,
		     map<str, dbl *>& p_dbl, map<str, bool *>& p_bool)
{
  str param_string = "\nPARAMETERS:\n\n";
  for (pair<str, str *> param : p_str) {
    str spaces = (((param.first).size() < 8) ? "\t\t\t=   " : "\t\t=   ");
    param_string += param.first + spaces + *(param.second) + "\n";
  }
  param_string += "\n";
  for (pair<str, int *> param : p_int) {
    str spaces = (((param.first).size() < 8) ? "\t\t\t=   " : "\t\t=   ");
    param_string += param.first + spaces + to_string(*(param.second)) + "\n";
  }
  param_string += "\n";
  for (pair<str, dbl *> param : p_dbl) {
    str spaces = (((param.first).size() < 8) ? "\t\t\t=   " : "\t\t=   ");
    param_string += param.first + spaces + to_string(*(param.second)) + "\n";
  }
  param_string += "\n";
  for (pair<str, bool *> param : p_bool) {
    str spaces = (((param.first).size() < 8) ? "\t\t\t=   " : "\t\t=   ");
    param_string += param.first + spaces + bool_to_str(*(param.second)) + "\n";
  }
  return param_string;
}

void record_horizon(PAR *p, const VD& f_ps, int ind, int itn, int t_itn)
{
  cout << "Horizon Found at:\nr[" << ind << "] = " << p->r[ind] << "  (r_areal = "
       << sq(f_ps[ind])*p->r[ind] << ")\nt[" << t_itn << "] = "
       << p->t << "  (itn " << itn << ")\n" << endl;
  str param_str = "\noutfile name = " + p->outfile + "\nlaspt (save_pt) = " +
    to_string(p->lastpt) + " (" + to_string(p->save_pt) + ")\nnsteps (save_step) = " +
    to_string(p->nsteps) + " (" + to_string(p->save_step) + ")\nrmin = " + to_string(p->rmin)
    + ";\trmax = " + to_string(p->rmax) + "\nlambda = " + to_string(p->lam) +
    "\ndt = " + to_string(p->dt) + "\ndissipation = " + to_string(p->dspn) +
    "\nell_up_weight = " + to_string(p->ell_up_weight) + "\nhyp_tol = " + to_string(p->tol) +
    "\nell_tol = " + to_string(p->ell_tol) + "\nmaxit = " + to_string(p->maxit) + "\nic_r0 = " +
    to_string(p->ic_r0) + "\nic_Amp = " + to_string(p->ic_Amp) + "\nic_Dsq = " +
    to_string(p->ic_Dsq) + "\ndr = " + to_string(p->dr);
  param_str += "\noptions:\nhyperbolic psi evolution = " + bool_to_str(p->psi_hyp) +
    "\nsommerfeld bc = " + bool_to_str(p->somm_cond)
    + "\ndissipation at bound = " + bool_to_str(p->dspn_bound) + "\nclean hyperbolic update functions = "
    + bool_to_str(p->clean_hyp) + "\nclean elliptic update functions = " + bool_to_str(p->clean_ell) + "\n";
  ofstream specs;
  str specs_name = "horizon-" + p->outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << "Horizon Found at:\nr[" << ind << "] = " << p->r[ind] << "  (r_areal = "
	<< sq(f_ps[ind])*p->r[ind] << ")\nt[" << t_itn << "] = "
	<< p->t << "  (itn " << itn << ")\nUsing:\ndr = " << p->dr << "\nic_Amp = " << p->ic_Amp
	<< "\nic_Dsq = " << p->ic_Dsq << "\nic_r0 = " << p->ic_r0 << "\n\n\nFULL PARAMETER DATA:\n";
  specs << param_str;
  specs.close();
  cout << param_str;
  return;
}

#endif


