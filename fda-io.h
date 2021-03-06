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
#include <cctype> // for isdigit()
#include "bbhutil.h" // for output to .sdf
#include "sim-structs.h"
#include "sim-header.h"
#include "fda-fns.h"
#include "ellis-fns.h"
#include "ellis-proc.h"

map<int, str> get_field_names()
{
  map<int, str> names { {0, "Al"}, {1, "Be"}, {2, "Ps"},
			{3, "Xi"}, {4, "Pi"}, {5, "Xi2"}, {6, "Pi2"}};
  return names;
}

map<int, str> get_diagnostic_names()
{
  map<int, str> names { {0, "ResAl"}, {1, "ResBe"}, {2, "ResPs"},
			{3, "ResXi"}, {4, "ResPi"}, {5, "ResXi2"}, {6, "ResPi2"},
			{7, "iresAl"}, {8, "iresBe"}, {9, "iresPs"},
			{10, "iresXi"}, {11, "iresPi"}, {12, "iresXi2"}, {13, "iresPi2"},
			{14, "outnull"}, {15, "revnull"}, {16, "maspect"}, {17, "ricci"},
			{18, "cHamiltonian"}, {19, "cMomentum"},
			{20, "cKext"}, {21, "cDtKext"},
			{22, "EEtt"}, {23, "EEtx"}, {24, "EExx"}, {25, "EEhh"} };
  return names;
}

inline str bool_to_str(bool bool_in)
{
  return ( (bool_in) ? "TRUE" : "FALSE" );
}

inline void write_sdf(const BBHP *bp, dbl t)
{
  gft_out_bbox(bp->file, t, bp->shape, bp->rank, bp->coords, bp->data);
}

inline void write_sdf_direct(char *file, dbl *data, PAR *p)
{
  gft_out_bbox(file, (p->t), &(p->wr_shape), 1, &(p->coord_lims[0]), data);
}

inline void prepare_write(const VD& fld, VD& wr_fld, PAR *p)
{
  for (auto ind : (p->inds)) {
    wr_fld[ind.first] = fld[ind.second];
  }
}

inline void write_bbhp(BBHP *bp, PAR *p)
{
  if (((p->save_pt) == 1) || ((bp->full_field) == NULL)) {
    gft_out_bbox(bp->file, p->t, bp->shape, bp->rank, bp->coords, bp->data);
  }
  else {
    prepare_write(*(bp->full_field), bp->wr_field, p);
    gft_out_bbox(bp->file, p->t, bp->shape, bp->rank, bp->coords, bp->data);
  }
}

inline void write_bbhp_vec(vector<BBHP *>& bp_vec, PAR *p)
{
  for (BBHP *bp : bp_vec) { write_bbhp(bp, p); }
}

inline int read_bbhp(BBHP *bp, int time)
{
  return gft_read_brief(bp->file, time, bp->data);
}
inline int read_bbhp_vec(const vector<BBHP *>& bpv, int time)
{
  for (BBHP *bp : bpv) {
    if (gft_read_brief(bp->file, time, bp->data) == 0) { return 0; }
  }
  return 1;
}

inline void compute_bbhp_vec(const vector<BBHP *>& bpv, SSV& s, int k)
{
  for (BBHP *bp : bpv) { (bp->wr_field)[k] = (bp->compute)(s); }
}

// read fields using bbhutil
int read_step(const vector<char *>& files, int times[], const vector<dbl *>& fields, int nfields)
{
  for (int k = 0; k < nfields; ++k) {
    if (gft_read_brief(files[k], times[k], fields[k]) == 0) { return 0; }
  }
  return 1;
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
   {"-check_step",&(p->check_step)}, {"-resn_factor",&(p->resn_factor)},
   {"-signal_code",&(p->signal_code)}};
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
   {"-clean_hyp",&(p->clean_hyp)}, {"-clean_ell",&(p->clean_ell)},
   {"-same_times",&(p->same_times)}, {"-same_grids",&(p->same_grids)},
   {"-sym_pert",&(p->sym_pert)}};
  return p_bool;
}

void param_collect(char **source, int num, map<str, str>& dest) {
  for (int arg = 1; arg < num; ++arg) {
    if (source[arg][0] == '-') {
      if (isdigit(source[arg][1]) == 0) { dest[source[arg]] = source[arg+1]; }
    }
  }
}

void file_param_collect(str filename, map<str, str>& dest) {
  vector<str> source;
  str entry;
  ifstream p_file(filename);
  if (p_file.is_open()) {
    while (getline(p_file, entry)) {
      if (entry == "TRUE") { source.push_back("1"); }
      else if (entry == "FALSE") { source.push_back("0"); }
      else { source.push_back(entry); }
    }
    for (unsigned int arg = 0; arg < source.size(); ++arg) {
      if (source[arg][0] == '-') {
	if (isdigit(source[arg][1]) == 0) { dest[source[arg]] = source[arg+1]; }
      }
    }
  }
  else { cout << "\n***UNABLE TO OPEN PARAMETER FILE***\n" << endl; }
}

void param_set(map<str, str>& p_all, map<str, str *>& p_str,
	       map<str, int *>& p_int, map<str, dbl *>& p_dbl,
	       map<str, bool *>& p_bool)
{
  for (pair<str, str> p : p_all) {
    if (p_str.count(p.first)) { *p_str[p.first] = p.second; }
    else if (p_int.count(p.first)) { *p_int[p.first] = atoi(&p.second[0]); }
    else if (p_dbl.count(p.first)) { *p_dbl[p.first] = atof(&p.second[0]); }
    else if (p_bool.count(p.first)) { *p_bool[p.first] = (bool) atoi(&p.second[0]); }
  }
}

inline void write_notice(str outname, str message)
{
  cout << "\n\n" << message << endl;
  ofstream notice;
  notice.open(outname, ofstream::out);
  notice << message << endl;
  notice.close();
  return;
}

void write_tseries(WRS *wr, FLDS *f, PAR *p)
{
  dbl coords[2] = {-(p->dt), p->t};
  if (p->write_outnull) {
    int shape = (wr->areal_min).size();
    str nm_armin = "min-" + (wr->p_areal).filename;
    str nm_xarmin = "x" + nm_armin;
    str nm_null0 = "zero-" + (wr->p_outnull).filename;
    str nm_rev0 = "zero-" + (wr->p_revnull).filename;
    gft_out_bbox(&nm_armin[0], 0, &shape, 1, &(coords[0]), &(wr->areal_min[0]));
    gft_out_bbox(&nm_xarmin[0], 0, &shape, 1, &(coords[0]), &(wr->xareal_min[0]));
    gft_out_bbox(&nm_null0[0], 0, &shape, 1, &(coords[0]), &(wr->outnull_0[0]));
    gft_out_bbox(&nm_rev0[0], 0, &shape, 1, &(coords[0]), &(wr->revnull_0[0]));
  }
  if (p->write_ricci) {
    int shape = (wr->ricci_max).size();
    str nm_ricmax = "max-" + (wr->p_ricci).filename;
    str nm_xricmax = "x" + nm_ricmax;
    gft_out_bbox(&nm_ricmax[0], 0, &shape, 1, &(coords[0]), &(wr->ricci_max[0]));
    gft_out_bbox(&nm_xricmax[0], 0, &shape, 1, &(coords[0]), &(wr->xricci_max[0]));
  }
  return;
}

void write_diagnostics(WRS *wr, FLDS *f, PAR *p)
{
  if ((p->write_res) && (p->write_abp)) {
    int npts = (p->npts);
    VD res_0(3*npts, 0);
    dbl res = get_res_abp_t0(res_0, f, p);
    if (res > p->dr) { cout << p->t << " res = " << res << endl; }
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
    get_ricci(f, p);
    int ind_max = distance((f->ricci).begin(), max_element((f->ricci).begin(), (f->ricci).end()));
    (wr->ricci_max).push_back(f->ricci[ind_max]);
    (wr->xricci_max).push_back(p->r[ind_max]);
    write_bbhp(&(wr->p_ricci), p);
  }
  // write outnull
  if (p->write_outnull) {
    get_nullex(f, p);
    int indL = (p->zeropt)/2;
    int indR = 3*indL;
    for (int k = indL; k < indR; ++k) {
      if ((f->revnull[k]) <= 0) {
	(wr->revnull_0).push_back(p->r[k]);
	break;
      }
    }
    for (int k = indR; k > indL; --k) {
      if ((f->outnull[k]) <= 0) {
	(wr->outnull_0).push_back(p->r[k]);
	break;
      }
    }
    get_areal(f, p);
    int ind_min = distance((f->areal).begin(), min_element((f->areal).begin(), (f->areal).end()));
    (wr->areal_min).push_back(f->areal[ind_min]);
    (wr->xareal_min).push_back(p->r[ind_min]);
    write_bbhp(&(wr->p_outnull), p);
    write_bbhp(&(wr->p_revnull), p);
    write_bbhp(&(wr->p_areal), p);
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

str get_paramFile_string(map<str, str *>& p_str, map<str, int *>& p_int,
			 map<str, dbl *>& p_dbl, map<str, bool *>& p_bool)
{
  str param_string = "-outfile\n" + *(p_str["-outfile"]) + "\n";
  for (pair<str, int *> param : p_int) {
    param_string += param.first + "\n" + to_string(*(param.second)) + "\n";
  }
  for (pair<str, dbl *> param : p_dbl) {
    param_string += param.first + "\n" + to_string(*(param.second)) + "\n";
  }
  for (pair<str, bool *> param : p_bool) {
    param_string += param.first + "\n" + bool_to_str(*(param.second)) + "\n";
  }
  param_string += "-hold_const\n" + *(p_str["-hold_const"]) + "\n";
  return param_string;
}

void record_horizon(PAR *p, const VD& f_ps, int ind, int itn, int t_itn)
{
  str horizon_msg = "Horizon Found at:\nr[" + to_string(ind) + "] = " + to_string(p->r[ind])
    + "  (r_areal = " + to_string(sq(f_ps[ind])*sqrt(r2(p,ind))) + ")\nt[" + to_string(t_itn)
    + "] = " + to_string(p->t) + "  (itn " + to_string(itn) + ")\nUsing:\ndr = " + to_string(p->dr)
    + "\nic_Amp = " + to_string(p->ic_Amp) + "\nic_Dsq = " + to_string(p->ic_Dsq) + "\nic_r0 = "
    + to_string(p->ic_r0) + "\n\nic2_Amp = " + to_string(p->ic2_Amp) + "\nic2_Dsq = "
    + to_string(p->ic2_Dsq) + "\nic2_r0 = " + to_string(p->ic2_r0) + "\n\nfile = " + p->outfile
    + "\nresolution = " + to_string(p->resn_factor) + "\npsi_hyp = " + bool_to_str(p->psi_hyp)
    + "\nsym_pert = " + bool_to_str(p->sym_pert);
  cout << horizon_msg << "\n\nAPPARENT HORIZON EXIT\n" << endl;
  ofstream specs;
  str specs_name = "horizon-" + p->outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << horizon_msg << endl;
  specs.close();
  return;
}

void record_antihorizon(PAR *p, const VD& f_ps, int ind, int itn, int t_itn)
{
  str horizon_msg = "Anti-Horizon Found at:\nr[" + to_string(ind) + "] = " + to_string(p->r[ind])
    + "  (r_areal = " + to_string(sq(f_ps[ind])*sqrt(r2(p,ind))) + ")\nt[" + to_string(t_itn)
    + "] = " + to_string(p->t) + "  (itn " + to_string(itn) + ")\nUsing:\ndr = " + to_string(p->dr)
    + "\nic_Amp = " + to_string(p->ic_Amp) + "\nic_Dsq = " + to_string(p->ic_Dsq) + "\nic_r0 = "
    + to_string(p->ic_r0) + "\n\nic2_Amp = " + to_string(p->ic2_Amp) + "\nic2_Dsq = "
    + to_string(p->ic2_Dsq) + "\nic2_r0 = " + to_string(p->ic2_r0) + "\n\nfile = " + p->outfile
    + "\nresolution = " + to_string(p->resn_factor) + "\npsi_hyp = " + bool_to_str(p->psi_hyp)
    + "\nsym_pert = " + bool_to_str(p->sym_pert);
  cout << horizon_msg << "\n\nANTI-HORIZON EXIT\n" << endl;
  ofstream specs;
  str specs_name = "antihorizon-" + p->outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << horizon_msg << endl;
  specs.close();
  return;
}

void record_exit(PAR *p, int code, int itn, int t_itn)
{
  str exit_msg = "EXIT without horizon at t[" + to_string(t_itn) + "] = " + to_string(p->t)
    + "\ncode = " + to_string(code) + "\nitn = " + to_string(itn) + ")\nUsing:\ndr = " + to_string(p->dr)
    + "\nic_Amp = " + to_string(p->ic_Amp) + "\nic_Dsq = " + to_string(p->ic_Dsq) + "\nic_r0 = "
    + to_string(p->ic_r0) + "\n\nic2_Amp = " + to_string(p->ic2_Amp) + "\nic2_Dsq = "
    + to_string(p->ic2_Dsq) + "\nic2_r0 = " + to_string(p->ic2_r0) + "]n\nfile = " + p->outfile
    + "\nresolution = " + to_string(p->resn_factor) + "\npsi_hyp = " + bool_to_str(p->psi_hyp)
    + "\nsym_pert = " + bool_to_str(p->sym_pert);
  cout << exit_msg << "\n\nEXIT WITHOUT HORIZON\n" << endl;
  ofstream specs;
  str specs_name = "EXIT-" + p->outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << exit_msg << endl;
  specs.close();
  return;
}

// for ellis-conv
vector<str> get_unames(str pre, vector<str>& fnames) {
  vector<str> unames;
  for (str name : fnames) { unames.push_back(pre + "-" + name + ".sdf"); }
  return unames;
}

#endif


