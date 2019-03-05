/*

this program will compute

q(t) = ||u_4h(x,t) - u_2h(x,t)||/||u_2h(x,t) - u_h(x,t)||

with (u_4h, u_2h, u_h) from resolution levels (0, 1, 2) 
and save to outfile (.csv) with the line format:
   step_number (= t/dt) , q_phi(t) , q_pi(t) , q_phi+pi(t)

input parameters in terminal as:
(do NOT include .sdf extension in outfile)

 ./p1-ctest <outfile> <lastpt> <save_pt> <nsteps> <save_step>
 <lam> <r2m> <rmin> <rmax> <dspn> <tol> <maxit> <ic_Dsq> 
 <ic_r0> <ic_Amp> <check_step> <zero_pi> <sommerfeld> 
 <dspn_bound> <write_ires> <write_itn> <hold_const>
 <same_times> <same_grids>

where the value of any <parameter> can be set with the
following ordered pair of command line arguments:

 -parameter parameter_value

default values can be found at the start of main()
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include "fda-io.h"
#include "fda-fns.h"

using namespace std;

vector<str> get_unames(str pre, vector<str>& fnames) {
  vector<str> unames;
  for (str name : fnames) { unames.push_back(pre + "-" + name + ".sdf"); }
  return unames;
}

int main(int argc, char **argv)
{
  // coarse simulation parameters
  str outfile = "ellis";
  str outname = "0";
  str pre1 = "0", pre2 = "0", pre3 = "0", pre4 = "0", pre5 = "0",
    pre6 = "0", pre7 = "0", pre8 = "0", pre9 = "0", pre10 = "0",
    pre11 = "0", pre12 = "0", pre13 = "0", pre14 = "0", pre15 = "0",
    pre16 = "0", pre17 = "0", pre18 = "0", pre19 = "0", pre20 = "0", pre21 = "0";
  str ipre1 = "0", ipre2 = "0", ipre3 = "0", ipre4 = "0", ipre5 = "0",
    ipre6 = "0", ipre7 = "0", ipre8 = "0", ipre9 = "0", ipre10 = "0",
    ipre11 = "0", ipre12 = "0", ipre13 = "0", ipre14 = "0", ipre15 = "0",
    ipre16 = "0", ipre17 = "0", ipre18 = "0", ipre19 = "0", ipre20 = "0", ipre21 = "0";
  int lastpt = 1200; // grid size
  int save_pt = 1; // write only every (save_pt)th grid point
  int nsteps = 1200; // time steps
  int save_step = 1; // write only every (save_step)th time step
  dbl lam = 0.25; // dt/dr
  dbl r2m = 0;
  dbl rmax = 60;
  dbl rmin = -rmax;
  dbl dspn = 0.5; // dissipation coefficient
  dbl tol = 0.0000000001; // iterative method tolerance
  dbl ell_tol = 0.01*tol;
  int maxit = 100; // max iterations for debugging
  dbl ic_Dsq = 25.0; // gaussian width
  dbl ic_r0 = 20.0; // gaussian center
  dbl ic_Amp = 0.01; // gaussian amplitude
  bool sommerfeld = true; // sommerfeld condition at outer bound?
  bool dspn_bound = false; // dissipate boundary points?
  bool psi_hyp = true; // psi evolved with hyperbolic eom?
  // variable to hold constant across resolutions
  str hold_const = "lambda"; // "lambda", "dt", or "dr"
  bool same_times = true;
  bool same_grids = true;
  // resolution factors
  int resn0 = 1;
  int resn1 = 2*resn0;
  int resn2 = 4*resn0; // 4h, 2h, and h

  // get parameters from command line
  map<str, str *> p_str {{"-outfile",&outfile}, {"-outname",&outname},
   {"-pre1",&pre1}, {"-pre2",&pre2}, {"-pre3",&pre3}, {"-pre4",&pre4},
   {"-pre5",&pre5}, {"-pre6",&pre6}, {"-pre7",&pre7}, {"-pre8",&pre8},
   {"-pre9",&pre9}, {"-pre10",&pre10}, {"-pre11",&pre11}, {"-pre12",&pre12},
   {"-pre13",&pre13}, {"-pre14",&pre14}, {"-pre15",&pre15}, {"-pre16",&pre16},
   {"-pre17",&pre17}, {"-pre18",&pre18}, {"-pre19",&pre19}, {"-pre20",&pre20},
   {"-pre21",&pre21}, {"-hold_const",&hold_const},
   {"-ipre1",&ipre1}, {"-ipre2",&ipre2}, {"-ipre3",&ipre3}, {"-ipre4",&ipre4},
   {"-ipre5",&ipre5}, {"-ipre6",&ipre6}, {"-ipre7",&ipre7}, {"-ipre8",&ipre8},
   {"-ipre9",&ipre9}, {"-ipre10",&ipre10}, {"-ipre11",&ipre11}, {"-ipre12",&ipre12},
   {"-ipre13",&ipre13}, {"-ipre14",&ipre14}, {"-ipre15",&ipre15}, {"-ipre16",&ipre16},
   {"-ipre17",&ipre17}, {"-ipre18",&ipre18}, {"-ipre19",&ipre19}, {"-ipre20",&ipre20},
   {"-ipre21",&ipre21}};
  map<str, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt", &save_pt},
   {"-nsteps", &nsteps}, {"-save_step",&save_step}, {"-maxit",&maxit},
   {"-resn0", &resn0}, {"-resn1", &resn1}, {"-resn2", &resn2}};
  map<str, dbl *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin},
   {"-rmax",&rmax}, {"-dspn",&dspn}, {"-tol",&tol}, {"-ic_Dsq",&ic_Dsq},
   {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}, {"-ell_tol",&ell_tol}};
  map<str, bool *> p_bool {
   {"-sommerfeld",&sommerfeld}, {"-dspn_bound",&dspn_bound},
   {"-same_times",&same_times}, {"-same_grids",&same_grids},
   {"-psi_hyp",&psi_hyp}};
  map<str, str> params;
  param_collect(argv, argc, params);
  param_set(params, p_str, p_int, p_dbl, p_bool);
  resn1 = 2*resn0;
  resn2 = 4*resn0; // 4h, 2h, and h

  vector<str> prefixes;
  vector<str> iprefixes;
  if (pre1 != "0") { prefixes.push_back(pre1); }
  if (pre2 != "0") { prefixes.push_back(pre2); }
  if (pre3 != "0") { prefixes.push_back(pre3); }
  if (pre4 != "0") { prefixes.push_back(pre4); }
  if (pre5 != "0") { prefixes.push_back(pre5); }
  if (pre6 != "0") { prefixes.push_back(pre6); }
  if (pre7 != "0") { prefixes.push_back(pre7); }
  if (pre8 != "0") { prefixes.push_back(pre8); }
  if (pre9 != "0") { prefixes.push_back(pre9); }
  if (pre10 != "0") { prefixes.push_back(pre10); }
  if (pre11 != "0") { prefixes.push_back(pre11); }
  if (pre12 != "0") { prefixes.push_back(pre12); }
  if (pre13 != "0") { prefixes.push_back(pre13); }
  if (pre14 != "0") { prefixes.push_back(pre14); }
  if (pre15 != "0") { prefixes.push_back(pre15); }
  if (pre16 != "0") { prefixes.push_back(pre16); }
  if (pre17 != "0") { prefixes.push_back(pre17); }
  if (pre18 != "0") { prefixes.push_back(pre18); }
  if (pre19 != "0") { prefixes.push_back(pre19); }
  if (pre20 != "0") { prefixes.push_back(pre20); }
  if (pre21 != "0") { prefixes.push_back(pre21); }
  if (ipre1 != "0") { iprefixes.push_back(ipre1); }
  if (ipre2 != "0") { iprefixes.push_back(ipre2); }
  if (ipre3 != "0") { iprefixes.push_back(ipre3); }
  if (ipre4 != "0") { iprefixes.push_back(ipre4); }
  if (ipre5 != "0") { iprefixes.push_back(ipre5); }
  if (ipre6 != "0") { iprefixes.push_back(ipre6); }
  if (ipre7 != "0") { iprefixes.push_back(ipre7); }
  if (ipre8 != "0") { iprefixes.push_back(ipre8); }
  if (ipre9 != "0") { iprefixes.push_back(ipre9); }
  if (ipre10 != "0") { iprefixes.push_back(ipre10); }
  if (ipre11 != "0") { iprefixes.push_back(ipre11); }
  if (ipre12 != "0") { iprefixes.push_back(ipre12); }
  if (ipre13 != "0") { iprefixes.push_back(ipre13); }
  if (ipre14 != "0") { iprefixes.push_back(ipre14); }
  if (ipre15 != "0") { iprefixes.push_back(ipre15); }
  if (ipre16 != "0") { iprefixes.push_back(ipre16); }
  if (ipre17 != "0") { iprefixes.push_back(ipre17); }
  if (ipre18 != "0") { iprefixes.push_back(ipre18); }
  if (ipre19 != "0") { iprefixes.push_back(ipre19); }
  if (ipre20 != "0") { iprefixes.push_back(ipre20); }
  if (ipre21 != "0") { iprefixes.push_back(ipre21); }
  int nwr = prefixes.size();
  int inwr = iprefixes.size();

  // derived parameters from coarse file
  int gs = lastpt / save_pt;
  int num_steps = nsteps / save_step;
  dbl dr = (rmax - rmin) / ((dbl) lastpt);
  dbl dt = lam * dr;
  int npts0 = gs + 1;
  int npts1 = (same_grids) ? npts0 : 2*gs + 1;
  int npts2 = (same_grids) ? npts0 : 4*gs + 1;

  VD zeros(2, 0.0);
  vector< VD > norms(nwr, VD(2, 0.0));
  vector< VD > u_4h(nwr, VD(npts0, 0.0)),
    u_2h(nwr, VD(npts1, 0.0)), u_h(nwr, VD(npts2, 0.0));
  vector< vector<dbl *> > field_arr;
  
  VD izeros(3, 0.0);
  vector< VD > inorms(inwr, VD(3, 0.0));
  vector< VD > ires_4h(inwr, VD(npts0, 0.0)),
    ires_2h(inwr, VD(npts1, 0.0)), ires_h(inwr, VD(npts2, 0.0));
  vector< vector<dbl *> > ires_arr;

  vector<str> fnames{to_string(resn0) + "-" + outfile,
	to_string(resn1) + "-" + outfile, to_string(resn2) + "-" + outfile};
  vector< vector<str> > unames;
  vector< vector<char *> > name_arr;
  
  vector< vector<str> > inames;
  vector< vector<char *> > iname_arr;
  // output file
  ofstream ofs;
  if (outname == "0") { outname = "conv-" + fnames[0] + ".csv"; }
  ofs.open(outname, ofstream::out);
  ofs <<  "coarse,mid,fine,constant,points,times,same_times,same_grids\n"
      <<  fnames[0]+","+fnames[1]+","+fnames[2]+","+hold_const+"," << npts0
      <<","<< num_steps <<","<< boolalpha << same_times <<","<< same_grids
      << "\n\ndspn,dspn_bound,sommerfeld,tol,maxit\n" << dspn <<","
      << dspn_bound <<","<< sommerfeld <<","<< tol <<","
      << maxit << "\n\nr2m,rmin,rmax,ic_Dsq,ic_r0,ic_Amp\n" << r2m <<","
      << rmin <<","<< rmax <<","<< ic_Dsq <<","<< ic_r0 <<","<< ic_Amp
      << "\n\ncoarse grid:\nlastpt,save_pt,nsteps,save_step,\n" << lastpt
      <<","<< save_pt <<","<< nsteps <<","<< save_step << "\n\nlam,dr,dt,"
      << "tmax\n" << lam <<","<< dr <<","<< dt <<","<< dt*nsteps
      << "\n\ntime";
  for (int k = 0; k < nwr; ++k) {
    ofs << ",Q" << prefixes[k] << "(t)";
    unames.push_back(get_unames(prefixes[k], fnames));
    name_arr.push_back({&unames[k][0][0], &unames[k][1][0], &unames[k][2][0]});
    field_arr.push_back({&u_4h[k][0], &u_2h[k][0], &u_h[k][0]});
  }
  for (int k = 0; k < inwr; ++k) {
    ofs << ",Q" << iprefixes[k] << "4/2(t)" << ",Q" << iprefixes[k] << "2/1(t)";
    inames.push_back(get_unames(iprefixes[k], fnames));
    iname_arr.push_back({&inames[k][0][0], &inames[k][1][0], &inames[k][2][0]});
    ires_arr.push_back({&ires_4h[k][0], &ires_2h[k][0], &ires_h[k][0]});
  }
  ofs << endl;  
  
  // iterate through time steps
  int t1 = ((same_times) ? 1 : 2);
  int t2 = ((same_times) ? 1 : 4);
  int r1 = ((same_grids) ? 1 : 2);
  int r2 = ((same_grids) ? 1 : 4);
  int times[3];

  gft_set_multi();
  for (int t = 0; t < num_steps; ++t) {
    // compute ||u_4h(x,t) - u_2h(x,t)||/||u_2h(x,t) - u_h(x,t)||    
    times[0] = t+1, times[1] = t1*t+1, times[2] = t2*t+1;
    for (int k = 0; k < nwr; ++k) {
      read_step(name_arr[k], times, field_arr[k], 3);
      norms[k] = zeros;
    }
    for (int k = 0; k < inwr; ++k) {
      read_step(iname_arr[k], times, ires_arr[k], 3);
      inorms[k] = izeros;
    }
    // iterate through grid points
    for (int j = 0; j < npts0; ++j) {
      for (int k = 0; k < nwr; ++k) {
	norms[k][0] += pow(u_4h[k][j] - u_2h[k][r1*j], 2);
	norms[k][1] += pow(u_2h[k][r1*j] - u_h[k][r2*j], 2);
      }
      for (int k = 0; k < inwr; ++k) {
	inorms[k][0] += pow(ires_4h[k][j], 2);
	inorms[k][1] += pow(ires_2h[k][r1*j], 2);
	inorms[k][2] += pow(ires_h[k][r2*j], 2);
      }
    }
    // write time, q_u,
    ofs << t*save_step*dt;
    for (int k = 0; k < nwr; ++k) { ofs <<","<< sqrt(norms[k][0]/norms[k][1]); }
    for (int k = 0; k < inwr; ++k) {
      ofs <<","<< sqrt(inorms[k][0]/inorms[k][1]) <<","<< sqrt(inorms[k][1]/inorms[k][2]);
    }
    ofs << endl;
  }
  gft_close_all();
  
  cout << outname << "  written with:" << endl;
  cout << "grid points used = " << npts0 << "  " << ((same_grids) ?
   "(same grids)" : "(dif grids)") << "\ntime steps used = " << num_steps
       << "  " << ((same_times) ? "(same times)" : "(dif times)") << endl;
  ofs.close();

  return 0;
}
