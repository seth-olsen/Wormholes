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
#include "sim-init.h"

int main(int argc, char **argv)
{
  PAR p;
  int err_code = params_init(&p, argc, argv);
  if (err_code != 0) {
    cout << "\nPARAM INIT error code = " << err_code << endl;
    return err_code;
  }
  int resn0 = p.resn_factor;
  int resn1 = 2*resn0;
  int resn2 = 4*resn0; // 4h, 2h, and h

  vector<str> prefixes;
  vector<str> iprefixes;
  if (p.write_xp) {
    prefixes.push_back("Xi");
    prefixes.push_back("Pi");
  }
  if (p.write_xp2) {
    prefixes.push_back("Xi2");
    prefixes.push_back("Pi2");
  }
  if (p.write_abp) {
    prefixes.push_back("Al");
    prefixes.push_back("Be");
    prefixes.push_back("Ps");
    if (p.write_res) { iprefixes.push_back("ResPs"); }
  }
  if (p.write_ires_xp) {
      iprefixes.push_back("iresXi");
      iprefixes.push_back("iresPi");
  }
  if (p.write_ires_xp2) {
      iprefixes.push_back("iresXi2");
      iprefixes.push_back("iresPi2");
  }
  if (p.write_ires_abp) {
      iprefixes.push_back("iresAl");
      iprefixes.push_back("iresBe");
      iprefixes.push_back("iresPs");
  }
  if (p.write_outnull) { prefixes.push_back("outnull"); }
  if (p.write_maspect) { prefixes.push_back("maspect"); }
  if (p.write_ricci) { prefixes.push_back("ricci"); }
  if (p.write_mtot) {
      iprefixes.push_back("EEtt");
      iprefixes.push_back("EEtx");
      iprefixes.push_back("EExx");
      iprefixes.push_back("EEhh");
      iprefixes.push_back("cHamiltonian");
      iprefixes.push_back("cMomentum");
      iprefixes.push_back("cKext");
      iprefixes.push_back("cDtKext");
  }
  int nwr = prefixes.size();
  int inwr = iprefixes.size();
  str outfile = p.outfile;
  outfile.erase(0, 2);

  // derived parameters from coarse file
  int gs = p.lastwr;
  int num_steps = p.nsteps / p.save_step;
  int npts0 = gs + 1;
  int npts1 = (p.same_grids) ? npts0 : 2*gs + 1;
  int npts2 = (p.same_grids) ? npts0 : 4*gs + 1;

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
  if (fnames[0] != p.outfile) { cout << "\n**ERROR IN FILENAME**\n" << endl; }
  vector< vector<str> > unames;
  vector< vector<char *> > name_arr;
  
  vector< vector<str> > inames;
  vector< vector<char *> > iname_arr;
  // output file
  ofstream ofs;
  str outname = "conv-" + fnames[0] + ".csv";
  ofs.open(outname, ofstream::out);
  ofs << p.outfile << "\ntime";
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
  int t1 = ((p.same_times) ? 1 : 2);
  int t2 = ((p.same_times) ? 1 : 4);
  int r1 = ((p.same_grids) ? 1 : 2);
  int r2 = ((p.same_grids) ? 1 : 4);
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
    ofs << t*(p.save_step)*(p.dt);
    for (int k = 0; k < nwr; ++k) { ofs <<","<< sqrt(norms[k][0]/norms[k][1]); }
    for (int k = 0; k < inwr; ++k) {
      ofs <<","<< sqrt(inorms[k][0]/inorms[k][1]) <<","<< sqrt(inorms[k][1]/inorms[k][2]);
    }
    ofs << endl;
  }
  gft_close_all();
  
  cout << outname << "  written with:" << endl;
  cout << "grid points used = " << npts0 << "  " << ((p.same_grids) ?
   "(same grids)" : "(dif grids)") << "\ntime steps used = " << num_steps
       << "  " << ((p.same_times) ? "(same times)" : "(dif times)") << endl;
  ofs.close();

  return 0;
}
