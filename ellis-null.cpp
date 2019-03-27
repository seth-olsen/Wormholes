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

  if (argc != 2) {
    cout << "\nINCORRECT NUMBER OF ARGUMENTS" << endl;
    return argc;
  }
  else {
    int err_code = params_init(&p, argc, argv);
    if (err_code != 0) {
      cout << "\nPARAM INIT error code = " << err_code << endl;
      return err_code;
    }
  }
  str outfile = p.outfile;

  // derived parameters from coarse file
  int gs = p.lastwr;
  int zeropt = p.zerowr;
  int num_steps = p.nsteps / p.save_step;
  int npts = gs + 1;
  if (zeropt != (gs/2)) {
    cout << "\nERROR WITH GRID SHAPE\n" << endl;
    return -1;
  }

  str c = ",";
  VD outnull(npts, 0.0);
  VD revnull(npts, 0.0);
  VD ps(npts, 0.0);
  vector<dbl *> field_arr {&outnull[0], &revnull[0], &ps[0]};

  str outnull_name = "outnull-" + outfile + ".sdf";
  str revnull_name = "revnull-" + outfile + ".sdf";
  str ps_name = "Ps-" + outfile + ".sdf";
  vector<char *> name_arr {&outnull_name[0], &revnull_name[0], &ps_name[0]};
  
  // output file
  ofstream ofs;
  str outname = "nullCheck-" + outfile + ".csv";
  ofs.open(outname, ofstream::out);
  ofs << p.outfile << "\nstep,time,point,x,r2ps4,outnull,revnull" << endl;
  
  // iterate through time steps
  int times[3];
  gft_set_multi();
  for (int t = 2; t < num_steps; ++t) {
    ofs << t << c << t*(p.save_step)*(p.dt) << endl;
    times[0] = t+1; times[1] = t+1; times[2] = t+1;
    if (read_step(name_arr, times, field_arr, 3) == 0) {
      ofs.close();
      gft_close_all();
      cout << outname << "  written up to " << t << endl;
      return t;
    }
    // iterate through grid points
    for (int j = 0; j < zeropt; ++j) {
      if ((outnull[j] > 0) || (revnull[j] < 0)) {
	ofs << " , ," << j << c << p.r[(p.save_pt)*j] << c << pw4(ps[j])*r2(&p,j);
	if (outnull[j] > 0) { ofs << c << outnull[j]; }
	else { ofs << ", "; }
	if (revnull[j] < 0) { ofs << c << revnull[j] << endl; }
	else { ofs << endl; }
      }
    }
    for (int j = zeropt + 1; j < npts; ++j) {
      if ((outnull[j] < 0) || (revnull[j] > 0)) {
	ofs << " , ," << j << c << p.r[(p.save_pt)*j] << c << pw4(ps[j])*r2(&p,j);
	if (outnull[j] < 0) { ofs << c << outnull[j]; }
	else { ofs << ", "; }
	if (revnull[j] > 0) { ofs << c << revnull[j] << endl; }
	else { ofs << endl; }
      }
    }
  }
  gft_close_all();
  
  cout << outfile << "  written with:" << endl;
  cout << "grid points used = " << npts << "  " << ((p.same_grids) ?
   "(same grids)" : "(dif grids)") << "\ntime steps used = " << num_steps
       << "  " << ((p.same_times) ? "(same times)" : "(dif times)") << endl;
  ofs.close();

  return 0;
}
