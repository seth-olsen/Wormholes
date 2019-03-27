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
#include <iterator> // for distance()

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
  VD maspect(npts, 0.0);
  VD ps(npts, 0.0);
  vector<dbl *> field_arr {&maspect[0], &ps[0]};

  str maspect_name = "maspect-" + outfile + ".sdf";
  str ps_name = "Ps-" + outfile + ".sdf";
  vector<char *> name_arr {&maspect_name[0], &ps_name[0]};
  
  // output file
  ofstream ofs;
  str outname = "massCheck-" + outfile + ".csv";
  ofs.open(outname, ofstream::out);
  ofs << p.outfile << "\nstep,time,m(-R),m(0),m(R),m_max,at,x,r2ps4,m_min,at,x,r2ps4" << endl;
  
  // iterate through time steps
  int times[2];
  gft_set_multi();
  for (int t = 2; t < num_steps; ++t) {
    times[0] = t+1; times[1] = t+1;
    if (read_step(name_arr, times, field_arr, 2) == 0) { return t; }
    ofs << t << c << t*(p.save_step)*(p.dt) << c
	<< maspect[0] << c << maspect[zeropt] << c << maspect[gs] << c;
    int ind_max = distance(maspect.begin(), max_element(maspect.begin(), maspect.end())) - 1;
    int ind_min = distance(maspect.begin(), min_element(maspect.begin(), maspect.end())) - 1;
    ofs << maspect[ind_max] << c << ind_max << c << p.r[(p.save_pt)*ind_max] << c
	<< pw4(ps[ind_max])*r2(&p,ind_max) << c
	<< maspect[ind_min] << c << ind_min << c << p.r[(p.save_pt)*ind_min] << c
	<< pw4(ps[ind_min])*r2(&p,ind_min) << endl;
  }
  gft_close_all();
  
  cout << outfile << "  written with:" << endl;
  cout << "grid points used = " << npts << "  " << ((p.same_grids) ?
   "(same grids)" : "(dif grids)") << "\ntime steps used = " << num_steps
       << "  " << ((p.same_times) ? "(same times)" : "(dif times)") << endl;
  ofs.close();

  return 0;
}
