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
  if ((argc < 3) || (argc > 5)) {
    cout << "\nINCORRECT NUMBER OF ARGUMENTS" << endl;
    return argc;
  }
  else {
    int new_argc = 2;
    char *new_argv[] = {argv[0], argv[1]};
    int err_code = params_init(&p, new_argc, new_argv);
    if (err_code != 0) {
      cout << "\nPARAM INIT error code = " << err_code << endl;
      return err_code;
    }
  }
  str c = ",";
  str outfile = p.outfile;
  bool write_mass = false;
  bool write_null = false;
  bool write_areal = false;
  bool write_areal_new = false;
  for (int k = 2; k < argc; ++k) {
    str out_type = argv[k];
    if ((out_type == "m") || (out_type == "M") || (out_type == "mass") || (out_type == "Mass")
	|| (out_type == "maspect") || (out_type == "Maspect")) {
      write_mass = true;
    }
    else if ((out_type == "n") || (out_type == "N") || (out_type == "null") || (out_type == "Null")
	|| (out_type == "outnull") || (out_type == "Outnull") || (out_type == "outgoing_null")
	|| (out_type == "revnull") || (out_type == "Revnull") || (out_type == "outnull_rev")) {
      write_null = true;
    }
    else if ((out_type == "a") || (out_type == "A") || (out_type == "areal") || (out_type == "Areal")
	|| (out_type == "r") || (out_type == "R") || (out_type == "rad") || (out_type == "Rad")
	|| (out_type == "radius") || (out_type == "Radius")) {
      write_areal = true;
      write_areal_new = false;
    }
    else if ((out_type == "a_new") || (out_type == "A_new") || (out_type == "areal_new") || (out_type == "Areal_new")
	|| (out_type == "r_new") || (out_type == "R_new") || (out_type == "rad_new") || (out_type == "Rad_new")
	|| (out_type == "radius_new") || (out_type == "Radius_new")) {
      write_areal_new = true;
      write_areal = false;
    }
  }
    
  // derived parameters from coarse file
  int gs = p.lastwr;
  int zeropt = p.zerowr;
  int num_steps = p.nsteps / p.save_step;
  int npts = gs + 1;
  if (zeropt != (gs/2)) {
    cout << "\nERROR WITH GRID SHAPE\n" << endl;
    return -1;
  }
  VD ps(npts, 0.0);
  str outname = outfile;
  vector<dbl *> field_arr {&ps[0]};
  str ps_name = "Ps-" + outfile + ".sdf";
  vector<char *> name_arr {&ps_name[0]};

  gft_set_multi();

  if (write_areal_new) {
    // output file 
    VD areal(npts, 0.0);
    str areal_name = "areal-" + outfile + ".sdf";
    ofstream ofs;
    str outname = "arealCheck-" + outfile + ".csv";
    ofs.open(outname, ofstream::out);
    ofs << p.outfile << "\nstep,time,R(0),R_min,at,x,R_Lmin,at,x,R_Rmin,at,x" << endl;
    // iterate through time steps
    int times[1];
    for (int t = 0; t < (num_steps + 1); ++t) {
      times[0] = t+1;
      if (read_step(name_arr, times, field_arr, 1) == 0) {
	cout << areal_name << "  written up to " << t << endl;
	break;
      }
      for (int j = 0; j < npts; ++j) { areal[j] = sq(ps[j])*sqrt(r2(&p,j)); }
      write_sdf_direct(&areal_name[0], &areal[0], &p);
      ofs << t << c << t*(p.save_step)*(p.dt) << c << areal[zeropt] << c;
      int ind_min = distance(areal.begin(), min_element(areal.begin(), areal.end()));
      int ind_Lmin = distance(areal.begin(), min_element(areal.begin(), areal.begin() + zeropt));
      int ind_Rmin = distance(areal.begin(), min_element(areal.begin() + zeropt + 1, areal.end()));
      ofs << areal[ind_min] << c << ind_min << c << p.r[(p.save_pt)*ind_min] << c
	  << areal[ind_Lmin] << c << ind_Lmin << c << p.r[(p.save_pt)*ind_Lmin] << c
	  << areal[ind_Rmin] << c << ind_Rmin << c << p.r[(p.save_pt)*ind_Rmin] << endl;
      p.t += ((p.dt) * (p.save_step));
    }
    ofs.close();
  }
  else if (write_areal) {
    // output file 
    VD areal(npts, 0.0);
    str areal_name = "areal-" + outfile + ".sdf";
    ofstream ofs;
    str outname = "arealCheck-" + outfile + ".csv";
    ofs.open(outname, ofstream::out);
    ofs << p.outfile << "\nstep,time,R(0),R_min,at,x,R_Lmin,at,x,R_Rmin,at,x" << endl;
    // iterate through time steps
    for (int t = 0; t < (num_steps + 1); ++t) {
      if (gft_read_brief(&areal_name[0], t+1, &areal[0]) == 0)  {
	cout << areal_name << "  written up to " << t << endl;
	break;
      }
      ofs << t << c << t*(p.save_step)*(p.dt) << c << areal[zeropt] << c;
      int ind_min = distance(areal.begin(), min_element(areal.begin(), areal.end()));
      int ind_Lmin = distance(areal.begin(), min_element(areal.begin(), areal.begin() + zeropt));
      int ind_Rmin = distance(areal.begin(), min_element(areal.begin() + zeropt + 1, areal.end()));
      ofs << areal[ind_min] << c << ind_min << c << p.r[(p.save_pt)*ind_min] << c
	  << areal[ind_Lmin] << c << ind_Lmin << c << p.r[(p.save_pt)*ind_Lmin] << c
	  << areal[ind_Rmin] << c << ind_Rmin << c << p.r[(p.save_pt)*ind_Rmin] << endl;
      p.t += ((p.dt) * (p.save_step));
    }
    ofs.close();
  }

  if (write_mass) {
    VD maspect(npts, 0.0);
    field_arr.push_back(&maspect[0]);
    str maspect_name = "maspect-" + outfile + ".sdf";
    name_arr.push_back(&maspect_name[0]);
    // output file
    ofstream ofs;
    str outname = "massCheck-" + outfile + ".csv";
    ofs.open(outname, ofstream::out);
    ofs << p.outfile << "\nstep,time,m(-R),m(0),m(R),m_max,at,x,r2ps4,m_min,at,x,r2ps4" << endl;  
    // iterate through time steps
    int times[2];
    for (int t = 2; t < (num_steps + 1); ++t) {
      times[0] = t+1; times[1] = t+1;
      if (read_step(name_arr, times, field_arr, 2) == 0) {
	cout << outname << "  written up to " << t << endl;
        break;
      }
      ofs << t << c << t*(p.save_step)*(p.dt) << c
	  << maspect[0] << c << maspect[zeropt] << c << maspect[gs] << c;
      int ind_max = distance(maspect.begin(), max_element(maspect.begin(), maspect.end()));
      int ind_min = distance(maspect.begin(), min_element(maspect.begin(), maspect.end()));
      ofs << maspect[ind_max] << c << ind_max << c << p.r[(p.save_pt)*ind_max] << c
	  << pw4(ps[ind_max])*r2(&p,ind_max) << c
	  << maspect[ind_min] << c << ind_min << c << p.r[(p.save_pt)*ind_min] << c
	  << pw4(ps[ind_min])*r2(&p,ind_min) << endl;
    }
    ofs.close();
  }
  
  if (write_null) {
    VD outnull(npts, 0.0);
    VD revnull(npts, 0.0);
    field_arr.push_back(&outnull[0]);
    field_arr.push_back(&revnull[0]);
    str outnull_name = "outnull-" + outfile + ".sdf";
    str revnull_name = "revnull-" + outfile + ".sdf";
    name_arr.push_back(&outnull_name[0]);
    name_arr.push_back(&revnull_name[0]);
    // output file
    ofstream ofs;
    str outname = "nullCheck-" + outfile + ".csv";
    ofs.open(outname, ofstream::out);
    ofs << p.outfile << "\nstep,time,point,x,r2ps4,outnull,revnull" << endl;    
    // iterate through time steps
    int times[3];
    for (int t = 2; t < (num_steps + 1); ++t) {
      times[0] = t+1; times[1] = t+1; times[2] = t+1;
      if (read_step(name_arr, times, field_arr, 3) == 0) {
	cout << outname << "  written up to " << t << endl;
        break;
      }
      ofs << t << c << t*(p.save_step)*(p.dt) << c << zeropt << ",0,"
	  << pw4(ps[zeropt])*r2(&p,zeropt) << c
	  << outnull[zeropt] << c << revnull[zeropt] << endl;
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
    ofs.close();
  }

  gft_close_all();
  str out_types = " ";
  if (write_mass) { out_types += "mass "; }
  if (write_null) { out_types += "null "; }
  if (write_areal) { out_types += "areal "; }
  cout << outfile << out_types << "written with:\ngrid points used = " << npts
       << "  "  << ((p.same_grids) ? "(same grids)" : "(dif grids)")
       << "\ntime steps used = " << num_steps
       << "  " << ((p.same_times) ? "(same times)" : "(dif times)") << endl;

  return 0;
}
