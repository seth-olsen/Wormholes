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
  time_t start_time = time(NULL); // time for rough performance measure
  // INITITALIZATION
  PAR p;
  FLDS f;
  int err_code = params_init(&p, argc, argv);
  if (err_code != 0) {
    cout << "\nPARAM INIT error code = " << err_code << endl;
    return err_code;
  }
  err_code = fields_init(&f, &p);
  if (err_code != 0) {
    cout << "\nFIELD INIT error code = " << err_code << endl;
    return err_code;
  }
  else { cout << "\nsimulation in progress..." << endl; }
  gft_set_multi();
  if ((p.same_times) || (p.same_grids)) {
    WRS wr;
    vector<BBHP *> writer_vec = writers_init(&wr, &f, &p);
    // DO SOME CHECKS HERE
    for (int i = 0; i < (p.nsteps); ++i) {
      // WRITING
      if ((i % (p.save_step)) == 0) {
	write_bbhp_vec(writer_vec, &p);
	if ((i % (p.check_diagnostics)) == 0) {
	  write_diagnostics(&wr, &f, &p);
	}
      }
      // SOLVE FOR NEXT STEP
      err_code = fields_step(&f, &p, i);
      if (err_code) {
	cout << "\nFIELD STEP error code = " << err_code << endl;
	gft_close_all();
	return err_code;
      }
    }
    // WRITE LAST STEP
    if (((p.nsteps) % (p.save_step)) == 0) {
      write_bbhp_vec(writer_vec, &p);
      if (((p.nsteps) % (p.check_diagnostics)) == 0) {
	write_diagnostics(&wr, &f, &p);
      }
    }
  }
  else {
    str al_nm = "Al-" + (p.outfile) + ".sdf";
    str be_nm = "Be-" + (p.outfile) + ".sdf";
    str ps_nm = "Ps-" + (p.outfile) + ".sdf";
    str xi_nm = "Xi-" + (p.outfile) + ".sdf";
    str pi_nm = "Pi-" + (p.outfile) + ".sdf";
    for (int i = 0; i < (p.nsteps); ++i) {
      // WRITING
      gft_out_bbox(&(al_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Al[0]));
      gft_out_bbox(&(be_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Be[0]));
      gft_out_bbox(&(ps_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Ps[0]));
      gft_out_bbox(&(xi_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Xi[0]));
      gft_out_bbox(&(pi_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Pi[0]));
      // SOLVE FOR NEXT STEP
      err_code = fields_step(&f, &p, i);
      if (err_code) {
	cout << "\nFIELD STEP error code = " << err_code << endl;
	gft_close_all();
	return err_code;
      }
    }
    // WRITE LAST STEP
    gft_out_bbox(&(al_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Al[0]));
    gft_out_bbox(&(be_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Be[0]));
    gft_out_bbox(&(ps_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Ps[0]));
    gft_out_bbox(&(xi_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Xi[0]));
    gft_out_bbox(&(pi_nm[0]), (p.t), &(p.npts), 1, &(p.coord_lims[0]), &(f.Pi[0]));
  }
  gft_close_all();
  cout << (p.outfile) + " written in "
       << difftime(time(NULL),start_time) << " seconds" << endl;
  return err_code;
}
  
