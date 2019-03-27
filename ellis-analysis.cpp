/*
  DIAGNOSTIC CODES:
  xi = ResXi
  pi = ResPi
  xi2 = ResXi2
  pi2 = ResPi2
  al = ResAl
  be = ResBe
  ps = ResPs
  ixi = iresXi
  ipi = iresPi
  ixi2 = iresXi2
  ipi2 = iresPi2
  ial = iresAl
  ibe = iresBe
  ips = iresPs
  mass = maspect
  null = outnull
  revnull = outnull_rev
  ric = ricci
  ham = hamiltonian = cHamiltonian
  mom = momentum = cMomentum
  kext = cKext
  dkext = cDtKext
  tt = EEtt
  tx = xt = EEtx
  xx = EExx
  hh = EEhh

 */

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
  FLDS f;
  WRS wr;
  gft_set_multi();
  time_t start_time = time(NULL); // time for rough performance measure
  // INITIALIZE PARAMETERS AND READERS/WRITERS
  vector<BBHP *> writer_vec = analysis_init(&wr, &f, &p, argc, argv);
  if (writer_vec.size() == 0) {
    cout << "\nANALYSIS INIT ERROR" << endl;
    return -1;
  }
  else { cout << "\nanalysis in progress..." << endl; }
  vector<BBHP *> reader_vec { &(wr.p_Al), &(wr.p_Be), &(wr.p_Ps),
      &(wr.p_Xi), &(wr.p_Pi) };
  if (p.write_xp2) {
    reader_vec.push_back(&(wr.p_Xi2));
    reader_vec.push_back(&(wr.p_Pi2));
  }
  SSV s;
  s.lsq = p.lsq;
  p.t = 0;
  for (int t = 1; t < ((p.nsteps) + 2); ++t) {
    // SET OLD FIELDS
    wr.set_old_fields(&f);
    // GET NEW TIME DATA
    if (read_bbhp_vec(reader_vec, t) == 0) { return t; }
    // SET BOUNDARY VALS and COMPUTE
    wr.get_bound_vals(s, &f, &p, 0);
    compute_bbhp_vec(writer_vec, s, 0);
    for (int k = 1; k < p.lastwr; ++k) {
      // SET SITE VALS and COMPUTE
      wr.get_site_vals(s, &f, &p, k);
      compute_bbhp_vec(writer_vec, s, k);
    }
    // SET BOUNDARY VALS and COMPUTE
    wr.get_bound_vals(s, &f, &p, p.lastwr);
    compute_bbhp_vec(writer_vec, s, p.lastwr);
    // WRITE DIAGNOSTICS and INCREMENT TIME
    write_bbhp_vec(writer_vec, &p);
    p.t += p.dt;
  }
  
  gft_close_all();
  cout << (p.outfile) + " analysis completed in "
       << difftime(time(NULL),start_time) << " seconds" << endl;
  
  return 0;
}
