#ifndef SIMULATION_HEADER_H
#define SIMULATION_HEADER_H

// libraries
#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <cmath> // for ICs
#include <vector> // for everything
#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include <iterator> // for distance()
#include "bbhutil.h" // for output to .sdf
#include "lapacke.h"

using namespace std;

// TYPEDEFS
typedef string str;
typedef double dbl;

typedef vector<dbl> VD;
typedef map<int, double> MAPID;
typedef map<int, double> MAPII;

typedef struct bbhutil_params BBHP;
typedef struct sim_fields FLDS;
typedef struct sim_writers WRS;
typedef struct sim_params PAR;
typedef struct sim_site_values SSV;

typedef dbl (*D_FN)(SSV&);
typedef void (*SSVSET)(SSV&, FLDS *, PAR *, int);


#endif
