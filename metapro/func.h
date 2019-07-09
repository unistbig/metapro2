/*
 * func.h
 *
 *  Created on: 2019. 6. 26.
 *      Author: node02
 */


#include <time.h>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include<fstream>
#include<vector>
#include<stdlib.h> /* exit, EXIT_FAILURE*/
#include <stdio.h>
#include<string>
#include<string.h>
#include<cstring>
#include<sstream>
#include <iomanip>
#include <math.h> /*erf*/
#include<cmath>
#include<iomanip>
#include <stdexcept>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <ginac/ginac.h> // Symbolic integration
#include <boost/algorithm/string/replace.hpp>
#include <limits>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <map>
using namespace std;

#ifndef ARTP_H
#define ARTP_H
void ARTP(vector<string> id, double ** pmat, double ** dmat, int tt, int nFeat, int nP, int nperm, string output_prefix);
#endif

#ifndef LANCASTER_H
#define LANCASTER_H
void lancaster(vector<string> id, double ** pmat, double ** dmat, int tt,  vector<int> SS, int nFeat, int nP, string output_prefix);
#endif

#ifndef ORDMETA_H
#define ORDMETA_H
void ordmeta(vector<string> id, double ** pmat, double ** dmat, int tt,  int nFeat, int nP, string output_prefix);
#endif

#ifndef WFISHER_H
#define wFISHER_H
void wFisher(vector<string> id, double ** pmat, double ** dmat,int tt,   vector<int> SS, int nFeat, int nP, string output_prefix);
#endif

#ifndef WZ_H
#define WZ_H
void wZ(vector<string> id, double ** pmat, double ** dmat, int tt,  vector<int> SS, int nFeat, int nP, string output_prefix);
#endif
