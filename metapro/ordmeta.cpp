/*
 * ordmeta.cpp
 *
 *  Created on: 2019. 6. 26.
 *      Author: node02
 */
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
#include <iterator>

using namespace std;
using std::setw;
using std::setprecision;
using std::numeric_limits;
using namespace GiNaC;

void ordmeta(vector<string> id, double ** pmat, double ** dmat, int tt,  int nFeat, int nP, string output_prefix)
{
	ofstream OUTPUT;
	double LB, LB2, baseP, baseP2, Pval, Pval2;
	vector<string> CohortVEC;
	symbol x("x"), y("y");
	double temp;
	ex templete, templete2, formula, formula2;
	int boost_error;
	int datasize;
	string TEMP_NAidx;
	vector<string> NAidx;
	int i,j;
	vector<double> P, Pinv;
	int nP2;
	double minP, minP2;
	double P_temp, P_temp2;
	string numCohort, numCohort2, nc;
	string output = output_prefix + "_ORDMETA.out";
	OUTPUT.open(output.c_str());

	double inf = std::numeric_limits<double>::infinity();
	for(i = 0; i < nFeat; i++)
	{
		TEMP_NAidx = "";
		boost_error = 0;
		Pval = 0.0;
		Pval2 = 0.0;
		datasize = 0;
		P.clear(); Pinv.clear();
		P.resize(nP,1.0); Pinv.resize(nP, 1.0);
		for(j = 0; j<nP; j++)
		{
			if(!isnan(dmat[i][j])){datasize++;}

			if(tt == 1)
				P[j] = pmat[i][j];
			else
			{
				if(dmat[i][j] > 0){
					P[j] = pmat[i][j]/2.0;
				}else if( dmat[i][j] < 0 ){
					P[j] = 1-pmat[i][j]/2.0;
				}else if(dmat[i][j]==0){P[j] = 0.5;
				}else{TEMP_NAidx = TEMP_NAidx +" "+to_string(j+1); temp=inf;}
			}
		}
		NAidx.push_back(TEMP_NAidx);

		vector<double> Psort;

		Psort.assign(P.begin(), P.end());
		cout << "C" << endl;
		sort(Psort.begin(), Psort.end());
		Psort.resize(datasize);

		for(j = 0; j<nP; j++)
			Pinv[j] = 1.0-P[j];

		vector<double> Pinv_sort;
		Pinv_sort.assign(Pinv.begin(), Pinv.end());
		sort(Pinv_sort.begin(), Pinv_sort.end());
		Pinv_sort.resize(datasize);

		nP2 = Psort.size();
		minP = 1.0;
		minP2 = 1.0;

		for(j = 0; j < nP2; j++)
		{
			P_temp = boost::math::ibeta(j+1, nP2-j, Psort[j]);
			if(P_temp < minP){minP = P_temp; baseP = Psort[j];}

			P_temp2 = boost::math::ibeta(j+1, nP2-j, Pinv_sort[j]);
			if(P_temp2 < minP2){minP2 = P_temp2; baseP2 = Pinv_sort[j];}
		}
		numCohort = "";
		numCohort2 = "";
		for(j = 0; j < nP; j++)
		{
			if(!isinf(P[j]) && P[j] <= baseP)
				numCohort = numCohort + to_string(j+1) + " " ;
			if(!isinf(Pinv[j]) && Pinv[j] <= baseP)
				numCohort2 = numCohort2 + to_string(j+1)+ " ";
		}

		formula = 1;
		formula2 = 1;

		for(j = 0; j < (nP2-1); j++)
		{
			LB = gsl_cdf_beta_Pinv(minP, j+1, nP2 - j);
			LB2 = gsl_cdf_beta_Pinv(minP2, j+1, nP2-j);
			if(isnan(LB) || isnan(LB2))
			{
				boost_error = 1;
				break;
			}

			formula = formula.subs(y==x);
			templete = integral(x, LB, y, formula);
			formula = templete.eval_integ();
			formula = formula*(j+1);
			formula2 = formula2.subs(y==x);
			templete2 = integral(x, LB2, y, formula2);
			formula2 = templete2.eval_integ();
			formula2 = formula2*(j+1);
		}
		if(boost_error == 0)
		{
			LB = gsl_cdf_beta_Pinv(minP, nP,1); //Fi_inv(nStudy, nStudy, minP);
			LB2 = gsl_cdf_beta_Pinv(minP2, nP, 1);
			if(isnan(LB)){
				CohortVEC.push_back(" ");
				boost_error = 1;
				continue;
			}else{
				formula = formula.subs(y==x);
				templete=integral(x, LB, 1, formula);
				formula = templete.eval_integ();
				formula = formula*nP2;
				double D = ex_to<numeric>(formula).to_double();
				Pval+=D;
				//PVEC.push_back(1-Pval);

				formula2 = formula2.subs(y==x);
				templete2 = integral(x, LB2, 1, formula2);
				formula2 = templete2.eval_integ();
				formula2 = formula2*nP2;
				double D2 = ex_to<numeric>(formula2).to_double();
				Pval2+=D2;

//				cout << "Pval =" << Pval << endl;
//				cout << "Pval2 =" << Pval2 << endl;
				//PVEC2.push_back(1-Pval2);
				double finalP;
				if(tt == 2)
				{
					if(Pval>Pval2){finalP = 2*(1.0-Pval); nc = numCohort;}
					if(Pval2>=Pval){finalP = 2*(1.0-Pval2);nc = numCohort2;}
				}else{
					if(Pval>Pval2){finalP = (1.0-Pval);  nc = numCohort;}
					if(Pval2>=Pval){finalP = (1.0-Pval2); nc = numCohort;}
				}
				if(finalP>1.0){finalP = 1.0;}
				OUTPUT << id[i] << "\t" << finalP << "\t" << nc << "\n";
			}
		}
	}
	OUTPUT.close();
}
