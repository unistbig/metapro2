/*
 * wZ.cpp
 *
 *  Created on: 2019. 6. 25.
 *      Author: node02
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <math.h>

using namespace std;

//wZ(id, pmat, dmat, nP, nFeat, direc, output);
void wZ(vector<string> id, double ** pmat, double ** dmat,  int tt, vector<int> SS_z, int nFeat_z, int nP, string output_prefix)
{
	int i, j;
	vector<double> P, Z, D;
	double num, sum_d, finalP;
	ofstream OUT;
	string output = output_prefix + "_wZ.out";
	OUT.open(output.c_str());

	sum_d = 0.0;
	for(i = 0; i < nP; i++)
	{
		sum_d+=pow(SS_z[i],2);
	}
	sum_d = sqrt(sum_d);
	for(i = 0; i < nFeat_z; i++)
	{
		P.clear();
		for(j = 0;  j < nP ; j++)
			P.push_back(pmat[i][j]);
		if(tt == 1)
		{
			double num_temp;
			num = 0.0;
			for(j = 0; j < nP; j++)
			{
				num_temp = gsl_cdf_gaussian_Qinv(P[j], 1.0) * SS_z[j];
				num += num_temp;
			}
			num /= sum_d;
			finalP = gsl_cdf_gaussian_Q(num, 1.0);
			OUT << id[i] << "\t" << finalP << "\n";
		}else{
			D.clear();
			for(j = 0; j < nP; j++)
				D.push_back(dmat[i][j]);
			num = 0.0;
			double num_temp;
			for(j = 0; j < nP; j++)
			{
				num_temp = gsl_cdf_gaussian_Qinv(D[j]*P[j]/2.0,1.0) * SS_z[j];
				num += num_temp;
			}
			num /= sum_d;
			finalP = 2.0*gsl_cdf_gaussian_Q(abs((-1)*num),1.0);
			OUT << id[i] << "\t" << finalP << "\n";
		}
	}
	OUT.close();
}
