/*
 * lancaster.cpp
 *
 *  Created on: 2019. 6. 26.
 *      Author: node02
 */

/*
 * wFisher.cpp
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

void lancaster(vector<string> id, double ** pmat, double ** dmat, int tt, vector<int> SS, int nFeat, int nP, string output_prefix)
{
	int i,j;
	double num, numpos, numneg;
	double Gpos, Gneg, Gsum_pos, Gsum_neg, Pval_pos, Pval_neg, finalP;
	double G, Gsum;
	int sum;
	vector<double> P, PPOS, PNEG;
	vector <int> D;
	vector <double> weight;
	ofstream OUT;
	string output = output_prefix + "_Lancaster_SS.out";
	OUT.open(output.c_str());

	sum = 0;
	for(i = 0; i < nP; i++)
		sum+=SS[i];
	weight.clear();
	for(i = 0; i < nP; i++)
	{
		num = (double) SS[i];
		weight.push_back(num);
	}

	for(i = 0; i < nFeat; i++)
	{
		P.clear();

		for(j = 0; j < nP; j++)
			P.push_back(pmat[i][j]);
		if(tt == 2)
		{
			D.clear();
			for(j = 0; j < nP; j++)
				D.push_back(dmat[i][j]);

			PPOS.clear();
			PNEG.clear();
			Gsum_pos = 0.0;
			Gsum_neg = 0.0;
			for(j = 0; j < nP; j++)
			{
				if(D[j] >= 0.0){
					PPOS.push_back(P[2]/2.0);
					PNEG.push_back(1-P[2]/2.0);
					Gpos = gsl_cdf_gamma_Qinv(PPOS[j], weight[j], 2);
					Gneg = gsl_cdf_gamma_Qinv(PNEG[j], weight[j], 2);
					Gsum_pos += Gpos;
					Gsum_neg += Gneg;
				}else if(D[j] < 0.0){
					PPOS.push_back(1-P[2]/2.0);
					PNEG.push_back(P[2]/2.0);
					Gpos = gsl_cdf_gamma_Qinv(PPOS[j], weight[j], 2);
					Gneg = gsl_cdf_gamma_Qinv(PNEG[j], weight[j], 2);
					Gsum_pos += Gpos;
					Gsum_neg += Gneg;
				}
			}
			Pval_pos = gsl_cdf_gamma_Q(Gsum_pos, sum, 2);
			Pval_neg = gsl_cdf_gamma_Q(Gsum_neg, sum, 2);
			if(Pval_pos <= Pval_neg){
				finalP = 2.0*Pval_pos;
			}else
				finalP = 2.0*Pval_neg;
		}else{
			Gsum = 0.0;
			for(j = 0; j < nP; j++)
			{
				G = gsl_cdf_gamma_Qinv(P[j], weight[j],2);
				Gsum+=G;
			}
			finalP = gsl_cdf_gamma_Q(Gsum, sum, 2);
		}
		// write result
		OUT << id[i] << "\t" << finalP << "\n";
	}
	OUT.close();
}



