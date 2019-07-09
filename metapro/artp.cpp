/*
 * artp.cpp
 *
 *  Created on: 2019. 6. 26.
 *      Author: node02
 */

/*
 * ARTP.cpp
 *
 *  Created on: 2019. 5. 29.
 *      Author: node02
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

vector<double> sequentialProduct(vector<double> v, int size)
{
	vector<double> res(1,size);
	unsigned int i;
	res[0] = v[0];
	for(i = 1; i < v.size(); i++)
		res.push_back(v[i]*res[i-1]);
	return res;
}

void ARTP(vector<string> id, double ** pmat, double ** dmat, int tt, int nFeat, int nP, int nperm, string output_prefix)
{
	std::srand (  std::time(0)  );
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	vector<double> pvec, rvec, wvec, svec;
	vector<double> pvec_pos, pvec_neg, wvec_pos, wvec_neg, svec_pos, svec_neg;
	vector<int> D;

	double ** wmat, ** wmat_b;
	double ** smat, ** smat_b;
	double number;
	double svec_min_pos, svec_min_neg;
	int nperm_basic = 1000;
	int f,i,j,k;
	double obsMinP, finalP;
	ofstream OUT;
	string output = output_prefix + "_ARTP.out";
	OUT.open(output.c_str());

	wmat  = (double **)malloc(sizeof(double *) * nperm);
	for(i = 0; i < nperm; i++)
		wmat[i] = (double *)malloc(sizeof(double) * nP);

	smat  = (double **)malloc( sizeof(double *) *nperm);
	for(i = 0; i < nperm; i++)
		smat[i] = (double *)malloc(sizeof(double) * nP);

	wmat_b  = (double **)malloc(sizeof(double *) * nperm_basic);
	for(i = 0; i < nperm_basic; i++)
		wmat_b[i] = (double *)malloc(sizeof(double) * nP);

	smat_b  = (double **)malloc(sizeof(double *) * nperm_basic);
	for(i = 0; i < nperm_basic; i++)
		smat_b[i] = (double *)malloc(sizeof(double) * nP);

	vector<int> rank_ref(nperm);
	vector<int> rank_ref_basic(nperm_basic);

	for(i = 0; i < nperm; i++)
	{
		rvec.clear();
		for(j = 0; j < nP; j++)
		{
			number = distribution(generator);
			rvec.push_back(number);
		}
		sort(rvec.begin(), rvec.end());

		for(j = 0; j < nP; j++)
		{
			if(j == 0){wmat[i][j] = rvec[j];
			}else{
				wmat[i][j] = rvec[j] * wmat[i][j-1];
			}
		}
	}
	cout << "A" << endl;
	for(i = 0; i < nperm_basic; i++)
	{
		rvec.clear();
		for(j = 0; j < nP; j++)
		{
			number = distribution(generator);
			rvec.push_back(number);
		}
		sort(rvec.begin(), rvec.end());
		vector<double> a;
		a = sequentialProduct(rvec, nP);
		for(j = 0; j < nP; j++)
			wmat_b[i][j] = a[j];
	}
	cout << "B" << endl;
	for(i = 0; i < nP; i++)
	{
		rank_ref_basic.clear();
		rank_ref_basic.resize(nperm_basic,1);
		for(j = 0; j < (nperm_basic-1); j++)
		{
			for(k = (j+1); k < nperm_basic; k++)
			{
				if(wmat_b[j][i] >= wmat_b[k][i])
				{
					rank_ref_basic[j]++;
				}else if(wmat_b[j][i] <= wmat_b[k][i])
				{
					rank_ref_basic[k]++;
				}
			}
		}
		for(j = 0; j < nperm_basic; j++)
			smat_b[j][i] = ((double) rank_ref_basic[j] / (double) nperm_basic);
	}

	for(i = 0; i < nP; i++)
	{
		rank_ref.clear();
		rank_ref.resize(nperm,1);
		for(j = 0; j < (nperm-1); j++)
		{
			for(k = (j+1); k < nperm; k++)
			{
				if(wmat[j][i] >= wmat[k][i])
				{
					rank_ref[j]++;
				}else if(wmat[j][i] <= wmat[k][i])
				{
					rank_ref[k]++;
				}
			}
		}
		for(j = 0; j < nperm; j++)
			smat[j][i] = (double) rank_ref[j] / (double) nperm;
	}

	vector<double> minP_basic;
	minP_basic.clear();
	for(i = 0; i < nperm_basic; i++)
	{
		double minv=1.0 ;
		for(j = 0; j < nP; j++)
			if(smat_b[i][j]<minv){minv = smat_b[i][j];}
		minP_basic.push_back(minv);
	}

	vector<double> minP;
	minP.clear();
	for(i = 0; i < nperm; i++)
	{
		double minv=1.0 ;
		for(j = 0; j < nP; j++)
			if(smat[i][j]<minv){minv = smat[i][j];}
		minP.push_back(minv);
	}

	cout << "nFeat: " << nFeat << endl;

	for(f = 0; f < nFeat; f++)
	{
		if(tt == 2)
		{
			pvec_pos.clear(); pvec_neg.clear();
			wvec_pos.clear(); wvec_neg.clear();
			for(j = 0; j < nP; j++)
			{
				if(dmat[f][j] > 0)
				{
					pvec_pos.push_back(pmat[f][j]/2.0);
					pvec_neg.push_back(1-pmat[f][j]/2.0);
				}else{
					pvec_pos.push_back(1-pmat[f][j]/2.0);
					pvec_neg.push_back(pmat[f][j]/2.0);
				}
			}
			sort(pvec_pos.begin(), pvec_pos.end());
			sort(pvec_neg.begin(), pvec_neg.end());

			wvec_pos = sequentialProduct(pvec_pos,nP);
			wvec_neg = sequentialProduct(pvec_neg,nP);

			svec_pos.clear(); svec_neg.clear();
			cout <<"C"<<endl;
			svec_pos.resize(nP, 1);
			svec_neg.resize(nP, 1);
			cout <<"D"<<endl;
			for(j = 0; j < nP; j++)
			{
				for(i = 0; i < nperm_basic; i++)
				{
					if(wvec_pos[j] >= wmat_b[i][j])
						svec_pos[j]++;
					if(wvec_neg[j] >= wmat_b[i][j])
						svec_neg[j]++;
				}
				svec_pos[j]/=(nperm_basic+1);
				svec_neg[j]/=(nperm_basic+1);
			}
			svec_min_pos = *min_element(svec_pos.begin(), svec_pos.end());
			svec_min_neg = *min_element(svec_neg.begin(), svec_neg.end());
			if(svec_min_pos <= svec_min_neg){svec = svec_pos;}else{svec = svec_neg;}

		}else{
			pvec.clear();
			wvec.clear();
			for(j = 0; j < nP; j++)
				pvec.push_back(pmat[f][j]);
			sort(pvec.begin(), pvec.end());
			wvec = sequentialProduct(pvec,nP);

			svec.clear();
			svec.resize(nP, 1);
			for(j = 0; j < nP; j++)
			{
				for(i = 0; i < nperm_basic; i++)
				{
					if(wvec[j] >= wmat_b[i][j])
						svec[j]++;
				}
				svec[j]/=(nperm_basic+1);
			}
		}

/*

		svec.clear();
		svec.resize(nP, 1);
		for(j = 0; j < nP; j++)
		{
			for(i = 0; i < nperm_basic; i++)
			{
				if(wvec[j] >= wmat_b[i][j])
					svec[j]++;
			}
			svec[j]/=(nperm_basic+1);
		}
*/

		obsMinP = 1.0;
		for(j = 0; j< nP; j++)
			if(obsMinP > svec[j])
				obsMinP = svec[j];

		// Significance
		finalP = 0.0;
		for(i = 0; i < nperm_basic ; i++)
			if(obsMinP > minP_basic[i])
				finalP+=1.0;
		finalP/=(nperm_basic+1);
		if(tt == 2){finalP *= 2.0;}

		if(finalP > (1.0/ (double) nperm_basic)){
			OUT << id[f] << "\t" << finalP << endl;
		}else{
			svec.clear();
			svec.resize(nP,1);
			for(j = 0; j < nP; j++)
			{
				for(i = 0; i < nperm; i++)
				{
					if(wvec[j] >= wmat[i][j])
						svec[j]++;
				}
				svec[j]/=(nperm+1);
			}

			obsMinP = 1.0;
			for(j = 0; j< nP; j++)
				if(obsMinP > svec[j])
					obsMinP = svec[j];

			// Significance
			finalP = 0.0;
			for(i = 0; i < nperm ; i++)
				if(obsMinP > minP[i])
					finalP+=1.0;
			finalP/=(nperm+1);
			if(tt == 2){finalP *= 2.0;}
			OUT << id[f] << "\t" << finalP << endl;
		}

	}
	OUT.close();

	for(i = 0; i < nFeat; i++)
		free(pmat[i]);
	free(pmat);

	for(i = 0; i < nperm; i++)
		free(wmat[i]);
	free(wmat);

	for(i = 0; i < nperm_basic; i++)
		free(wmat_b[i]);
	free(wmat_b);

	for(i = 0; i < nperm; i++)
		free(smat[i]);
	free(smat);

	for(i = 0; i < nperm_basic; i++)
		free(smat_b[i]);
	free(smat_b);

}
