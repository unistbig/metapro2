/*
 * main.cpp
 *
 *  Created on: 2019. 6. 25.
 *      Author: Sora Yoon
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
#include "func.h"

using namespace std;

string buf;
vector<string> tokens;
string pvalue, direc, method, samplesize, output, perm;
int testtype;

void split(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	string::size_type pos = str.find_first_of(delimiters, lastPos);
	while(string::npos!=pos || string::npos != lastPos){
		tokens.push_back(str.substr(lastPos, pos-lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);
		pos = str.find_first_of(delimiters, lastPos);
	}
}

void PrintHelp()
{
	std::cout <<
			"--method <m>	P-value combination method. One of wZ, wFisher, ordmeta, lancaster and artp. Required.\n"
			"--pvalue <p>	Input P-value matrix file. Required.\n"
			"--output <o>	Output file name. The directory must exist. Required.\n"
			"--testtype <t>	1 for one-tailed test and 2 for two-tailed test. Required. \n"
			"--direc <d>	Effect direction matrix file. Required if --testtype is 2 \n"
			"--samplesize <s>	Sample size file. Optional for wZ/wFisher and required for lancaster. \n"
			"--perm <r>		The number of permutation for ARTP. Default = 10000. Optional.\n"
			"--help <h>";
	exit(1);
}

static int verbose_flag;

void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = ":m:p:t:d:o:s:r:h";
	const option long_opts[]={
			{"verbose",no_argument, &verbose_flag,1},
			{"brief", no_argument, &verbose_flag, 0},
			{"method", required_argument, 0, 'm'},
			{"pvalue", required_argument, 0, 'p'},
			{"testtype", required_argument, 0, 't'},
			{"direc", required_argument, 0, 'd'},
			{"output", required_argument, 0, 'o'},
			{"samplesize", required_argument, 0, 's'},
			{"perm", optional_argument, 0, 'r'},
			{"help", no_argument, 0, 'h'},
			{0,0,0,0}
	};

	int option_index = 0;

	while(true)
	{
		const auto opt = getopt_long(argc, argv, short_opts, long_opts, &option_index);
		if(-1 == opt)
			break;

		switch(opt)
		{
			case 0:
				if(long_opts[option_index].flag != 0)
					break;
				printf("option %s", long_opts[option_index].name);
				if(optarg)
					printf(" with arg %s", optarg);
				printf("\n");
				break;

			case 'm':
				method = string(optarg);
				std::cout << "Meta-analysis method: " << method << endl;
				break;

			case 'p':
				pvalue = string(optarg);
				std::cout << "P-value file: " <<pvalue << endl;
				break;
			case 't':
				testtype = stoi(optarg);
				cout << testtype << endl;
				if(testtype == 1){cout << "One-tailed test" << endl;
				}else if(testtype == 2){
					cout << "Two-tailed test" << endl;
				}else{
					cout << "Invalid test type." << endl;
					exit(0);
				}
				break;
			case 'd':
				direc = string(optarg);
				std::cout<< "Effect size file: " << direc << endl;
				break;
			case 's':
				samplesize = string(optarg);
				std::cout<< "Sample size file: " << samplesize << endl;
				break;
			case 'o':
				output = string(optarg);
				std::cout<< "Output file " << output << endl;
				break;
			case 'r':
				if(optarg){
					perm = string(optarg);
				}
				else
				{
					perm = "10000";
				}
				break;
			case 'h':
				PrintHelp();
				break;
			case '?':
				PrintHelp();
				break;
		}
	}
}

/*long long GetTimeDiff(unsigned int nFlag)
{
  const long long NANOS = 1000000000LL;
  static struct timespec startTS, endTS;
  static long long retDiff = 0;
  if(nFlag == 0)
  {
    retDiff = 0;
    if( -1 == clock_gettime(CLOCK_MONOTONIC, &startTS))
      printf("Failed to call clock_gettime\n");
  }
  else{
    if(-1 == clock_gettime(CLOCK_MONOTONIC, &endTS))
      printf("Failed to call clock_gettime\n");
     retDiff = NANOS * (endTS.tv_sec - startTS.tv_sec) + (endTS.tv_nsec - startTS.tv_nsec);
  }
  return retDiff/1000000;
}*/

int main(int argc, char** argv)
{
	ProcessArgs(argc, argv);

	ifstream IN;
	ofstream OUT;
	vector<string> id;
	double ** pmat;
	int nFeat, nP; // The number of genes and p-values, respectively.
	int i,j;
	string sep;
	vector<string> tokens;
	string buf;
	int tt = testtype;
	if (argc<5) {
	  printf("Mandatory argument(s) missing\n\n");
	  PrintHelp();
	  exit(1);
	}

	nFeat = 0;
	IN.open(pvalue.c_str());
	while(getline(IN, buf))
	{

		nFeat++;
	}
	IN.close();

	IN.open(pvalue.c_str());
	getline(IN, buf);
	tokens.clear();
	split(buf, tokens, "\t");
	sep = "\t";
	if(tokens.size() == 1){
		split(buf, tokens);
		sep = " ";
		if(tokens.size() == 1)
		{
			split(buf, tokens, ",");
			sep = ",";
			if(tokens.size() == 1)
			{
				cout<< "Invalid p-value separator." << endl;
				exit(0);
			}
		}
	}
	nP = tokens.size()-1; // Eliminate gene name column
	IN.close();

	pmat  = (double **)malloc(sizeof(double *) * nFeat);
	for(i = 0; i < nFeat; i++)
		pmat[i] = (double *)malloc(sizeof(double) * nP);

	double ** dmat;
	dmat  = (double **)malloc(sizeof(double *) * nFeat);
	for(i = 0; i < nFeat; i++)
		dmat[i] = (double *)malloc(sizeof(double) * nP);

	IN.open(pvalue.c_str());
	for(i = 0 ; i < nFeat; i++)
	{
		getline(IN, buf);
		tokens.clear();
		split(buf,tokens,sep);
		id.push_back(tokens[0]);
		for(j = 1; j < (nP+1); j++)
			pmat[i][(j-1)] = stod(tokens[j]);
	}
	IN.close();


	if(tt == 2)
	{
		ifstream IN2;
		IN2.open(direc.c_str());
		for(i = 0; i < nFeat; i++)
		{
			getline(IN2,buf);
			tokens.clear();
			split(buf, tokens, sep);
			for(j = 1; j<(nP-1); j++)
				dmat[i][j-1] = stod(tokens[j]);
		}
		IN2.close();
	}

	vector<int> SS(nP);
	ifstream IN3;
	IN3.open(samplesize.c_str());
	for(i = 0; i < nFeat; i++)
	{
		getline(IN3, buf);
		SS[i] = stoi(buf);
	}
	IN3.close();

	if(method == "wZ")
	{
		wZ(id, pmat, dmat, tt, SS, nFeat, nP, output);
	}else if(method == "wFisher")
	{
		//run wFisher
		wFisher(id, pmat, dmat, tt, SS, nFeat, nP, output);
	}else if(method == "lancaster")
	{
		//run lancaster
		lancaster(id, pmat, dmat,tt, SS, nFeat, nP, output);
	}else if(method == "ordmeta")
	{
		ordmeta(id, pmat, dmat,tt, nFeat, nP, output);
	}else if(method == "artp")
	{
		//run artp
		string aa = string(perm);
		int nperm = stoi(aa);
		ARTP(id, pmat, dmat,tt, nFeat, nP, nperm, output);
	}

	for(i = 0; i < nFeat; i++)
		free(pmat[i]);
	free(pmat);

	for(i = 0; i < nFeat; i++)
		free(dmat[i]);
	free(dmat);

	return 0;
}

