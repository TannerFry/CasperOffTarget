//
//  OffScoring.hpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 5/29/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

using namespace std;

#ifndef OffScoring_h
#define OffScoring_h

#include <stdio.h>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <math.h>
#include "gRNA.h"
#include "csprRef.h"
#include "CSeqTranslate.h"
#include "FileOp.h"
#include <algorithm>
#include <sstream>

class OffScoring {
public:
	int sequence_length;
	int seed_length;
	void loadCspr(csprRef*, string, string, int, int);
	void setOutputFile(std::string fn) 
	{
		output.openWrite(fn); 
		if(IS_AVERAGE)
		{ 
			output.write("AVG OUTPUT\n");
		}
		else
		{
			output.write("DETAILED OUTPUT\n");
		}
	}
	void settings(int m, double t, bool d, bool a) {
		MISMATCHES = m;
		THRESHOLD = t;
		IS_DETAILED = d;
		IS_AVERAGE = a;

	}

	//Key algorithm is called below (a bunch of loadPutative calls):
	void score(gRNA*);

	/* Parser and Temporary storage containers for storing the uncompressed putative off sequences */
private:
	double scorePutatives(std::set<long>, gRNA*); // single target parser (all putatives)
	struct offtarget {
		int chromscaff;
		long position;
		std::string sequence;
		double on_score;
	};

	struct missid {
		char on;
		char off;
	};

	double scoreStruct(offtarget, gRNA*);

private:
	double THRESHOLD;
	int MISMATCHES;
	bool IS_DETAILED;
	bool IS_AVERAGE;

	std::vector<std::vector<double>> Off_Matrix;
	/* Below are two hardcoded scoring functions.  In CASPER 2.0, make these not hardcoded but from the off settings. */
	double Ss_score(std::vector<int>);
	double St_score(std::vector<int>);
	double Sh_score(std::vector<int>, std::vector<std::string>);

	SeqTranslate S;
	WriteFile output;
	csprRef* ref;

	/* Below is the hardcoded matrix elements from the Hsu matrix.  In CASPER 2.0 we would like to change these to more generalizeable features. */
private:
	std::map<std::string, std::vector<double>> Hsu_Matrix;
	void fillHsumatrix(string cspr_file, string casper_file) 
	{
		ifstream file;
		string endo_name = cspr_file.substr(cspr_file.find_last_of('\\') + 1, cspr_file.size());
		endo_name = endo_name.substr(0, endo_name.find('.'));
		endo_name = endo_name.substr(endo_name.find('_')+1, endo_name.size());
		file.open(casper_file);
		string str = "";
		int i = 0;
		while (getline(file, str)) 
		{
			if (str.find("MATRIX-" + endo_name) != string::npos)
			{
				str.erase(remove(str.begin(), str.end(), '\n'), str.end());
				getline(file,str);
				while (str != "")
				{
					float n;
					stringstream stream(str);
					while (stream >> n) 
					{
						if (i == 0) 
						{
							Hsu_Matrix["GT"].push_back(n);
						}
						else if(i == 1)
						{
							Hsu_Matrix["AC"].push_back(n);
						}
						else if (i == 2)
						{
							Hsu_Matrix["GG"].push_back(n);
						}
						else if (i == 3)
						{
							Hsu_Matrix["TG"].push_back(n);
						}
						else if (i == 4)
						{
							Hsu_Matrix["TT"].push_back(n);
						}
						else if (i == 5)
						{
							Hsu_Matrix["CA"].push_back(n);
						}
						else if (i == 6)
						{
							Hsu_Matrix["CT"].push_back(n);
						}
						else if (i == 7)
						{
							Hsu_Matrix["GA"].push_back(n);
						}
						else if (i == 8)
						{
							Hsu_Matrix["AA"].push_back(n);
						}
						else if (i == 9)
						{
							Hsu_Matrix["AG"].push_back(n);
						}
						else if (i == 10)
						{
							Hsu_Matrix["TC"].push_back(n);
						}
						else if (i == 11)
						{
							Hsu_Matrix["CC"].push_back(n);
						}
					}
					
					i += 1;
					getline(file, str);
				
				}
			}
		}
		
		file.close();

		
		/*
		Hsu_Matrix["GT"] = { 1, 1.6533,0.9030,1.5977,0.9235,0.8070,0.9632,1.0163,0.2658,0.7119,1.5211,0.6942,1.0434,0.5255,0.8981,0.7164,0.8399,0.5376,0.2821,0.6898 };
		Hsu_Matrix["AC"] = { 1, 1.5142,1.1597,1.6582,0.9924,0.0247,0.5522,1.8687,0.7737,0.9270,0.7292,0.4842,0.4824,0.7060,1.0221,0.0181,0.3496,0.1811,0.1362,0.2700 };
		Hsu_Matrix["GG"] = { 1, 1.3234,1.4157,1.2967,1.2060,0.9524,0.2304,1.0163,0.8100,1.1559,0.7075,1.5791,0.3490,0.0899,0.0497,0.0045,0.2267,0.2153,0.5250,0.4965 };
		Hsu_Matrix["TG"] = { 1, 1.5366,1.2771,1.2689,1.2197,1.0645,0.7791,1.2445,0.9885,1.5319,0.7184,1.7381,0.4166,0.1285,0.0720,0.0549,0.2261,0.3119,0.1343,0.0601 };
		Hsu_Matrix["TT"] = { 1, 1.7347,0.8215,1.1579,1.1816,0.7380,0.9004,0.8368,0.2997,0.6210,1.1400,0.3561,0.6192,0.1799,0.2665,0.2793,0.2613,0.1152,0.1680,0.3372 };
		Hsu_Matrix["CA"] = { 1, 1.0186,0.9649,2.2504,1.6222,0.2405,0.7561,1.0651,0.1102,1.4293,0.3533,0.6178,0.7269,0.1165,0.0367,0.5013,0.4147,0.1786,0.5315,0.1664 };
		Hsu_Matrix["CT"] = { 1, 1.6719,0.9688,1.0732,1.0869,0.6475,1.0142,0.8635,0.3059,0.4487,0.9046,0.4327,0.5576,0.1379,0.0722,0.3279,0.2420,0.0433,0.1351,0.4403 };
		Hsu_Matrix["GA"] = { 1, 1.1662,0.4544,2.7867,1.0461,0.6036,0.8132,0.7875,0.6882,1.3655,0.1240,0.1953,0.2497,0.0132,0.0227,0.0478,0.3682,0.3175,0.5621,0.4588 };
		Hsu_Matrix["AA"] = { 1, 1.1916,1.0954,2.8709,1.1310,0.5160,0.6439,1.0322,0.5356,1.2868,0.0780,0.2592,0.2675,0.0469,0.0252,0.0052,0.0218,0.1718,0.6970,0.2720 };
		Hsu_Matrix["AG"] = { 1, 1.4786,1.0820,1.2952,0.7450,0.9763,0.4912,0.9272,0.6022,1.0375,0.3047,0.8210,0.0303,0.0365,0.0592,0.0253,0.1553,0.1006,0.2175,0.0275 };
		Hsu_Matrix["TC"] = { 1, 1.0400,0.9954,1.6466,1.3410,0.0102,0.5428,2.3401,0.4367,0.2143,0.3405,0.2640,0.0935,0.0462,0.0688,0.0165,0.3659,0.0546,0.0857,0.2543 };
		Hsu_Matrix["CC"] = { 1, 1.0345,1.0478,1.0507,1.4075,0.0540,0.6396,2.0810,0.4585,0.1555,0.1369,0.1026,0.0417,0.0105,0.0458,0.0099,0.2114,0.0552,0.0253,0.0596 };
		*/
	}
};

#endif /* OffScoring_hpp */

