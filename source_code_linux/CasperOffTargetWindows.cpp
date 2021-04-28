// CasperOffTargetWindows.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
//  main.cpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 3/8/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

#include "pch.h"
#include <iostream>
#include <thread>
#include <string>
#include <vector>
#include "csprRef.h"
#include "RefSequences.h"
#include "gRNA.h"
#include <iomanip>
#include <chrono>
using namespace std;
/* NOTE:  Connect this to a python overlay in order to run a subprocess.  That way you can call on the CASPER_Seq_Finder for any off-target NAG sequences if needed.
 * The NAG sequences don't have to be stored, just temporarily generated and scanned.  If there are multiple sequences to be checked then it can be stored.
 * Algorithm for NAGs: store them in files for common organisms and endonuclease combinations.  Otherwise just scan the genome for potential off targets irrespective of PAM.
 */

vector<string> split(string strToSplit, char delimeter)
{
	stringstream ss(strToSplit);
	string item;
	vector<string> splittedStrings;
	while (getline(ss, item, delimeter))
	{
		splittedStrings.push_back(item);
	}
	return splittedStrings;
}

//int argc, const char * argv[]
int main(int argc, const char * argv[]) {
	
	//argv: query file, is_compressed (T/F), cspr file, db file, output file, CASPERinfo, num_mismatches, threshold, detailed ouptut (T/F), avg output (T/F)
	//Convert all input into string objects:
	string query_file = string(argv[1]);
	string cspr_reference = string(argv[2]);
	string sql_file = string(argv[3]);
	string output_file = string(argv[4]);
	string settings_file = string(argv[5]);
	//Convert settings inputs into appropriate values:
	int NUM_MISMATCHES = stoi(string(argv[6]));
	double THRESHOLD = stod(string(argv[7]));
	bool detailed = false;
	bool average = false;
	if (string(argv[8]).find('T') != string::npos) {
		detailed = true;
	}
	if (string(argv[9]).find('T') != string::npos) {
		average = true;
	}

	//get seed_length and seq_length from CASPERinfo file using endo on .cspr file
	ifstream file;
	string endo_name = cspr_reference.substr(cspr_reference.find_last_of('\\') + 1, cspr_reference.size());
	endo_name = endo_name.substr(0, endo_name.find('.'));
	endo_name = endo_name.substr(endo_name.find('_') + 1, endo_name.size());
	file.open(settings_file);
	string str = "";
	int i = 0;
	int seed_length = 0;
	int sequence_length = 0;
	while (getline(file, str))
	{
		if (str.find(endo_name + ';') != string::npos)
		{
			vector<string> data = split(str, ';');
			seed_length = stoi(string(data[2]));
			sequence_length = stoi(string(data[3]));
			break;
		}
	}
	file.close();
	cout << "Running off-target analysis for: " << cspr_reference << "\n";

	// Object for running actual off-target algorithm:
	
	OnTargets otr;
	otr.se_l = seed_length;
	otr.seq_l = sequence_length;

	
	otr.LoadTargetQuery(query_file);
	otr.loadData(cspr_reference, sql_file);


	otr.ScoreSettings(settings_file, output_file, NUM_MISMATCHES, THRESHOLD, detailed, average, cspr_reference);

	
	otr.run_off_algorithm(16); // input should be the number of threads you want to generate and the output file
	
	cout << "All scores generated. Output file located at: " << output_file << endl;
	

	return 0;
}


