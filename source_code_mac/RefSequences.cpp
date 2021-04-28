//
//  RefSequences.cpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 3/16/18.
//  Copyright � 2018 University of Tennessee. All rights reserved.
//

#include "pch.h"
#include "RefSequences.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <functional>
#include <assert.h>
#include "csprRef.h"
#include "OffScoring.h"
#include <math.h>
#include <chrono>
#include <iomanip>

using namespace std;
char* r;

/* CONTAINER AND DATA MANIPULATION FUNCTIONS
 * loadData
 * LoadTargetQuery */


void bitap_bitwise_search(const char *text, const char *pattern, set<long> &seq, int seq_l)
{
	int m = strlen(pattern);
	unsigned long R;
	unsigned long pattern_mask[CHAR_MAX + 1];
	int i;
	/* Initialize the bit array R */
	R = ~1;

	/* Initialize the pattern bitmasks */
	for (int i = 0; i <= CHAR_MAX; ++i)
		pattern_mask[i] = ~0;
	for (int i = 0; i < m; ++i)
		pattern_mask[pattern[i]] &= ~(1UL << i);

	for (int i = 0; text[i] != '\0'; ++i) {
		/* Update the bit array */
		R |= pattern_mask[text[i]];
		R <<= 1;

		if (0 == (R & (1UL << m)))
		{
			seq.insert((i - m + 1) / seq_l);
		}
	}
}


int factors(int num) {
	int i;
	int x;
	for (i = 1; i <= num; i++) 
	{
		if (num % i == 0 && i < 6)
		{
			x = i;
		}
	}
	return x;
}


void OnTargets::loadData(string f_name, string dbfile) {
	//Open the file and load all the targets from the unique and repeats sections
	ref.LoadcsprFile(f_name, dbfile, seq_l);
}


/* This function generates off target scores for all of the sequences that were identifed in the inital searching function: FindSimilars */
void OnTargets::ScoreSettings(string settings_filename, string output_filename, int mismatches, double thres, bool det, bool avg, string cspr_file) {
	//FileOp sfile;
	fileop.open(settings_filename);
	scoreGenerator.settings(mismatches, thres, det, avg);
	scoreGenerator.loadCspr(&ref, cspr_file, settings_filename, seq_l, se_l);
	scoreGenerator.setOutputFile(output_filename);
	fileop.closeFile();
}


/* Grab the ontargets from the OFF_QUERY file */
void OnTargets::LoadTargetQuery(string query_file) {
	ifstream fileop(query_file);
	while (fileop) 
	{
		gRNA *query = new gRNA;
		string gcsline;
		getline(fileop, gcsline);
		if (gcsline.size() <= 1) 
		{
			break;
		}
		
		vector<string> gcs = Msplit(gcsline, ';');
		string sequence;
		query->set_location(gcs[0]);
		query->set_seq(gcs[1]);
		query->set_score(gcs[3]); // index 2 is the PAM, skip this and go to the score
		//cout << gcs[0] << ";" << gcs[1] << ";" << gcs[3] << endl;

		base_seqs.push_back(query);
	}
	fileop.close();
}


/* ALGORITHMIC FUNCTIONS: run_off_algorithm calls instances of findSimilars which searches and calls the scoring function from OffScoring */
void OnTargets::run_off_algorithm(int thr) 
{
	std::cout << "Running Off Target Algorithm for " << base_seqs.size() << " sequences... " << endl;
	vector<gRNA*> base = base_seqs;
	/* Run 16 threads to get through all of the gRNAs in question */
	r = strdup(ref.AccessRefString()->c_str());
	int i = 0;
	int total_size = base.size();
	scoreGenerator.seed_length = se_l;
	scoreGenerator.sequence_length = seq_l;

	//single process
	/*
	for (int j = 0; j < base.size(); j++)
	{
		cout << "Percentage of sequences scored: " << (double(j) / double(base.size())) * 100 << "%" << endl;
		findSimilars(base[j]);
	}
	*/

	//multi process
	while ((base.size() - i) / thr > 0) 
	{
		std::cout << "Percentage of sequences scored: " << (double(i) / double(base.size())) * 100 << "%" << std::endl;
		std::vector<std::thread*> running_threads;

		std::thread t0([this, &base, &i]() 
		{
			findSimilars(base[i]);
		});
		std::thread t1([this, &base, &i]() {
			findSimilars(base[i + 1]);
		});
		std::thread t2([this, &base, &i]() {
			findSimilars(base[i + 2]);
		});
		std::thread t3([this, &base, &i]() {
			findSimilars(base[i + 3]);
		});
		std::thread t4([this, &base, &i]() {
			findSimilars(base[i + 4]);
		});
		std::thread t5([this, &base, &i]() {
			findSimilars(base[i + 5]);
		});
		std::thread t6([this, &base, &i]() {
			findSimilars(base[i + 6]);
		});
		std::thread t7([this, &base, &i]() {
			findSimilars(base[i + 7]);
		});
		std::thread t8([this, &base, &i]() {
			findSimilars(base[i + 8]);
		});
		std::thread t9([this, &base, &i]() {
			findSimilars(base[i + 9]);
		});
		std::thread t10([this, &base, &i]() {
			findSimilars(base[i + 10]);
		});
		std::thread t11([this, &base, &i]() {
			findSimilars(base[i + 11]);
		});
		std::thread t12([this, &base, &i]() {
			findSimilars(base[i + 12]);
		});
		std::thread t13([this, &base, &i]() {
			findSimilars(base[i + 13]);
		});
		std::thread t14([this, &base, &i]() {
			findSimilars(base[i + 14]);
		});
		std::thread t15([this, &base, &i]() {
			findSimilars(base[i + 15]);
		});
		
		running_threads.push_back(&t0);
		running_threads.push_back(&t1);
		running_threads.push_back(&t2);
		running_threads.push_back(&t3);
		running_threads.push_back(&t4);
		running_threads.push_back(&t5);
		running_threads.push_back(&t6);
		running_threads.push_back(&t7);
		running_threads.push_back(&t8);
		running_threads.push_back(&t9);
		running_threads.push_back(&t10);
		running_threads.push_back(&t11);
		running_threads.push_back(&t12);
		running_threads.push_back(&t13);
		running_threads.push_back(&t14);
		running_threads.push_back(&t15);

		// Join all the running threads before continuing into the next iteration of the loop
	
		for (int j = 0;j < running_threads.size();j++) {
			running_threads[j]->join();
		}
		i += thr;
	}

	for (int k = i; k < base.size(); k++) 
	{
		findSimilars(base[k]);
	}
	
	
	std::cout << "Percentage of sequences scored: " << "100.00%" << endl;

	std::cout << "Done searching for putative matches. Beginning scoring process.                 " << endl; //spaces are to delete random numbers
}


// This is the target of the thread function
void OnTargets::findSimilars(gRNA* seq) 
{
	//See if any subset of the sequence appears in the csprRef.reftargets object:
	// break input sequence into 4 subsequences:
	int x = factors(seq_l - 9);
	for (int s = 0; s + 9 < seq_l + 1; s += x) 
	{
		bitap_bitwise_search(r, seq->get_sequence().substr(s, 9).c_str(), seq->putative_off_seq_ids, seq_l);
	}
	
	// check to see if this particular query has any putative off sequence hits:
	if (seq->putative_off_seq_ids.size() > 0) 
	{
		//now run the scoring algorithm by getting all the information
		scoreGenerator.score(seq);
	}
	delete seq;
}


vector<string> OnTargets::Msplit(const string &text, char sep) {
	vector<string> tokens;
	size_t start = 0, end = 0;
	while ((end = text.find(sep, start)) != string::npos) {
		tokens.push_back(text.substr(start, end - start));
		start = end + 1;
	}
	tokens.push_back(text.substr(start));
	return tokens;
}


