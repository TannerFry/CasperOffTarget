//
//  gRNA.hpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 3/19/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

#ifndef gRNA_h
#define gRNA_h

#include <stdio.h>
#include <vector>
#include <set>
#include <string>
#include <thread>
#include "CSeqTranslate.h"

class gRNA {
public:
	int seed_l = 0;
	int seq_l = 0;
	// setting information:
	void set_location(string s) { location = s; }
	void set_score(string s) { score = stod(s); }
	void set_seq(string s) { seq = s; }
	// retrieving information:
	string get_Clocation();
	int get_score() { return double(score); } // Convert to double so that the return value can be a fraction.
	//string get_Csequence() { return seq; }
	string get_sequence() { return seq; }
	long get_location() { return stol(location); }
	// setting off target information
	void addMatch(long id) { putative_off_seq_ids.insert(id); }
	std::set<long> returnPutativeMatches() { return putative_off_seq_ids; }

	//Tells you whether there were any putative hits found:
	bool hasHits();
	int numHits() { return putative_off_seq_ids.size(); }

	//long get10seq() { return S.decompress_location(compressed_seq); }

private:
	//Class containing decompression algorithms
	SeqTranslate S;

	struct offtargetwithscore {
		int chromscaff;
		long position;
		string sequence;
		double off_score;
	};

	vector<offtargetwithscore> offhits;  // Stores all the offtarget scores for the putative off_sequences

private:
	string location;
	int score;
	string seq;

	

public:
	// Below is the container that keeps track of all the matches that the algorithm finds for this target
	std::set<long> putative_off_seq_ids;
	// Adding to off target vector offscores:
	void addOffScore(double, int, long, string);
	string offtargetscores();

};

#endif /* gRNA_hpp */
