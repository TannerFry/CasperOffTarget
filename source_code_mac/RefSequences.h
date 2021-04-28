//
//  RefSequences.hpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 3/16/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

#ifndef RefSequences_h
#define RefSequences_h

#include <stdio.h>
#include <vector>
#include <string>
#include <thread>
#include <stdlib.h>
#include <stdio.h>
#include "csprRef.h"
#include "FileOp.h"
#include "gRNA.h"
#include "OffScoring.h"
#include "CSeqTranslate.h"

using namespace std;

class OnTargets {
public:
	//OnTargets();
	int se_l = 0;
	int seq_l = 0;
	void loadData(string, string); //pass filename of .cspr file for reference
	
	void LoadTargetQuery(string);  // pass OFF_QUERY file to the class
	void set_base_seqs(vector<gRNA*> x) { base_seqs = x; }

	// calls the scoring algorithm iterating through putative_off_seqs (passed object is settings file)
	void ScoreSettings(string, string, int, double, bool, bool, string);

	//generic algorithm that loops through all targets calling findSimilars. Initiates threads. Iterates through base_seqs and scores them with the scoring algorithm
	void run_off_algorithm(int);
	
private:
	void findSimilars(gRNA*);  //runs algorithm on individual target
	vector<string> Msplit(const string &text, char sep);
	vector<gRNA*> base_seqs;  //stores all queries loaded by the public function set_base_seqs

private:
	csprRef ref;
	OffScoring scoreGenerator;
	SeqTranslate S;
	FileOp fileop;
};

#endif /* RefSequences_hpp */
