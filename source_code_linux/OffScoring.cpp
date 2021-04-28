//
//  OffScoring.cpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 5/29/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

#include "pch.h"
#include "OffScoring.h"
#include <chrono>
#include <iomanip>
using namespace std;

/* Essentially the initializer.  Called right after the object is initialized in OnTargets object. */
void OffScoring::loadCspr(csprRef *r, string cspr_file, string casper_file, int seq_l, int seed_l) {
	ref = r;
	sequence_length = seq_l;
	seed_length = seed_l;
	fillHsumatrix(cspr_file, casper_file);
}


/* This is the MAIN code for the scoring algorithm */
void OffScoring::score(gRNA* myseq) {
	//Need to get the full uncompressed sequence for each of the putative off sequences
	double avgscore = 1 - scorePutatives(myseq->returnPutativeMatches(), myseq);
	std::string lineoutput = myseq->get_sequence() + ":";
	if (IS_AVERAGE)
	{
		lineoutput += std::to_string(avgscore);
		output.write(lineoutput + "\n");
	}
	if (IS_DETAILED) 
	{
		/* ORIGINAL LINES
			lineoutput += myseq->offtargetscores();
			output.write(lineoutput + "\n")
		*/
		
		output.write(lineoutput + std::to_string(avgscore) + "\n");
		output.write(myseq->offtargetscores());
		//output.write("\n");
	}
}


/* The following functions define the scoring algorithm. */

/* scorePutatives takes the ids of all putative off target sequences, and loads all of the uncompressed information into
 * an offtarget struct.  It then passes this struct to the individual scoring algorithm.*/
double OffScoring::scorePutatives(set<long> offs, gRNA* onseq) 
{
	
	std::vector<offtarget> decomposed_offs;
	for (auto it = offs.begin(); it != offs.end(); ++it) 
	{
		offtarget myoff;
		long id = *it;
		// CHECK TO SEE IF MATCH IS IN THE REPEATS SECTION, THIS REQUIRES MULTIPLE FINAL SEQUENCES TO BE RETURNED:
		if (id > ref->multiStart()) 
		{
			
			string c_base_seq = ref->AccessRefString()->substr(id * sequence_length, seed_length);
			// get the information from the multilocs vector
			long multirelloc = id - ref->multiStart();  //finds the relative location of the sequence for getting multi information
			vector<string> curmultis = ref->getMultis(multirelloc);
			// Loop to go through all the permutations in the curmultis vector to then put them into myoff and then into decomposed offs:
			for (int j = 0; j < curmultis.size(); j++) 
			{
				std::vector<std::string> ind_multi = S.Msplit(curmultis[j], ',');
				myoff.chromscaff = std::stoi(ind_multi[0]);
				myoff.position = stol(ind_multi[1]);
				myoff.on_score = stod(ind_multi[3]);

				// concatenate the base sequence and the tail and add it to the sequence:
				std::string the_tailseq;
				if (ind_multi[2].find('+') != string::npos) 
				{
					the_tailseq = ind_multi[2].substr(0, ind_multi[2].find('+'));
				}
				else 
				{
					the_tailseq = ind_multi[2].substr(0, ind_multi[2].find('-'));
				}

				
				myoff.sequence = the_tailseq + c_base_seq;
				
				// Check to make sure that the putative off target is not a self-match:
				if (myoff.sequence != onseq->get_sequence()) 
				{
					// Put the offtarget object into the decomposed offs vector:
					decomposed_offs.push_back(myoff);
				}
			}
		}
		else 
		{
			//myoff.sequence = S.decompress(ref->AccessRefString()->substr(id, 0), sequence_length);  // find out the actual location of the id and get the string there
			myoff.sequence = ref->AccessRefString()->substr(id * sequence_length, sequence_length);
			myoff.chromscaff = ref->getChrScaf(id);
			myoff.position = ref->getLoc(id);
			myoff.on_score = ref->getScore(id);
			
			// Check to make sure that the putative off target is not a self-match:
			if (myoff.sequence != onseq->get_sequence()) {
				// Put the offtarget object into the decomposed offs vector:
				decomposed_offs.push_back(myoff);
			}
		}
	}
	// These variables are to keep track of the scores counted so that an average can be reported.
	double score_tot = 1;
	int score_num = 1;

	// Now go through all of the decomposed offs, score them and add them to the onseq offscore vector object
	for (int i = 0; i < decomposed_offs.size(); i++) 
	{
		double singleScore = scoreStruct(decomposed_offs[i], onseq);
		if (singleScore < 1.0) 
		{
			if (singleScore > THRESHOLD) 
			{
				onseq->addOffScore(singleScore, decomposed_offs[i].chromscaff, decomposed_offs[i].position, decomposed_offs[i].sequence);
				
			}
			score_tot += 1 - singleScore;
			score_num++;
		}
	}
	return score_tot / score_num;
}

/* The most basic scoring function.  Identifies mismatches of the decompressed strings.
 * Later this will also be the place where locational effects are taken into account. */
double OffScoring::scoreStruct(offtarget oid, gRNA* on) {
	on->seed_l = seed_length;
	on->seq_l = sequence_length;
	// Get the R ratio out of the way:
	double Rratio = on->get_score() / oid.on_score;
	// Score the mismatches of the two sequences:
	std::vector<int> mismatches;
	std::vector<std::string> mismatch_id;
	std::string on_seq = on->get_sequence();
	// Compare the sequences of the on and offtarget.  Go in reverse order to match seed. Need to include bulges in 2.0
	for (int i = oid.sequence.size() - 1; i >= 0; i--) {
		if (oid.sequence[i] != on_seq[i]) {
			mismatches.push_back(i);
			std::string mstr = std::string() + on_seq[i] + S.revcom(oid.sequence[i]);
			mismatch_id.push_back(mstr);
		}
		if (mismatches.size() > MISMATCHES) {
			return 1.0;
		}
	}
	/* From here the mismatch vector is sent to the subscoring functions and the score is tallied with the appropriate algorithmic combination. */
	return ((sqrt(Sh_score(mismatches, mismatch_id)) + St_score(mismatches)) * pow(Ss_score(mismatches), 6) * pow(Rratio, 2)) / 4;


}

/* SUBSCORE SPECIFIC FUNCTIONS */

/* Matrix function, currently set up by downloading the matrix from the CASPERinfo file. */
/* TEMPORARY: Until we get the new matrix to work, just hardcode in the Hsu matrix. */
double OffScoring::Sh_score(vector<int> mismatches, vector<string> mid) {
	double tot_sh = 1.0;
	for (int i = 0; i < mismatches.size(); i++) 
	{	
		tot_sh *= Hsu_Matrix[mid[i]][mismatches[i]];
	}
	return tot_sh;
}

/* Step function, currently set for spCas9 here */
double OffScoring::Ss_score(vector<int> mismatches) {
	double tot_ss = 1.0;
	for (int i = 0;i < mismatches.size(); i++) 
	{
		if (mismatches[i] < 6) 
		{
			tot_ss -= 0.1;
		}
		else if (mismatches[i] < 12) 
		{
			tot_ss -= 0.05;
		}
		else 
		{
			tot_ss -= 0.0125;
		}
	}
	return tot_ss;
}

/* Generic inverse distance score */
double OffScoring::St_score(vector<int> mismatches) {
	double tot_st = 3.5977;
	for (int i = 0; i < mismatches.size(); i++) 
	{
		tot_st -= 1.0 / (mismatches[i] + 1);
	}
	return tot_st / 3.5977;
}


