//
//  gRNA.cpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 3/19/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

#include "pch.h"
#include "gRNA.h"

/*
std::string gRNA::get_Clocation() {
	// need to uncompress location into long object
	return compressed_location;
}
*/

//returns whether or not any putative hits were found somewhere on the genome.
bool gRNA::hasHits() {
	if (putative_off_seq_ids.size() > 0) {
		return true;
	}
	else { return false; }
}

void gRNA::addOffScore(double d, int c, long l, std::string s) {
	offtargetwithscore myoff;
	myoff.chromscaff = c;
	myoff.position = l;
	myoff.off_score = d;
	myoff.sequence = s;
	offhits.push_back(myoff);
}

std::string gRNA::offtargetscores() {
	std::string return_string;
	for (int i = 0;i < offhits.size();i++) {
		return_string += std::to_string(offhits[i].off_score) + "," + std::to_string(offhits[i].chromscaff) + "," + std::to_string(offhits[i].position) + "," + offhits[i].sequence + "\n";
	}
	return return_string;
}
