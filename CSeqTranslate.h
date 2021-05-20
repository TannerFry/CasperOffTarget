//
//  CSeqTranslate.hpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 8/31/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

#ifndef CSeqTranslate_h
#define CSeqTranslate_h

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>

using namespace std;

class SeqTranslate {
public:
	std::string compress(std::string);
	std::string compress(long);
	std::string decompress(std::string, int);
	long long decompress_location(std::string);
	//two functions to get revcoms of strings or chars:
	std::string revcom(std::string);
	char revcom(char);
	// For simplicity the split function has been added to this class so it can be shared with all the other classes as needed:
	std::vector<std::string> Msplit(const std::string &text, char sep);

private:
	std::string compress64string = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789=/";
	/* To get the nucleotide use nucleotide.at, to get the base4 conversion use nucleotide.find */
	std::string nucleotides = "ATCG";
};


#endif /* CSeqTranslate_hpp */
