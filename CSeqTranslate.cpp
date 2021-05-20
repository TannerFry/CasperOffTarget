//
//  CSeqTranslate.cpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 8/31/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

#include "pch.h"
#include "CSeqTranslate.h"

std::string SeqTranslate::compress(std::string input) {
	std::string output;
	long base10compression = 0;
	for (int i = 0; i < input.size(); i++) {
		base10compression += nucleotides.find(input[i])*pow(4, i);
	}
	output = compress(base10compression);
	return output;
}

std::string SeqTranslate::compress(long b10) {
	std::string output;
	while (b10 >= 64) {
		int rem = b10 % 64;
		b10 = b10 / 64;
		output = compress64string[rem] + output;  // Need to structure this way and not += because of the reverse order of compression used.
	}
	output = compress64string[b10] + output;
	return output;
}


std::string SeqTranslate::decompress(std::string input, int seq_size) {
	std::string output;
	long long b10input = decompress_location(input);
	while (b10input >= 4) {
		int rem = b10input % 4;
		b10input = b10input / 4;
		output += nucleotides[rem];
	}
	output += nucleotides[b10input];
	// The for loop fixes the problem of A's not being added to the end because they are removed on compression
	for (int i = output.size(); i < seq_size; i++) {
		output += 'A';
	}
	return output;
}


long long SeqTranslate::decompress_location(std::string input) {
	long long output = 0;
	for (int i = 0; i < input.size(); i++) {
		int power = input.size() - (i + 1);
		int index = compress64string.find(input[i]);
		if (index != std::string::npos) {
			//std::cout << index*pow(64,power) << std::endl;
			output += index * pow(64, power);
		}
	}
	//std::cout << output << std::endl;
	return output;
}



/* string output for revcom.  This takes in a string and outputs the reverse complement. */
std::string SeqTranslate::revcom(std::string in) {
	std::string out;
	for (int i = in.size() - 1;i >= 0;i--) {
		out += revcom(char(in[i]));
	}
	return out;
}

/*character output for revcom.  Also called by the string revcom function. */
char SeqTranslate::revcom(char c) {
	char n = 'N';
	switch (c) {
	case 'A':
		n = 'T';
		break;
	case 'T':
		n = 'A';
		break;
	case 'G':
		n = 'C';
		break;
	case 'C':
		n = 'G';
		break;
	}
	return n;
}

/* A method for splitting the repeats section. */
std::vector<std::string> SeqTranslate::Msplit(const std::string &text, char sep) {
	std::vector<std::string> tokens;
	std::size_t start = 0, end = 0;
	while ((end = text.find(sep, start)) != std::string::npos) {
		tokens.push_back(text.substr(start, end - start));
		start = end + 1;
	}
	tokens.push_back(text.substr(start));
	return tokens;
}
