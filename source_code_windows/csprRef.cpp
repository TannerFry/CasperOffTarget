//
//  csprRead.cpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 3/8/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

#include "pch.h"
#include "csprRef.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <assert.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "sqlite3.h"

using namespace std;

string seeds = "";
vector<vector<string>> multi_locs_global;
int sequence_length = 0;

vector<string> split(string s, const string delimiter)
{
	size_t pos = 0;
	string token;
	vector <string> tails;
	while ((pos = s.find(delimiter)) != string::npos) 
	{
		token = s.substr(0, pos);
		tails.push_back(token);
		s.erase(0, pos + delimiter.length());
	}
	tails.push_back(s);
	return tails;
}

int callback(void *NotUsed, int argc, char **argv, char **azColName)
{
	seeds += string(argv[0]);
	vector<string> chroms = split(string(argv[1]), ",");
	vector<string> locs = split(string(argv[2]), ",");
	vector<string> three = split(string(argv[3]), ",");
	vector<string> five = split(string(argv[4]), ",");
	vector<string> pam = split(string(argv[5]), ",");
	vector<string> scores = split(string(argv[6]), ",");
	vector<string> multi;
	
	//fill three and five
	for (int i = 0; i < chroms.size(); i++)
	{
		if (three.size() != chroms.size())
		{
			three.push_back("");
		}
		if (five.size() != chroms.size())
		{
			five.push_back("");
		}
	}
	for (int i = 0; i < chroms.size(); i++)
	{
		multi.push_back(chroms[i] + "," + locs[i] + "," + five[i] + "," + three[i] + "," + pam[i] + "," + scores[i]);
	}
	multi_locs_global.push_back(multi);
	return 0;
}

/* This function loads the reference targets. This is the import function for the .cspr file that is used as a reference
 * Function: This file takes in all the lines of the files and first stores it in a basic string vector.  Then the function
 * parses this string into appropriate containers for easier access in the OnTargets class. */
void csprRef::LoadcsprFile(string cspr, string dbfile, int seq_l) 
{
	long chromCounter = 0;
	bool first = true;
	sequence_length = seq_l;
	ifstream file(cspr, ios_base::in | ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
	inbuf.push(boost::iostreams::gzip_decompressor());
	inbuf.push(file);

	string line;
	//Convert streambuf to istream
	istream instream(&inbuf);
	
	while (true) 
	{
		string line;
		if (first == true)
		{
			for (int i = 0; i < 3; i++)
			{
				getline(instream, line);
			}
			first = false;
		}

		getline(instream, line);  //No line should be longer than 100 characters
		//Determine if we have moved to the repeats section of the file:
		if (line.empty()) 
		{
			// consolidate the id of the last unique by getting the length of the Locs vector:
			multistart = Locs.size();
			cout << "Reached Repeats section. Moving on to process repeat structure." << endl;
			Chrpos.push_back(chromCounter); // end of Chrpos tells you the start of the Repeats section.
			processMultis(dbfile);
			break;
			//Determine whether the new line sets the start of a new Chromosome/Scaffold:
		}
		else if (line.find('>') != string::npos)
		{
			Chrpos.push_back(chromCounter);
			cout << "Sequences from beginning to end of this Chromosome/Scaffold: " << chromCounter + 1 << endl;
		//Passes all tests that it is truly a unique sequence line
		}
		else 
		{
			vector<string> target = Msplit(line, ',');  //Split up the line
			Locs.push_back(stol(target[0]));  //Add compressed location to Locs vector
			Scores.push_back(stoi(target[3]));  //Add decompressed score to Scores vector
			RefTargets += target[1];  //Add the sequence to the long sequence string
			chromCounter++; //keeps track of how many targets are in the chromosome/scaffold
		}
	}
	
	file.close();
}

/* Code for processing the Repeats section which requires different processing of the locations and sequences */

void csprRef::processMultis(string dbfile) 
{
	sqlite3 *db;
	char *zErrMsg = 0;
	int rc;
	
	rc = sqlite3_open(dbfile.c_str(), &db);
	string sql = "SELECT seed, chromosome, location, three, five, pam, score FROM 'REPEATS';";
	rc = sqlite3_exec(db, sql.c_str(), callback, 0, &zErrMsg);

	RefTargets += seeds;
	RefTargets = RefTargets.c_str();
	multiLocs = multi_locs_global;

	sqlite3_close(db);
}


/* getChromScaf: Returns which chromosome/scaffold based on the index of the location */
int csprRef::getChrScaf(long id) {
	int i = 0;
	while (id > Chrpos[i]) {
		i++;
	}
	return i;
}


/* A method for splitting the repeats section. */
vector<string> csprRef::Msplit(const string &text, char sep) {
	vector<string> tokens;
	size_t start = 0, end = 0;
	while ((end = text.find(sep, start)) != string::npos) {
		tokens.push_back(text.substr(start, end - start));
		start = end + 1;
	}
	tokens.push_back(text.substr(start));
	return tokens;
}
