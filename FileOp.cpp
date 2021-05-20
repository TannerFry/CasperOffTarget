//
//  FileOp.cpp
//  CasperOffTarget
//
//  Created by Brian Mendoza on 6/4/18.
//  Copyright © 2018 University of Tennessee. All rights reserved.
//

#include "pch.h"
#include "FileOp.h"

FileOp::FileOp() {};

void FileOp::open(std::string fname)
{
	filename = fname;
	FileStream = new std::ifstream;
	FileStream->open(fname.c_str(), std::ifstream::in);
}

bool FileOp::newLine()
{
	if (!FileStream->eof()) return true;
	return false;
}

bool FileOp::endFile()
{
	if (FileStream) return true;
	return false;
}

/* Grabs a line from the file.  You can specify how long the character stream with the passed through
 integer. */
std::string FileOp::getLine(int l)
{
	char nts[5000000];
	FileStream->getline(nts,l);
	return std::string(nts);
}

/* Closes the file stream */
void FileOp::closeFile()
{
	FileStream->close();
}

/* This function is used if you want to get a line separated by a given token */
std::vector<std::string> FileOp::getSepLine(int count, char sep) {
	return Msplit(getLine(count), sep);
}

/* A method for splitting the repeats section. */
std::vector<std::string> FileOp::Msplit(const std::string &text, char sep) {
	std::vector<std::string> tokens;
	std::size_t start = 0, end = 0;
	while ((end = text.find(sep, start)) != std::string::npos) {
		tokens.push_back(text.substr(start, end - start));
		start = end + 1;
	}
	tokens.push_back(text.substr(start));
	return tokens;
}

/* BELOW ARE ALL THE FUNCTIONS FOR WRITING TO A FILE WriteFile Class. */
void WriteFile::openWrite(std::string fn) {
	outputfile.open(fn);
}

void WriteFile::write(std::string stuff) {
	outputfile << stuff;
}

void WriteFile::write(long stuff) {
	outputfile << stuff;
}

void WriteFile::write(double stuff) {
	outputfile << stuff;
}

void WriteFile::writeLine(std::string line) {
	outputfile << line << "\n";
}

void WriteFile::close() {
	outputfile.close();
}

/* Process settings.  This is a specific function used for generating an OffScoring class and filling out the parameters needs to have forward declaration.*/
//OffScoring FileOp::processSettings() {
	//do nothing right now
	/*// Get to the line with the HSU_MATRIX
	bool atmatrix = false;
	while (!atmatrix) {
		// check if it is the matrix title line
		if (getLine(20).find("MATRIX")) {
			atmatrix = true;
		}
	}
	// Start processing the items in the matrix
	std::string matrixline = getLine(1000);
	while (matrixline.find("---") != 0) {
		std::vector<std::string> Trow = Msplit(matrixline,'\t');
		std::vector<double> row;
		// Go through the Trow and create a vector with doubles instead of strings
		for (int i=0;i<Trow.size();i++) {
			row.push_back(std::stod(Trow[i]));
		}
		matrixline = getLine(1000);
	}*/
	//}
