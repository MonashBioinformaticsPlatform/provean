//============================================================================
// Copyright 2012-2014 J. Craig Venter Institute
// This file is part of PROVEAN.  PROVEAN is free software: you may
// redistribute it and/or modify it under the terms of the GNU General Public
// License version 3.  PROVEAN is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU
// General Public License along with this program. If not, see
// http://www.gnu.org/licenses/gpl.txt.
//============================================================================
// Name        : Sequence.cpp
// Author      : Yongwook Choi
//============================================================================

#include <iostream>
#include "Common.h"
#include "Sequence.h"

using namespace std;

char Sequence::buf_[BUF_SIZE_LARGE];

Sequence::Sequence() {
	cluster_id_ = -1;
}


Sequence::~Sequence() {
	// TODO Auto-generated destructor stub
}

bool Sequence::SetSequence(char* id, string blast_db_file, string blastdbcmd)
{
	FILE* fpipe;
	char command[BUF_SIZE_MED];

	sprintf(command, "%s -db %s -dbtype 'prot' -entry '%s' -line_length 1000", blastdbcmd.c_str(), blast_db_file.c_str(), id);

	if (!(fpipe = popen(command, "r"))) {
		cerr << "error in pipe" << endl;
		exit(1);
	}

	int len;

	memset(buf_, 0x00, sizeof(buf_));
	len = fread(buf_, sizeof(char), sizeof(buf_), fpipe);

	pclose(fpipe);

	id_ = id;

	return SetSequence();
}

bool Sequence::SetSequence()
{
	char* strptr;

	strptr = strtok(buf_, "\r\n\0");
	if (strptr == NULL) {
		return false;
	}

	def_.assign(strptr, 100);

	while (1) {
		strptr = strtok(NULL, " \r\n");
		if (strptr == NULL)
			break;
		seq_.append(strptr);
	}

	return true;
}

bool Sequence::SetSequenceFromFastaFile(const char* id, string filename)
{
	FILE* fp = fopen(filename.c_str(), "r");
	if (fp == NULL) {
		Log("File open error\n", true);
		cerr << filename << endl;
		exit(1);
		return false;
	}

	id_ = id;

	char* strptr;
	char buf[BUF_SIZE_LARGE];
	while (!feof(fp)) {
		if (fgets(buf, BUF_SIZE_LARGE, fp) == NULL)
			break;

		if (buf[0] == '>') {
			// def line
			def_.assign(buf, 100);
			while (strchr(buf, '\n') == NULL) {
				fgets(buf, sizeof(buf), fp);
			}
		} else {
			strptr = strtok(buf, " \r\n");
			if (strptr == NULL)
				break;
			seq_.append(strptr);
		}
	}

	fclose(fp);
	return true;
}

void Sequence::Print(FILE* fp)
{
	int segment_len = 60;
	fprintf(fp, "%s\n", id_.c_str());
	for (int i=0; i<(int)seq_.length(); i+=segment_len) {
		fprintf(fp, "%5d  %s\n", i+1, seq_.substr(i, segment_len).c_str());
	}

	fflush(fp);
}


