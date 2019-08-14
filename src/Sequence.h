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
// Name        : Sequence.h
// Author      : Yongwook Choi
//============================================================================

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>
#include "Common.h"

using namespace std;

class Sequence {
public:
	Sequence();
	virtual ~Sequence();

	bool SetSequence();
	bool SetSequence(char* id, string blast_db_file, string blastdbcmd);
	bool SetSequenceFromFastaFile(const char* id, string filename);
	void Print(FILE* fp);

	string id_;
	string seq_;
	string def_;
	int cluster_id_;
	double e_value_;
	double bit_score_;
	double weight_;

	static char buf_[BUF_SIZE_LARGE];
};

#endif /* SEQUENCE_H_ */
