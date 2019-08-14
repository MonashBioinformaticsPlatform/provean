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
// Name        : ScoreMatrix.cpp
// Author      : Yongwook Choi
//============================================================================

#include <iostream>
#include "ScoreMatrix.h"

using namespace std;

ScoreMatrix::ScoreMatrix() {
	// TODO Auto-generated constructor stub

}

ScoreMatrix::~ScoreMatrix() {
	// TODO Auto-generated destructor stub
}

void ScoreMatrix::SetScoreMatrix(string matrix_name, int gap_open, int gap_extend)
{
	matrix_name_ = matrix_name;
	gap_open_ = gap_open;
	gap_extend_ = gap_extend;

	if (matrix_name_.compare("BLOSUM62") == 0) {
		for (int i=0; i<MAX_AA; i++) {
			for (int j=0; j<MAX_AA; j++) {
				matrix_[i][j] = ::matrix_blosum_62[i][j];
			}
		}
	} else if (matrix_name_.compare("BLOSUM80") == 0) {
		for (int i=0; i<MAX_AA; i++) {
			for (int j=0; j<MAX_AA; j++) {
				matrix_[i][j] = ::matrix_blosum_80[i][j];
			}
		}
	} else {	// use BLOSUM62 as default
		cerr << "Unidentified substitution matrix name: " << matrix_name << endl;
		cerr << "BLOSUM62 will be used" << endl;
		for (int i=0; i<MAX_AA; i++) {
			for (int j=0; j<MAX_AA; j++) {
				matrix_[i][j] = ::matrix_blosum_62[i][j];
			}
		}
	}
}
