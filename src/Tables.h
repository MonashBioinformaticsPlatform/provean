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
// Name        : Tables.h
// Author      : Yongwook Choi
//============================================================================

#ifndef TABLES_H_
#define TABLES_H_

#include <string>
#include "Common.h"
#include "ScoreMatrix.h"

using namespace std;

class Tables {
private:
	bool SetTable(int len_q, int len_s);
	void FillForwardTable_Full();
	void FillForwardTable_Partial(int pos_from, int pos_to);
	int FillBackwardTable_Full();
	int FillBackwardTable_Partial(int pos_from, int pos_to);
	int FillBackwardTable_Stud();
	bool SetTable();
	void ComputeAlignment_Full();
	void ComputeAlignment_Partial(int pos_to);
	void ComputeAlignment_Partial();
	int query_aligned_start_pos_;
	int subject_aligned_start_pos_;

	// for moving tables
	int constraint_status_;
	int pos_query_;
	int pos_subject_;
	int score_;

	int** f_score_;	// forward scores
	int** f_score_q_;	// forward score with a gap at the end of query
	int** f_score_s_;	// forward score with a gap at the end of subject
	int** b_score_;	// backward scores
	int** b_score_q_;	// backward score with a gap at the end of query
	int** b_score_s_;	// backward score with a gap at the end of subject

	int** b_score_stud_;
	int** b_score_stud_q_;
	int** b_score_stud_s_;

	int* max_score_forward_;		// max score among [0][subject_len] .. [i][subject_len] is stored at i-th element
	int* max_score_backward_;

	int len_x_;
	int len_y_;

	string query_;
	string subject_;
	string query_aligned_;
	string subject_aligned_;

	ScoreMatrix* p_mat_;

	bool moving_table_;
	unsigned int stud_size_;
	int margin_;

	int pos_start_;
	int pos_end_;

public:
	Tables();
	virtual ~Tables();

	bool SetTable(string query, string subject, ScoreMatrix* p_mat);
	void FreeTables();

	void FillForwardTable();
	int FillBackwardTable();	// return reference score

	int GetRefAlignmentScore();
	int GetRefAlignmentScore_B();
	int GetVarAlignmentScore(int pos_left, int pos_right, string str_ins);

	void ComputeAlignment();
	void PrintAlignment(int segment_len);
	double GetIdentity();
};

#endif /* TABLES_H_ */
