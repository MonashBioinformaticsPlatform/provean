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
// Name        : Tables.cpp
// Author      : Yongwook Choi
//============================================================================

#include "Tables.h"
#include "Common.h"
#include <iostream>
#include <math.h>

#define NINFTY -100000

#define MEM_LIMIT 512000000	// 512M

#define CONSTRAINT_NONE 0
#define CONSTRAINT_SUBJECT_GAP 1
#define CONSTRAINT_QUERY_GAP 2


Tables::Tables() {
	f_score_ = NULL;
	f_score_q_ = NULL;
	f_score_s_ = NULL;
	b_score_ = NULL;
	b_score_q_ = NULL;
	b_score_s_ = NULL;
	b_score_stud_ = NULL;
	b_score_stud_q_ = NULL;
	b_score_stud_s_ = NULL;

	max_score_forward_ = NULL;
	max_score_backward_ = NULL;

	margin_ = 100;
}

Tables::~Tables() {
	FreeTables();
}


bool Tables::SetTable(int len_q, int len_s)
{
	if (len_q * len_s > 80000000) // 80 M
	{
		Log("Too much memory required!!\n", true);
		return false;
	}



	len_x_ = len_q+2;
	len_y_ = len_s+2;

	f_score_ = new int*[len_x_];
	f_score_q_ = new int*[len_x_];
	f_score_s_ = new int*[len_x_];
	b_score_ = new int*[len_x_];
	b_score_q_ = new int*[len_x_];
	b_score_s_ = new int*[len_x_];

	for (int i=0; i<len_x_; i++) {
		f_score_[i] = new int[len_y_];
		f_score_q_[i] = new int[len_y_];
		f_score_s_[i] = new int[len_y_];
		b_score_[i] = new int[len_y_];
		b_score_q_[i] = new int[len_y_];
		b_score_s_[i] = new int[len_y_];
	}

	return true;
}

bool Tables::SetTable()
{
	f_score_ = new int*[len_x_];
	f_score_q_ = new int*[len_x_];
	f_score_s_ = new int*[len_x_];
	b_score_ = new int*[len_x_];
	b_score_q_ = new int*[len_x_];
	b_score_s_ = new int*[len_x_];

	for (int i=0; i<len_x_; i++) {
		f_score_[i] = new int[len_y_];
		f_score_q_[i] = new int[len_y_];
		f_score_s_[i] = new int[len_y_];
		b_score_[i] = new int[len_y_];
		b_score_q_[i] = new int[len_y_];
		b_score_s_[i] = new int[len_y_];
	}

	return true;
}

bool Tables::SetTable(string query, string subject, ScoreMatrix* p_mat)
{
	query_ = query;
	subject_ = subject;
	p_mat_ = p_mat;

	this->FreeTables();

	max_score_forward_ = new int[query_.length()+2];
	max_score_backward_ = new int[query_.length()+2];

	if (query.length()*subject.length()*24 <= MEM_LIMIT) {
		pos_start_ = 0;
		pos_end_ = query_.length()+1;
		moving_table_ = false;
		len_x_ = query_.length()+2;
		len_y_ = subject_.length()+2;
		return this->SetTable();
	} else {
		stud_size_ = max((int)((double)MEM_LIMIT/(24*subject.length())), 500);
		pos_start_ = 0;
		pos_end_ = stud_size_;

		int num_stud = query_.length()/stud_size_;

		b_score_stud_ = new int*[num_stud];
		b_score_stud_q_ = new int*[num_stud];
		b_score_stud_s_ = new int*[num_stud];

		for (int i=0; i<num_stud; i++) {
			b_score_stud_[i] = new int[subject_.length()+2];
			b_score_stud_q_[i] = new int[subject_.length()+2];
			b_score_stud_s_[i] = new int[subject_.length()+2];
		}
		moving_table_ = true;
		len_x_ = stud_size_+margin_;
		len_y_ = subject_.length()+2;
		return this->SetTable();
	}

}

void Tables::FreeTables()
{
	if (f_score_ != NULL) {
		for (int i=0; i<len_x_; i++) {
			if (f_score_[i] != NULL) {
				delete[] f_score_[i];
			}
			if (f_score_q_[i] != NULL) {
				delete[] f_score_q_[i];
			}
			if (f_score_s_[i] != NULL) {
				delete[] f_score_s_[i];
			}
			if (b_score_[i] != NULL) {
				delete[] b_score_[i];
			}
			if (b_score_q_[i] != NULL) {
				delete[] b_score_q_[i];
			}
			if (b_score_s_[i] != NULL) {
				delete[] b_score_s_[i];
			}
		}
		delete[] f_score_;
		delete[] f_score_q_;
		delete[] f_score_s_;
		delete[] b_score_;
		delete[] b_score_q_;
		delete[] b_score_s_;

		f_score_ = NULL;
		f_score_q_ = NULL;
		f_score_s_ = NULL;
		b_score_ = NULL;
		b_score_q_ = NULL;
		b_score_s_ = NULL;
	}

	if (b_score_stud_ != NULL) {
		for (unsigned int i=0; i<query_.length()/stud_size_; i++) {
			delete[] b_score_stud_[i];
			delete[] b_score_stud_s_[i];
			delete[] b_score_stud_q_[i];
		}
		delete[] b_score_stud_;
		delete[] b_score_stud_s_;
		delete[] b_score_stud_q_;
		b_score_stud_ = NULL;
	}

	if (max_score_forward_ != NULL) {
		delete[] max_score_forward_;
		max_score_forward_ = NULL;
	}

	if (max_score_backward_ != NULL) {
		delete[] max_score_backward_;
		max_score_backward_ = NULL;
	}
}


void Tables::FillForwardTable()
{
	if (moving_table_) {
		this->FillForwardTable_Partial(0, stud_size_);
	} else {
		this->FillForwardTable_Full();
	}
}

void Tables::FillForwardTable_Full()
{
#ifdef DEBUG
	fprintf(stderr, "...FillForwardTable_Full\n");
#endif
	f_score_[0][0] = 0;
	f_score_s_[0][0] = NINFTY;
	f_score_q_[0][0] = NINFTY;

	for (unsigned int i=1; i<=subject_.length(); i++) {
		f_score_[0][i] = 0;
		f_score_q_[0][i] = 0;
		f_score_s_[0][i] = NINFTY;
	}

	for (unsigned int i=1; i<=query_.length(); i++) {
		f_score_[i][0] = 0;
		f_score_s_[i][0] = 0;
		f_score_q_[i][0] = NINFTY;
	}

	max_score_forward_[0] = 0;

	int query_aa_idx, subject_aa_idx;

	for (unsigned int i=1; i<=query_.length(); i++) {
		for (unsigned int j=1; j<=subject_.length(); j++) {
			f_score_s_[i][j] = Max3(f_score_[i-1][j]+p_mat_->gap_open_,
					f_score_s_[i-1][j]+p_mat_->gap_extend_,
					f_score_q_[i-1][j]+p_mat_->gap_open_);
			f_score_q_[i][j] = Max3(f_score_[i][j-1]+p_mat_->gap_open_,
					f_score_s_[i][j-1]+p_mat_->gap_open_,
					f_score_q_[i][j-1]+p_mat_->gap_extend_);

			query_aa_idx = aa2idx[(int)query_[i-1]];
			subject_aa_idx = aa2idx[(int)subject_[j-1]];
			if (query_aa_idx == NA || subject_aa_idx == NA) {
				fprintf(stderr, "Unknown AA character, query:[%c], subject:[%c]\n", query_[i-1], subject_[j-1]);
				exit(1);
			}
			f_score_[i][j] = Max3(f_score_q_[i-1][j-1], f_score_s_[i-1][j-1], f_score_[i-1][j-1])
					+ p_mat_->matrix_[query_aa_idx][subject_aa_idx];
		}

		max_score_forward_[i] = max(Max3(f_score_q_[i][subject_.length()], f_score_s_[i][subject_.length()], f_score_[i][subject_.length()]), max_score_forward_[i-1]);
	}

	return;
}



int Tables::FillBackwardTable()
{
	int ref_score;
	if (moving_table_) {
		FillBackwardTable_Stud();
		ref_score = FillBackwardTable_Partial(0, stud_size_);
	} else {
		ref_score = FillBackwardTable_Full();
	}

	return ref_score;
}

int Tables::FillBackwardTable_Full()
{
#ifdef DEBUG
	fprintf(stderr, "...FillBackwardTable_Full\n");
#endif
	// x : query
	// y : subject

	b_score_[query_.length()+1][subject_.length()+1] = 0;
	b_score_s_[query_.length()+1][subject_.length()+1] = NINFTY;
	b_score_q_[query_.length()+1][subject_.length()+1] = NINFTY;

	for (unsigned int i=1; i<=subject_.length(); i++) {
		b_score_[query_.length()+1][i] = 0;
		b_score_q_[query_.length()+1][i] = 0;
		b_score_s_[query_.length()+1][i] = NINFTY;
	}

	for (unsigned int i=1; i<=query_.length(); i++) {
		b_score_[i][subject_.length()+1] = 0;
		b_score_s_[i][subject_.length()+1] = 0;
		b_score_q_[i][subject_.length()+1] = NINFTY;
	}

	max_score_backward_[query_.length()+1] = 0;

	int query_aa_idx, subject_aa_idx;

	for (unsigned int i=query_.length(); i>=1; i--) {
		for (unsigned int j=subject_.length(); j>=1; j--) {
			b_score_s_[i][j] = Max3(b_score_[i+1][j]+p_mat_->gap_open_,
					b_score_s_[i+1][j]+p_mat_->gap_extend_,
					b_score_q_[i+1][j]+p_mat_->gap_open_);
			b_score_q_[i][j] = Max3(b_score_[i][j+1]+p_mat_->gap_open_,
					b_score_s_[i][j+1]+p_mat_->gap_open_,
					b_score_q_[i][j+1]+p_mat_->gap_extend_);

			query_aa_idx = aa2idx[(int)query_[i-1]];
			subject_aa_idx = aa2idx[(int)subject_[j-1]];
			if (query_aa_idx == NA || subject_aa_idx == NA) {
				fprintf(stderr, "Unknown AA character, query:[%c], subject:[%c]\n", query_[i-1], subject_[j-1]);
				exit(1);
			}
			b_score_[i][j] = Max3(b_score_q_[i+1][j+1], b_score_s_[i+1][j+1], b_score_[i+1][j+1])
					+ p_mat_->matrix_[query_aa_idx][subject_aa_idx];
		}

		max_score_backward_[i] = max(Max3(b_score_q_[i][1], b_score_s_[i][1], b_score_[i][1]), max_score_backward_[i+1]);
	}

	int max_score = NINFTY;

	for (unsigned int i=query_.length(); i>=1; i--) {
		int max3 = Max3(b_score_q_[i][1], b_score_s_[i][1], b_score_[i][1]);
		if (max3 > max_score) {
			max_score = max3;

			query_aligned_start_pos_ = i;
			subject_aligned_start_pos_ = 1;
		}
	}
	for (unsigned int i=subject_.length(); i>=1; i--) {
		int max3 = Max3(b_score_q_[1][i], b_score_s_[1][i], b_score_[1][i]);
		if (max3 > max_score) {
			max_score = max3;

			query_aligned_start_pos_ = 1;
			subject_aligned_start_pos_ = i;
		}
	}

	return max_score;
}

void Tables::ComputeAlignment_Partial(int pos_to)
{
	score_ = Max3(b_score_q_[pos_query_-pos_start_][pos_subject_], b_score_s_[pos_query_-pos_start_][pos_subject_], b_score_[pos_query_-pos_start_][pos_subject_]);
	if (score_ == b_score_s_[pos_query_-pos_start_][pos_subject_]) {
		constraint_status_ = CONSTRAINT_SUBJECT_GAP;
	} else if (score_ == b_score_q_[pos_query_-pos_start_][pos_subject_]) {
		constraint_status_ = CONSTRAINT_QUERY_GAP;
	}

	while (pos_query_ <= pos_to && pos_subject_ <= (int)subject_.length()) {
		if (constraint_status_ == CONSTRAINT_QUERY_GAP) {
			if (score_ == b_score_q_[pos_query_-pos_start_][pos_subject_+1] + p_mat_->gap_extend_) {
				score_ -= p_mat_->gap_extend_;
				constraint_status_ = CONSTRAINT_QUERY_GAP;
			} else {
				score_ -= p_mat_->gap_open_;
				constraint_status_ = CONSTRAINT_NONE;
			}
			query_aligned_.append(1, '-');
			subject_aligned_.append(subject_, pos_subject_-1, 1);
			pos_subject_++;
		} else if (constraint_status_ == CONSTRAINT_SUBJECT_GAP) {
			if (score_ == b_score_s_[pos_query_+1-pos_start_][pos_subject_]+p_mat_->gap_extend_) {
				score_ -= p_mat_->gap_extend_;
				constraint_status_ = CONSTRAINT_SUBJECT_GAP;
			} else {
				score_ -= p_mat_->gap_open_;
				constraint_status_ = CONSTRAINT_NONE;
			}
			query_aligned_.append(query_, pos_query_-1, 1);
			subject_aligned_.append(1, '-');
			pos_query_++;
		} else {	//  CONSTRAINT_NONE
			if (score_ == b_score_q_[pos_query_-pos_start_][pos_subject_+1] + p_mat_->gap_extend_
					|| score_ == b_score_[pos_query_-pos_start_][pos_subject_+1] + p_mat_->gap_open_
					|| score_ == b_score_s_[pos_query_-pos_start_][pos_subject_+1] + p_mat_->gap_open_) {
				if (score_ == b_score_q_[pos_query_-pos_start_][pos_subject_+1] + p_mat_->gap_extend_) {
					score_ -= p_mat_->gap_extend_;
					constraint_status_ = CONSTRAINT_QUERY_GAP;
				} else {
					score_ -= p_mat_->gap_open_;
					constraint_status_ = CONSTRAINT_NONE;
				}
				query_aligned_.append(1, '-');
				subject_aligned_.append(subject_, pos_subject_-1, 1);
				pos_subject_++;
			} else if (score_ == b_score_[pos_query_+1-pos_start_][pos_subject_]+p_mat_->gap_open_
					|| score_ == b_score_q_[pos_query_+1-pos_start_][pos_subject_]+p_mat_->gap_open_
					|| score_ == b_score_s_[pos_query_+1-pos_start_][pos_subject_]+p_mat_->gap_extend_) {
				if (score_ == b_score_s_[pos_query_+1-pos_start_][pos_subject_]+p_mat_->gap_extend_) {
					score_ -= p_mat_->gap_extend_;
					constraint_status_ = CONSTRAINT_SUBJECT_GAP;
				} else {
					score_ -= p_mat_->gap_open_;
					constraint_status_ = CONSTRAINT_NONE;
				}
				query_aligned_.append(query_, pos_query_-1, 1);
				subject_aligned_.append(1, '-');
				pos_query_++;
			} else {
				score_ -= p_mat_->matrix_[aa2idx[(int)query_[pos_query_-1]]][aa2idx[(int)subject_[pos_subject_-1]]];
				constraint_status_ = CONSTRAINT_NONE;
				query_aligned_.append(query_, pos_query_-1, 1);
				subject_aligned_.append(subject_, pos_subject_-1, 1);
				pos_query_++;
				pos_subject_++;
			}
		}
	}

}

void Tables::ComputeAlignment_Full()
{
	int index_aligned = 0;
	int pos_query = query_aligned_start_pos_;
	int pos_subject = subject_aligned_start_pos_;

	int score;
	int constraint = CONSTRAINT_NONE;

	char* query_aligned = new char[query_.length() + subject_.length()+1];
	char* subject_aligned = new char[query_.length() + subject_.length()+1];

	if (pos_query == 1) {
		for (int i=1; i<pos_subject; i++) {
			query_aligned[index_aligned] = '*';
			subject_aligned[index_aligned] = subject_[i-1];//tolower(subject_[i-1]);
			index_aligned++;
		}
	} else if (pos_subject == 1) {
		for (int i=1; i<pos_query; i++) {
			query_aligned[index_aligned] = query_[i-1];//tolower(query_[i-1]);
			subject_aligned[index_aligned] = '*';
			index_aligned++;
		}
	}

	score = Max3(b_score_q_[pos_query][pos_subject], b_score_s_[pos_query][pos_subject], b_score_[pos_query][pos_subject]);
	if (score == b_score_s_[pos_query][pos_subject]) {
		constraint = CONSTRAINT_SUBJECT_GAP;
	} else if (score == b_score_q_[pos_query][pos_subject]) {
		constraint = CONSTRAINT_QUERY_GAP;
	}

	while (pos_query <= (int)query_.length() && pos_subject <= (int)subject_.length()) {
		if (constraint == CONSTRAINT_QUERY_GAP) {
			if (score == b_score_q_[pos_query][pos_subject+1] + p_mat_->gap_extend_) {
				score -= p_mat_->gap_extend_;
				constraint = CONSTRAINT_QUERY_GAP;
			} else {
				score -= p_mat_->gap_open_;
				constraint = CONSTRAINT_NONE;
			}
			query_aligned[index_aligned] = '-';
			subject_aligned[index_aligned] = subject_[pos_subject-1];
			index_aligned++;
			pos_subject++;
		} else if (constraint == CONSTRAINT_SUBJECT_GAP) {
			if (score == b_score_s_[pos_query+1][pos_subject]+p_mat_->gap_extend_) {
				score -= p_mat_->gap_extend_;
				constraint = CONSTRAINT_SUBJECT_GAP;
			} else {
				score -= p_mat_->gap_open_;
				constraint = CONSTRAINT_NONE;
			}
			query_aligned[index_aligned] = query_[pos_query-1];
			subject_aligned[index_aligned] = '-';
			index_aligned++;
			pos_query++;
		} else {	//  CONSTRAINT_NONE
			if (score == b_score_q_[pos_query][pos_subject+1] + p_mat_->gap_extend_
					|| score == b_score_[pos_query][pos_subject+1] + p_mat_->gap_open_
					|| score == b_score_s_[pos_query][pos_subject+1] + p_mat_->gap_open_) {
				if (score == b_score_q_[pos_query][pos_subject+1] + p_mat_->gap_extend_) {
					score -= p_mat_->gap_extend_;
					constraint = CONSTRAINT_QUERY_GAP;
				} else {
					score -= p_mat_->gap_open_;
					constraint = CONSTRAINT_NONE;
				}
				query_aligned[index_aligned] = '-';
				subject_aligned[index_aligned] = subject_[pos_subject-1];
				index_aligned++;
				pos_subject++;
			} else if (score == b_score_[pos_query+1][pos_subject]+p_mat_->gap_open_
					|| score == b_score_q_[pos_query+1][pos_subject]+p_mat_->gap_open_
					|| score == b_score_s_[pos_query+1][pos_subject]+p_mat_->gap_extend_) {
				if (score == b_score_s_[pos_query+1][pos_subject]+p_mat_->gap_extend_) {
					score -= p_mat_->gap_extend_;
					constraint = CONSTRAINT_SUBJECT_GAP;
				} else {
					score -= p_mat_->gap_open_;
					constraint = CONSTRAINT_NONE;
				}
				query_aligned[index_aligned] = query_[pos_query-1];
				subject_aligned[index_aligned] = '-';
				index_aligned++;
				pos_query++;
			} else {
				score -= p_mat_->matrix_[aa2idx[(int)query_[pos_query-1]]][aa2idx[(int)subject_[pos_subject-1]]];
				constraint = CONSTRAINT_NONE;
				query_aligned[index_aligned] = query_[pos_query-1];
				subject_aligned[index_aligned] = subject_[pos_subject-1];
				index_aligned++;
				pos_query++;
				pos_subject++;
			}
		}
	}

	while (pos_query <= (int)query_.length()) {
		query_aligned[index_aligned] = query_[pos_query-1];//tolower(query_[pos_query-1]);
		subject_aligned[index_aligned] = '*';
		index_aligned++;
		pos_query++;
	}

	while (pos_subject <= (int)subject_.length()) {
		query_aligned[index_aligned] = '*';
		subject_aligned[index_aligned] = subject_[pos_subject-1];//tolower(subject_[pos_subject-1]);
		index_aligned++;
		pos_subject++;
	}

	query_aligned[index_aligned] = '\0';
	subject_aligned[index_aligned] = '\0';

	query_aligned_.assign(query_aligned);
	subject_aligned_.assign(subject_aligned);

	delete[] query_aligned;
	delete[] subject_aligned;
}

double Tables::GetIdentity()
{
	int count = 0;
	int total = 0;
	for (unsigned int i=0; i<query_aligned_.length(); i++) {
		if (query_aligned_[i] != '*' && subject_aligned_[i] != '*') {
			if (query_aligned_[i] == subject_aligned_[i]) {
				count++;
			}
			total++;
		}
	}

	return (double)count/total;
}

// get alignment & store aligned sequences using BackwardTable
void Tables::ComputeAlignment()
{
	if (moving_table_) {
		this->ComputeAlignment_Partial();
		pos_start_ = 0;
		FillBackwardTable_Partial(0, stud_size_);
	} else {
		this->ComputeAlignment_Full();
	}
}

void Tables::ComputeAlignment_Partial()
{
	// The first table should already be computed
	pos_query_ = query_aligned_start_pos_;
	pos_subject_ = subject_aligned_start_pos_;
	constraint_status_ = CONSTRAINT_NONE;

	this->query_aligned_.clear();
	this->subject_aligned_.clear();

	if (pos_query_ == 1) {
		query_aligned_.append(pos_subject_-1, '*');
		subject_aligned_.append(subject_, 0, pos_subject_-1);// maybe better to make lower cases
	} else if (pos_subject_ == 1) {
		query_aligned_.append(query_, 0, pos_query_-1);
		subject_aligned_.append(pos_query_-1, '*');
	}

	int pos_from = stud_size_*(int)floor((double)pos_query_/stud_size_);
	int pos_to = pos_from + stud_size_;

	while (pos_query_ <= (int)query_.length() && pos_subject_ <= (int)subject_.length()) {
		if (pos_to > (int)query_.length()+1) {
			pos_to = query_.length()+1;
		}
		pos_start_ = pos_from;
		this->FillBackwardTable_Partial(pos_from, pos_to);
		this->ComputeAlignment_Partial(pos_to-1);
		pos_from = pos_to-1;
		pos_to += stud_size_;
	}

	if (pos_query_ <= (int)query_.length()) {
		query_aligned_.append(query_, pos_query_-1, query_.length()-pos_query_+1);
		subject_aligned_.append(query_.length()-pos_query_+1, '*');
	}

	if (pos_subject_ <= (int)subject_.length()) {
		query_aligned_.append(subject_.length()-pos_subject_+1, '*');
		subject_aligned_.append(subject_, pos_subject_-1, subject_.length()-pos_subject_+1);
	}

}

void Tables::PrintAlignment(int segment_len)
{
	int len = query_aligned_.length();

	string med(len, '\0');

	for (int i=0; i<len; i++) {
		if (query_aligned_[i] == '-' || subject_aligned_[i] == '-' || query_aligned_[i] == '*' || subject_aligned_[i] == '*') {
			med[i] = ' ';
		} else if (query_aligned_[i] == subject_aligned_[i]) {
			med[i] = query_aligned_[i];
		} else if (p_mat_->matrix_[aa2idx[(int)query_aligned_[i]]][aa2idx[(int)subject_aligned_[i]]] > 0) {
			med[i] = '+';
		} else {
			med[i] = ' ';
		}
	}

	for (int i=0; i<len; i+=segment_len) {
		fprintf(stdout, "%s\n", query_aligned_.substr(i, segment_len).c_str());
		fprintf(stdout, "%s\n", med.substr(i, segment_len).c_str());
		fprintf(stdout, "%s\n", subject_aligned_.substr(i, segment_len).c_str());
		fprintf(stdout, "\n");
	}
}

void Tables::FillForwardTable_Partial(int pos_from, int pos_to)
{
#ifdef DEBUG
	fprintf(stderr, "...FillForwardTable_Partial(%d,%d)\n", pos_from, pos_to);
#endif
	int query_aa_idx, subject_aa_idx;

	if (pos_to > (int)query_.length()) {
		pos_to = query_.length();
	}

	for (int i=pos_from; i<=pos_to; i++) {
		if (i == 0) {
			f_score_[0][0] = 0;
			f_score_s_[0][0] = NINFTY;
			f_score_q_[0][0] = NINFTY;
			for (unsigned int j=1; j<=subject_.length(); j++) {
				f_score_[0][j] = 0;
				f_score_q_[0][j] = 0;
				f_score_s_[0][j] = NINFTY;
			}
			max_score_forward_[0] = 0;
		} else {
			f_score_[i-pos_start_][0] = 0;
			f_score_s_[i-pos_start_][0] = 0;
			f_score_q_[i-pos_start_][0] = NINFTY;

			for (unsigned int j=1; j<=subject_.length(); j++) {
				f_score_s_[i-pos_start_][j] = Max3(f_score_[i-pos_start_-1][j]+p_mat_->gap_open_,
						f_score_s_[i-pos_start_-1][j]+p_mat_->gap_extend_,
						f_score_q_[i-pos_start_-1][j]+p_mat_->gap_open_);
				f_score_q_[i-pos_start_][j] = Max3(f_score_[i-pos_start_][j-1]+p_mat_->gap_open_,
						f_score_s_[i-pos_start_][j-1]+p_mat_->gap_open_,
						f_score_q_[i-pos_start_][j-1]+p_mat_->gap_extend_);

				query_aa_idx = aa2idx[(int)query_[i-1]];
				subject_aa_idx = aa2idx[(int)subject_[j-1]];
				if (query_aa_idx == NA || subject_aa_idx == NA) {
					fprintf(stderr, "Unknown AA character, query:[%c], subject:[%c]\n", query_[i-1], subject_[j-1]);
					exit(1);
				}
				f_score_[i-pos_start_][j] = Max3(f_score_q_[i-pos_start_-1][j-1], f_score_s_[i-pos_start_-1][j-1], f_score_[i-pos_start_-1][j-1])
						+ p_mat_->matrix_[query_aa_idx][subject_aa_idx];
			}

			max_score_forward_[i] = max(Max3(f_score_q_[i-pos_start_][subject_.length()], f_score_s_[i-pos_start_][subject_.length()], f_score_[i-pos_start_][subject_.length()]), max_score_forward_[i-1]);
		}
	}

	return;
}


int Tables::FillBackwardTable_Partial(int pos_from, int pos_to)
{
#ifdef DEBUG
	fprintf(stderr, "...FillBackwardTable_Partial(%d,%d)\n", pos_from, pos_to);
#endif
	// x : query
	// y : subject
	int query_aa_idx, subject_aa_idx;

	// copy from stud
	if (pos_to == (int)query_.length()+1) {
		b_score_[pos_to-pos_start_][subject_.length()+1] = 0;
		b_score_s_[pos_to-pos_start_][subject_.length()+1] = NINFTY;
		b_score_q_[pos_to-pos_start_][subject_.length()+1] = NINFTY;
		for (unsigned int j=1; j<=subject_.length(); j++) {
			b_score_[pos_to-pos_start_][j] = 0;
			b_score_q_[pos_to-pos_start_][j] = 0;
			b_score_s_[pos_to-pos_start_][j] = NINFTY;
		}
	} else {
		memcpy(b_score_[pos_to-pos_start_], b_score_stud_[pos_to/stud_size_-1], sizeof(int)*len_y_);
		memcpy(b_score_q_[pos_to-pos_start_], b_score_stud_q_[pos_to/stud_size_-1], sizeof(int)*len_y_);
		memcpy(b_score_s_[pos_to-pos_start_], b_score_stud_s_[pos_to/stud_size_-1], sizeof(int)*len_y_);
	}

	if (pos_from == 0) {
		pos_from = 1;
	}

	for (int i=pos_to-1; i>=pos_from; i--) {
		b_score_[i-pos_start_][subject_.length()+1] = 0;
		b_score_s_[i-pos_start_][subject_.length()+1] = 0;
		b_score_q_[i-pos_start_][subject_.length()+1] = NINFTY;

		for (unsigned int j=subject_.length(); j>=1; j--) {
			b_score_s_[i-pos_start_][j] = Max3(b_score_[i-pos_start_+1][j]+p_mat_->gap_open_,
					b_score_s_[i-pos_start_+1][j]+p_mat_->gap_extend_,
					b_score_q_[i-pos_start_+1][j]+p_mat_->gap_open_);
			b_score_q_[i-pos_start_][j] = Max3(b_score_[i-pos_start_][j+1]+p_mat_->gap_open_,
					b_score_s_[i-pos_start_][j+1]+p_mat_->gap_open_,
					b_score_q_[i-pos_start_][j+1]+p_mat_->gap_extend_);

			query_aa_idx = aa2idx[(int)query_[i-1]];
			subject_aa_idx = aa2idx[(int)subject_[j-1]];
			if (query_aa_idx == NA || subject_aa_idx == NA) {
				fprintf(stderr, "Unknown AA character, query:[%c], subject:[%c]\n", query_[i-1], subject_[j-1]);
				exit(1);
			}
			b_score_[i-pos_start_][j] = Max3(b_score_q_[i-pos_start_+1][j+1], b_score_s_[i-pos_start_+1][j+1], b_score_[i-pos_start_+1][j+1])
					+ p_mat_->matrix_[query_aa_idx][subject_aa_idx];
		}


		// for the 1st table
		if (pos_from == 1) {
			max_score_backward_[i] = max(Max3(b_score_q_[i][1], b_score_s_[i][1], b_score_[i][1]), max_score_backward_[i+1]);


		}
	}

	int max_score = NINFTY;
	if (pos_from == 1) {
		max_score = max_score_backward_[1];
		for (unsigned int i=subject_.length(); i>=1; i--) {
			int max3 = Max3(b_score_q_[1][i], b_score_s_[1][i], b_score_[1][i]);
			if (max3 > max_score) {
				max_score = max3;

				query_aligned_start_pos_ = 1;
				subject_aligned_start_pos_ = i;
			}
		}

		if (max_score == max_score_backward_[1]) {
			int i = 2;
			while (max_score == max_score_backward_[i]) {
				i++;
			}
			query_aligned_start_pos_ = i-1;
			subject_aligned_start_pos_ = 1;
		}

	}

	return max_score;
}

int Tables::FillBackwardTable_Stud()
{
#ifdef DEBUG
	fprintf(stderr, "...FillBackwardTable_Stud\n");
#endif
	// x : query
	// y : subject
	b_score_[(query_.length()+1)%2][subject_.length()+1] = 0;
	b_score_s_[(query_.length()+1)%2][subject_.length()+1] = NINFTY;
	b_score_q_[(query_.length()+1)%2][subject_.length()+1] = NINFTY;

	for (unsigned int i=1; i<=subject_.length(); i++) {
		b_score_[(query_.length()+1)%2][i] = 0;
		b_score_q_[(query_.length()+1)%2][i] = 0;
		b_score_s_[(query_.length()+1)%2][i] = NINFTY;
	}

	max_score_backward_[query_.length()+1] = 0;

	int query_aa_idx, subject_aa_idx;

	for (unsigned int i=query_.length(); i>=stud_size_; i--) {
		b_score_[i%2][subject_.length()+1] = 0;
		b_score_s_[i%2][subject_.length()+1] = 0;
		b_score_q_[i%2][subject_.length()+1] = NINFTY;
		for (unsigned int j=subject_.length(); j>=1; j--) {
			b_score_s_[i%2][j] = Max3(b_score_[(i+1)%2][j]+p_mat_->gap_open_,
					b_score_s_[(i+1)%2][j]+p_mat_->gap_extend_,
					b_score_q_[(i+1)%2][j]+p_mat_->gap_open_);
			b_score_q_[i%2][j] = Max3(b_score_[i%2][j+1]+p_mat_->gap_open_,
					b_score_s_[i%2][j+1]+p_mat_->gap_open_,
					b_score_q_[i%2][j+1]+p_mat_->gap_extend_);

			query_aa_idx = aa2idx[(int)query_[i-1]];
			subject_aa_idx = aa2idx[(int)subject_[j-1]];
			if (query_aa_idx == NA || subject_aa_idx == NA) {
				fprintf(stderr, "Unknown AA character, query:[%c], subject:[%c]\n", query_[i-1], subject_[j-1]);
				exit(1);
			}
			b_score_[i%2][j] = Max3(b_score_q_[(i+1)%2][j+1], b_score_s_[(i+1)%2][j+1], b_score_[(i+1)%2][j+1])
					+ p_mat_->matrix_[query_aa_idx][subject_aa_idx];
		}

		max_score_backward_[i] = max(Max3(b_score_q_[i%2][1], b_score_s_[i%2][1], b_score_[i%2][1]), max_score_backward_[i+1]);

		if (i%stud_size_ == 0) {
			memcpy(b_score_stud_[i/stud_size_-1], b_score_[i%2], sizeof(int)*len_y_);
			memcpy(b_score_stud_q_[i/stud_size_-1], b_score_q_[i%2], sizeof(int)*len_y_);
			memcpy(b_score_stud_s_[i/stud_size_-1], b_score_s_[i%2], sizeof(int)*len_y_);
		}
	}

	return 0;
}

int Tables::GetRefAlignmentScore()
{
	int max_score = NINFTY;
	int query_end_pos = 0, subject_end_pos = 0;

	for (unsigned int i=1; i<=query_.length(); i++) {
		int max3 = Max3(f_score_q_[i][subject_.length()], f_score_s_[i][subject_.length()], f_score_[i][subject_.length()]);
		if (max3 > max_score) {
			max_score = max3;
			query_end_pos = subject_.length();
			subject_end_pos = i;
		}
	}
	for (unsigned int i=1; i<=subject_.length(); i++) {
		int max3 = Max3(f_score_q_[query_.length()][i], f_score_s_[query_.length()][i], f_score_[query_.length()][i]);
		if (max3 > max_score) {
			max_score = max3;
			query_end_pos = i;
			subject_end_pos = query_.length();
		}
	}

	return max_score;
}

int Tables::GetRefAlignmentScore_B()
{
	int max_score = NINFTY;

	for (unsigned int i=1; i<=query_.length(); i++) {
		int max3 = Max3(b_score_q_[i][1], b_score_s_[i][1], b_score_[i][1]);
		if (max3 > max_score) {
			max_score = max3;
		}
	}
	for (unsigned int i=1; i<=subject_.length(); i++) {
		int max3 = Max3(b_score_q_[1][i], b_score_s_[1][i], b_score_[1][i]);
		if (max3 > max_score) {
			max_score = max3;
		}
	}

	return max_score;
}


int Tables::GetVarAlignmentScore(int pos_left, int pos_right, string str_ins)
{
	while (pos_right > pos_end_) {
		// move table
		if (pos_left <= pos_end_) {
			for (int i=pos_left; i<=pos_end_; i++) {
				memcpy(f_score_[i-pos_left], f_score_[i-pos_start_], sizeof(int)*len_y_);
				memcpy(f_score_s_[i-pos_left], f_score_s_[i-pos_start_], sizeof(int)*len_y_);
				memcpy(f_score_q_[i-pos_left], f_score_q_[i-pos_start_], sizeof(int)*len_y_);
			}
			int len_copied = pos_end_ - pos_left +1;
			pos_start_ = pos_left;
			pos_end_ = pos_end_+stud_size_;
			if (pos_end_ > (int)query_.length()+1) {
				pos_end_ = query_.length()+1;
			}
			this->FillForwardTable_Partial(pos_start_+len_copied, pos_end_);
			//this->FillBackwardTable_Partial(pos_start_, pos_end_);
		} else {
			memcpy(f_score_[0], f_score_[pos_end_-pos_start_], sizeof(int)*len_y_);
			memcpy(f_score_s_[0], f_score_s_[pos_end_-pos_start_], sizeof(int)*len_y_);
			memcpy(f_score_q_[0], f_score_q_[pos_end_-pos_start_], sizeof(int)*len_y_);
			pos_start_ = pos_end_;
			pos_end_ = pos_end_+stud_size_;
			if (pos_end_ > (int)query_.length()+1) {
				pos_end_ = query_.length()+1;
			}
			this->FillForwardTable_Partial(pos_start_+1, pos_end_);
			//this->FillBackwardTable_Partial(pos_start_, pos_end_);
		}

		if (pos_right <= pos_end_) {
			this->FillBackwardTable_Partial(pos_start_, pos_end_);
		}
	}



	int max_score = NINFTY;

	int len_added = str_ins.length();

	if (len_added == 0) {
		// deletion
		for (unsigned int j=0; j<=subject_.length(); j++) {
			int max3 = Max3(
					Max3(f_score_[pos_left-pos_start_][j], f_score_q_[pos_left-pos_start_][j], f_score_s_[pos_left-pos_start_][j])
					+ Max3(b_score_[pos_right-pos_start_][j+1], b_score_q_[pos_right-pos_start_][j+1], b_score_s_[pos_right-pos_start_][j+1]),
					f_score_q_[pos_left-pos_start_][j] + b_score_q_[pos_right-pos_start_][j+1] + (p_mat_->gap_extend_ - p_mat_->gap_open_),
					f_score_s_[pos_left-pos_start_][j] + b_score_s_[pos_right-pos_start_][j+1] + (p_mat_->gap_extend_ - p_mat_->gap_open_));

			if (max3 > max_score) {
				max_score = max3;
			}
		}

		if (max_score_forward_[pos_left] > max_score) {
			max_score = max_score_forward_[pos_left];
		}
		if (max_score_backward_[pos_right] > max_score) {
			max_score = max_score_backward_[pos_right];
		}

	} else {
		int** t_score = new int*[len_added];
		int** t_score_s = new int*[len_added];
		int** t_score_q = new int*[len_added];
		for (int i=0; i<len_added; i++) {
			t_score[i] = new int[subject_.length()+2];
			t_score_s[i] = new int[subject_.length()+2];
			t_score_q[i] = new int[subject_.length()+2];
		}

		int query_aa_idx, subject_aa_idx;

		// boundary
		for (int i=0; i<len_added; i++) {
			t_score[i][0] = 0;
			t_score_s[i][0] = 0;
			t_score_q[i][0] = NINFTY;
		}

		// 1st column
		for (unsigned int j=1; j<=subject_.length(); j++) {
			char var_aa = str_ins[0];

			t_score_s[0][j] = Max3(f_score_[pos_left-pos_start_][j]+p_mat_->gap_open_,
					f_score_s_[pos_left-pos_start_][j]+p_mat_->gap_extend_,
					f_score_q_[pos_left-pos_start_][j]+p_mat_->gap_open_);
			t_score_q[0][j] = Max3(t_score[0][j-1]+p_mat_->gap_open_,
					t_score_s[0][j-1]+p_mat_->gap_open_,
					t_score_q[0][j-1]+p_mat_->gap_extend_);

			query_aa_idx = aa2idx[(int)var_aa];
			subject_aa_idx = aa2idx[(int)subject_[j-1]];
			if (query_aa_idx == NA || subject_aa_idx == NA) {
				fprintf(stderr, "Unknown AA character, query:[%c], subject:[%c]\n", var_aa, subject_[j-1]);
				exit(1);
			}
			t_score[0][j] = Max3(f_score_q_[pos_left-pos_start_][j-1], f_score_s_[pos_left-pos_start_][j-1], f_score_[pos_left-pos_start_][j-1])
					+ p_mat_->matrix_[query_aa_idx][subject_aa_idx];
		}

		// from 2nd to last column
		for (int i=1; i<len_added; i++) {
			for (unsigned int j=1; j<=subject_.length(); j++) {
				t_score_s[i][j] = Max3(t_score[i-1][j]+p_mat_->gap_open_,
						t_score_s[i-1][j]+p_mat_->gap_extend_,
						t_score_q[i-1][j]+p_mat_->gap_open_);
				t_score_q[i][j] = Max3(t_score[i][j-1]+p_mat_->gap_open_,
						t_score_s[i][j-1]+p_mat_->gap_open_,
						t_score_q[i][j-1]+p_mat_->gap_extend_);

				query_aa_idx = aa2idx[(int)str_ins[i]];
				subject_aa_idx = aa2idx[(int)subject_[j-1]];
				if (query_aa_idx == NA || subject_aa_idx == NA) {
					fprintf(stderr, "Unknown AA character, query:[%c], subject:[%c]\n", str_ins[i], subject_[j-1]);
					exit(1);
				}
				t_score[i][j] = Max3(t_score_q[i-1][j-1], t_score_s[i-1][j-1], t_score[i-1][j-1])
						+ p_mat_->matrix_[query_aa_idx][subject_aa_idx];
			}
		}


		// find maximum
		//
		for (unsigned int j=0; j<=subject_.length(); j++) {
			int max3 = Max3(
					Max3(t_score[len_added-1][j], t_score_q[len_added-1][j], t_score_s[len_added-1][j])
					+ Max3(b_score_[pos_right-pos_start_][j+1], b_score_q_[pos_right-pos_start_][j+1], b_score_s_[pos_right-pos_start_][j+1]),
					t_score_q[len_added-1][j] + b_score_q_[pos_right-pos_start_][j+1] + (p_mat_->gap_extend_ - p_mat_->gap_open_),
					t_score_s[len_added-1][j] + b_score_s_[pos_right-pos_start_][j+1] + (p_mat_->gap_extend_ - p_mat_->gap_open_));

			if (max3 > max_score) {
				max_score = max3;
			}
		}

		for (int i=0; i<len_added; i++) {
			int max3 = Max3(t_score_q[i][subject_.length()], t_score_s[i][subject_.length()], t_score[i][subject_.length()]);
			if (max3 > max_score) {
				max_score = max3;
			}
		}

		if (max_score_forward_[pos_left] > max_score) {
			max_score = max_score_forward_[pos_left];
		}
		if (max_score_backward_[pos_right] > max_score) {
			max_score = max_score_backward_[pos_right];
		}

		for (int i=0; i<len_added; i++) {
			delete[] t_score[i];
			delete[] t_score_s[i];
			delete[] t_score_q[i];
		}
		delete[] t_score;
		delete[] t_score_s;
		delete[] t_score_q;
	}

	return max_score;
}
