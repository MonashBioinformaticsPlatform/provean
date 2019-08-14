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
// Name        : Options.h
// Author      : Yongwook Choi
//============================================================================

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include "Common.h"
#include "ScoreMatrix.h"

using namespace std;

class Options {
public:
	Options();
	virtual ~Options();

	double cluster_threshold_;
	ScoreMatrix score_matrix_;
	string blast_output_file_name_;
	string query_file_name_;
	string blast_db_file_name_;	// for BLAST
	string variation_;
	string psiblast_command_;
	string cdhit_command_;
	string blastdbcmd_command_;
	string supporting_set_file_name_;	// supporting sequence set
	int flag_save_supporting_set_;	// save supporting sequence set
	string save_supporting_set_file_name_;
	int num_cluster_;
	int flag_quiet_;
	int flag_all_sap_;
	int flag_save_blastout_;
	int num_threads_;
	string save_blast_output_file_name_;
	string tmp_dir_;
	int tmp_dir_given_;

	int percent_length_query_;
	int percent_length_rep_;

	int SetOptions(int argc, char** argv);
	void PrintOptions(FILE* out);
	void PrintUsage();

	string subject_sequences_fasta_file_name_;

};

#endif /* OPTIONS_H_ */
