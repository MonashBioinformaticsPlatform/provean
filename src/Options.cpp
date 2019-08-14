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
// Name        : Options.cpp
// Author      : Yongwook Choi
//============================================================================

#include "Options.h"
#include <getopt.h>

char usage[] =
		"USAGE:\
		\n  provean [Options]\
		\n\nOptions:\
		\n  -q <string>\
		\n  -d <string>\
		\n  -v <string>\
		\n  --psiblast <string>\
		\n  --cdhit <string>\
		\n  --blastdbcmd <string>\
		\n  --quiet\
		\n  --supporting_set <string>\
		\n  --save_supporting_set <string>\
		\n  --subject_sequences <string>\
		\n  --num_threads <integer, >=1> (default=1)\
		\n  --tmp_dir <string>\
		\n";


Options::Options() {
	flag_all_sap_ = 0;
	flag_quiet_ = 0;
	flag_save_blastout_ = 0;
	flag_save_supporting_set_ = 0;
	num_cluster_ = 30;
	cluster_threshold_ = 0.75;
	score_matrix_.SetScoreMatrix("BLOSUM62", -10, -1);
	num_threads_ = 1;
	tmp_dir_given_ = 0;

	this->percent_length_query_ = 30;
	this->percent_length_rep_ = 30;
}

Options::~Options() {
	// TODO Auto-generated destructor stub
}


// return values
// 0 for successful setting
//
int Options::SetOptions(int argc, char** argv) {
	int c;

	while (1)
	{
		static struct option long_options[] =
		{
				{"quiet", no_argument, &flag_quiet_, 1},
				{"all_sap", no_argument, &flag_all_sap_, 1},
				{"blastout", required_argument, 0, 1},
				{"gap_open", required_argument, 0, 2},
				{"gap_extend", required_argument, 0, 3},
				{"num_cluster", required_argument, 0, 4},
				{"clustering_threshold", required_argument, 0, 5},
				{"save_blastout", required_argument, 0, 6},
				{"psiblast", required_argument, 0, 7},
				{"cdhit", required_argument, 0, 8},
				{"blastdbcmd", required_argument, 0, 9},
				{"save_supporting_set", required_argument, 0, 10},
				{"supporting_set", required_argument, 0, 11},
				{"subject_sequences", required_argument, 0, 12},
				{"num_threads", required_argument, 0, 13},
				{"tmp_dir", required_argument, 0, 14},
				{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long(argc, argv, "q:d:v:",
				long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			break;
		case 1:
			blast_output_file_name_ = optarg;
			break;
		case 2:
			this->score_matrix_.gap_open_ = -atoi(optarg);
			break;
		case 3:
			this->score_matrix_.gap_extend_ = -atoi(optarg);
			break;
		case 4:
			this->num_cluster_ = atoi(optarg);
			break;
		case 5:
			this->cluster_threshold_ = (double)atoi(optarg)/100.0;
			break;
		case 6:
			flag_save_blastout_ = 1;
			save_blast_output_file_name_ = optarg;
			break;
		case 7:
			psiblast_command_ = optarg;
			break;
		case 8:
			cdhit_command_ = optarg;
			break;
		case 9:
			blastdbcmd_command_ = optarg;
			break;
		case 10:
			flag_save_supporting_set_ = 1;
			save_supporting_set_file_name_ = optarg;
			break;
		case 11:
			supporting_set_file_name_ = optarg;
			break;
		case 12:
			subject_sequences_fasta_file_name_ = optarg;
			break;
		case 13:
			num_threads_ = atoi(optarg);
			break;
		case 14:
			tmp_dir_ = optarg;
			tmp_dir_given_ = 1;
			break;
		case 'q':
			query_file_name_ = optarg;
			break;
		case 'd':
			blast_db_file_name_ = optarg;
			break;
		case 'v':
			variation_ = optarg;
			break;
		case 'b':
			blast_output_file_name_ = optarg;
			break;

		default: /* '?' */
			this->PrintUsage();
			exit(EXIT_FAILURE);
		}
	}

	if (query_file_name_.empty()
			|| blast_db_file_name_.empty()
			|| psiblast_command_.empty()
			|| cdhit_command_.empty()
			|| (variation_.empty() && !this->flag_all_sap_)) {
		this->PrintUsage();
		exit(EXIT_FAILURE);
	}

	if (num_threads_ < 1) {
		this->PrintUsage();
		exit(EXIT_FAILURE);
	}

	return 0;
}

void Options::PrintOptions(FILE* out)
{
	if (!this->flag_quiet_) {
		fprintf(out, "\n## Parameters ##\n");
	}

	fprintf(out, "# Query sequence file:\t%s\n", query_file_name_.data());
	fprintf(out, "# Variation file:\t%s\n", variation_.data());
	fprintf(out, "# Protein database:\t%s\n", blast_db_file_name_.data());

	if (!this->flag_quiet_) {
//		fprintf(out, "# Precomputed BLAST output file (optional):\t");
//		if (blast_output_file_name_.empty()) {
//			fprintf(out, "Not provided\n");
//		} else {
//			fprintf(out, "%s\n", blast_output_file_name_.data());
//		}
//
//		fprintf(out, "# BLAST output file for storing (optional):\t");
//		if (save_blast_output_file_name_.empty()) {
//			fprintf(out, "Not provided\n");
//		} else {
//			fprintf(out, "%s\n", save_blast_output_file_name_.data());
//		}

		fprintf(out, "# Supporting sequence set file (optional):\t");
		if (supporting_set_file_name_.empty()) {
			fprintf(out, "Not provided\n");
		} else {
			fprintf(out, "%s\n", supporting_set_file_name_.data());
		}

		fprintf(out, "# Supporting sequence set file for storing (optional):\t");
		if (save_supporting_set_file_name_.empty()) {
			fprintf(out, "Not provided\n");
		} else {
			fprintf(out, "%s\n", save_supporting_set_file_name_.data());
		}

		fprintf(out, "# Substitution matrix:\t%s\n", score_matrix_.matrix_name_.data());
		fprintf(out, "# Gap costs:\t%d, %d\n", -score_matrix_.gap_open_, -score_matrix_.gap_extend_);
		fprintf(out, "# Clustering threshold:\t%.3f\n", cluster_threshold_);
		fprintf(out, "# Maximum number of clusters:\t%d\n", num_cluster_);
		fprintf(out, "\n");
	}
}

void Options::PrintUsage()
{
	fprintf(stderr, "%s", usage);
}

