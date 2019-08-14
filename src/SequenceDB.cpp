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
// Name        : SequenceDB.cpp
// Author      : Yongwook Choi
//============================================================================


#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "Common.h"
#include "SequenceDB.h"

#define NINFTY -100000

#define CONSTRAINT_NONE 0
#define CONSTRAINT_SUBJECT_GAP 1
#define CONSTRAINT_QUERY_GAP 2

using namespace std;

SequenceDB::SequenceDB() {
	max_seq_len_ = 0;

	p_tables_ = new Tables;
	num_seqs_in_cluster_ = NULL;
	num_variants_ = 0;

	b_supporting_itself = false;
	p_options_ = NULL;
}

SequenceDB::~SequenceDB() {
	if (!temp_dir_.empty()) {
		char rm_temp_dir_command[BUF_SIZE_MED];
		sprintf(rm_temp_dir_command, "rm -rf %s", temp_dir_.c_str());
		int ret;
		ret = system(rm_temp_dir_command);
		if (ret != 0) {
			fprintf(stderr, "removing temporary directory failed (%d)\n", ret);
		}
	}
}

void SequenceDB::SetLoggingOption(bool b)
{
	this->b_logging_ = b;
}


bool SequenceDB::IsTooShortComparedToQuery(Sequence* pSeq)
{
	if (pSeq->seq_.length()*100 <= this->query_seq_.seq_.length()*this->p_options_->percent_length_query_) {
		return true;
	} else {
		return false;
	}
}

//
// return : number of sequences added
//
int SequenceDB::SetSequencesFromBlastOut(string blast_file, string blast_db_file, int max_num_seqs)
{
	if (blast_file.empty())
		return -1;

	FILE* fp = fopen(blast_file.c_str(), "r");

	if (fp == NULL) {
		cerr << "cannot open: " << blast_file << endl;
		return -1;
	}

	subject_seqs_.clear();

	char buf[BUF_SIZE_MED];

	if (b_logging_) {
		Log("retrieving subject sequences...\n", true);
	}

	char* db_id;
	char* gi;
	char* strptr;
	int num_seqs = 0;

	// get subject seqs
	while (!feof(fp)) {
		if (fgets(buf, BUF_SIZE_MED, fp) == NULL)
			break;

		db_id = strtok(buf, "\t\n");
		gi = strtok(NULL, "\t\n");
		char id[BUF_SIZE_MED]="";

		if (db_id == NULL || gi == NULL) {
			Log("Parsing Error\n", true);
			return -1;
		}

//		if (atoi(gi) == 0) {
//			// use DB ID
//			strcpy(id, db_id);
//		} else {
//			// use GI number
//			strcpy(id, "gi|");
//			strcpy(id+3, gi);
//		}
		strcpy(id, db_id);

		if (db_id != NULL && gi != NULL) {
			if (!subject_seqs_.empty() && subject_seqs_.back().id_.compare(id) == 0) {
				continue;
			}
			Sequence seq;

			strptr = strtok(NULL, "\t\n");

			seq.e_value_ = atof(strptr);

			strptr = strtok(NULL, "\t\n");
			seq.bit_score_ = atof(strptr);

			if (seq.SetSequence(id, blast_db_file, p_options_->blastdbcmd_command_) == false) {
				Log("No subject sequence info\n", true);
				fprintf(stdout, "%s\n", id);
				continue;
			}

			if (this->IsTooShortComparedToQuery(&seq)) {
				cerr << "[too short sequence]Id: " << seq.id_ << ", length: " << seq.seq_.length() << endl;
				continue;
			}

			subject_seqs_.push_back(seq);
			if (max_seq_len_ < seq.seq_.length()) {
				max_seq_len_ = seq.seq_.length();
			}
			num_seqs++;
			if (num_seqs == max_num_seqs) {
				break;
			}
		} else {
			Log("Wrong entry in blast out\n", true);
			break;
		}
	}

	fclose(fp);

	return num_seqs;
}

//
// return : number of sequences added
//
int SequenceDB::SetSequenceInfoFromBlastOut(string blast_file, int max_num_seqs)
{
	if (blast_file.empty())
		return -1;

	FILE* fp = fopen(blast_file.c_str(), "r");

	if (fp == NULL) {
		cerr << "cannot open: " << blast_file << endl;
		return -1;
	}

	subject_seqs_.clear();

	char buf[BUF_SIZE_MED];

	if (b_logging_) {
		Log("retrieving subject sequence information...\n", true);
	}

	char* db_id;
	char* gi;
	char* strptr;
	int num_seqs = 0;

	// get subject seqs
	while (!feof(fp)) {
		if (fgets(buf, BUF_SIZE_MED, fp) == NULL)
			break;

		db_id = strtok(buf, "\t\n");
		gi = strtok(NULL, "\t\n");
		char id[BUF_SIZE_MED]="";

		if (db_id == NULL || gi == NULL) {
			Log("Parsing Error\n", true);
			return -1;
		}

		strcpy(id, db_id);

		if (db_id != NULL && gi != NULL) {
			if (!subject_seqs_.empty() && subject_seqs_.back().id_.compare(id) == 0) {
				continue;
			}
			Sequence seq;

			seq.id_.assign(id);
			seq.id_.erase(seq.id_.find_last_not_of(" ")+1);
			strptr = strtok(NULL, "\t\n");

			seq.e_value_ = atof(strptr);

			strptr = strtok(NULL, "\t\n");
			seq.bit_score_ = atof(strptr);

			subject_seqs_.push_back(seq);

			num_seqs++;
			if (num_seqs == max_num_seqs) {
				break;
			}
		} else {
			Log("Wrong entry in blast out\n", true);
			break;
		}
	}

	fclose(fp);

	this->UpdateMaxSeqLength();

	return num_seqs;
}

int SequenceDB::SetSequencesFromFastaFile(string fasta_file)
{
	if (fasta_file.empty())
			return -1;

	FILE* fp = fopen(fasta_file.c_str(), "r");

	if (fp == NULL) {
		cerr << "cannot open: " << fasta_file << endl;
		return -1;
	}

	if (b_logging_) {
		Log("loading subject sequences from a FASTA file...\n", true);
	}

	int num_seq = 0;
	char* strptr;
	char buf[BUF_SIZE_LARGE];
	Sequence seq;
	while (!feof(fp)) {
		if (fgets(buf, BUF_SIZE_LARGE, fp) == NULL)
			break;

		if (buf[0] == '>') {
			if (!seq.id_.empty() && !this->IsTooShortComparedToQuery(&seq)) {
				subject_seqs_.push_back(seq);
				num_seq++;
				if (max_seq_len_ < seq.seq_.length()) {
					max_seq_len_ = seq.seq_.length();
				}
				//seq.seq_.clear();
			}
			// def line
			string def(buf, 100);
			size_t pos = def.find_first_of('\n', 0);
			if (pos != string::npos) {
				seq.def_ = def.substr(0, pos);
			} else {
				seq.def_ = def + " ...";
			}

			string id(buf, 10000);
			size_t pos1 = id.find_first_not_of(" ", 1);
			size_t pos2 = id.find_first_of(" \n", pos1);
			if (pos2 != string::npos) {
				seq.id_ = id.substr(pos1, pos2-pos1);
			} else {
				seq.id_ = id.substr(pos1, 10000-1-pos1);
			}

			while (strchr(buf, '\n') == NULL) {
				fgets(buf, sizeof(buf), fp);
			}

			//strptr = strtok(buf, " \n");
			//seq.id_.assign(++strptr);
			seq.seq_.clear();
		} else {
			strptr = strtok(buf, " \n");
			if (strptr == NULL)
				break;
			seq.seq_.append(strptr);
		}

		if (num_seq > 100000) {
			cerr << fasta_file << " has too many sequences. Please check the file." << endl;
			exit(-1);
		}
	}

	if (!seq.id_.empty() && !this->IsTooShortComparedToQuery(&seq)) {
		subject_seqs_.push_back(seq);
		num_seq++;
		if (max_seq_len_ < seq.seq_.length()) {
			max_seq_len_ = seq.seq_.length();
		}
	}

	fclose(fp);

	return num_seq;
}

bool SequenceDB::SetQuerySequenceFromFastaFile(const char* id, string fasta_file)
{
	if (!p_options_->flag_quiet_) {
		Log("loading query sequence from a FASTA file...\n", true);
	}

	return query_seq_.SetSequenceFromFastaFile(id, fasta_file);
}


int SequenceDB::GetNumberOfVariants()
{
	return this->variants_.size();
}

void SequenceDB::SetOptions(Options* p_options)
{
	this->p_options_ = p_options;
}

void SequenceDB::SetSubjectSequences()
{
	int num_seqs;
	bool b_temp_blast_out_file = false;

	//
	subject_seqs_.clear();

	if (!p_options_->subject_sequences_fasta_file_name_.empty()) {
		// --subject_sequences option is given, sequences in fasta format
		if (this->SetSequencesFromFastaFile(p_options_->subject_sequences_fasta_file_name_) == 0) {
			this->AddQueryAsSupportingSequence();
		} else {
			this->ClusterAndSelectFromFastaFile(p_options_->subject_sequences_fasta_file_name_);
		}

	} else if (!p_options_->supporting_set_file_name_.empty()) {
		// --supporting_set option is given
		this->LoadSupportingSet(p_options_->supporting_set_file_name_);
	} else {
		// no supporting set file
		// need to blast, cluster, and select
		if (p_options_->blast_output_file_name_.empty()) {
			// blast output file is NOT given
			b_temp_blast_out_file = true;
			p_options_->blast_output_file_name_ = temp_dir_ + "/blast_out";
			this->RunBlast(p_options_->query_file_name_, p_options_->blast_db_file_name_, p_options_->blast_output_file_name_);

			if (p_options_->flag_save_blastout_) {
				char cp_command[BUF_SIZE_MED];
				sprintf(cp_command, "cp %s %s", p_options_->blast_output_file_name_.c_str(), p_options_->save_blast_output_file_name_.c_str());
				// will replace system
				int ret;
				ret = system(cp_command);
				if (ret != 0) {
					fprintf(stderr, "saving blastout file failed (%d)\n", ret);
				}
			}
		}

		num_seqs = this->SetSequenceInfoFromBlastOut(p_options_->blast_output_file_name_, 100000);

		if (num_seqs < 0) {
			// error
			fprintf(stderr, "error in creating fasta file from blast output\n");
			return;
		} else if (num_seqs == 0) {
			Log("no subject sequences found!\n", true);
			this->AddQueryAsSupportingSequence();
		} else {
			string tmp_fasta_file = temp_dir_ + "/tmp.fasta";
			this->CreateFastaFileOfSubjectSequences(tmp_fasta_file);
			if (this->ClusterAndSelectFromFastaFile(tmp_fasta_file) == 0) {
				this->AddQueryAsSupportingSequence();
			} else {
				this->SetSequences();
			}
		}

		if (b_temp_blast_out_file) {
			char rm_command[BUF_SIZE_MED];
			sprintf(rm_command, "rm %s", p_options_->blast_output_file_name_.c_str());
			int ret;
			ret = system(rm_command);
			if (ret != 0) {
				fprintf(stderr, "removing temporary files failed (%d)\n", ret);
			}
		}
	}

	if (p_options_->flag_save_supporting_set_) {
		this->SaveSupportingSet(p_options_->save_supporting_set_file_name_);
	}
}

void SequenceDB::SaveSupportingSet(string out_file)
{
	FILE* fp = fopen(out_file.c_str(), "w");
	if (fp == NULL) {
		fprintf(stderr, "file open error: %s\n", out_file.c_str());
		exit(-1);
	}

	if (this->b_supporting_itself) {
		//Log("supporting sequence set was not saved because no supporting sequences found.\n", true);
		Log("scores were computed based on the query sequence itself.\n", true);

		fprintf(fp, "# PROVEAN version:\t%s\n", VERSION);
		fprintf(fp, "# No supporting sequences were found in the database used.\n");
		fprintf(fp, "# The scores were computed based on the query sequence itself.\n");
	} else {
		char tmp_file_ids[BUF_SIZE_MED];

		sprintf(tmp_file_ids, "%s/sss.ids", temp_dir_.c_str());
		FILE* fp_ids = fopen(tmp_file_ids, "w");
		if (!fp_ids) {
			Log("File open error!!\n", true);
			cerr << tmp_file_ids << endl;
			return;
		}


		fprintf(fp, "# PROVEAN version:\t%s\n", VERSION);
		fprintf(fp, "# Sequence_ID\tCluster_ID\tE-value\tBit_score\n");

		for (vector<Sequence>::iterator it_seq=subject_seqs_.begin(); it_seq<subject_seqs_.end(); it_seq++) {
			fprintf(fp, "%s\t%d\t%.0e\t%.2f\n", (*it_seq).id_.c_str(), (*it_seq).cluster_id_, (*it_seq).e_value_, (*it_seq).bit_score_);
			fprintf(fp_ids, "%s\n", (*it_seq).id_.c_str());
		}

		fclose(fp_ids);

		Log("supporting sequence set was saved at " + out_file + "\n", true);

		// save fasta file
		this->CreateFastaFileUsingBlastdbcmd(tmp_file_ids, out_file + ".fasta");

		Log("supporting sequences were saved in FASTA format at " + out_file + ".fasta\n", true);
	}

	fclose(fp);
}

void SequenceDB::LoadSupportingSet(string in_file)
{
	FILE* fp = fopen(in_file.c_str(), "r");
	if (fp == NULL) {
		fprintf(stderr, "file open error: %s\n", in_file.c_str());
		exit(-1);
	}

	if (this->b_logging_) {
		Log("retrieving supporting sequences & cluster information...\n", true);
	}

	subject_seqs_.clear();
	max_seq_len_ = 0;

	char buf[BUF_SIZE_MED];
	char* id;
	char* strptr;
	int num_seqs = 0;
	int max_cluster_id = 0;

	while (!feof(fp)) {
		if (fgets(buf, BUF_SIZE_MED, fp) == NULL)
			break;

		if (buf[0] == '#') {
			if (strncmp(buf, "# PROVEAN version:\t", 19) == 0) {
				if (atof(&buf[19]) != atof(VERSION)) {
					fprintf(stderr, "error: incompatible supporting sequence set file (%s)\n", in_file.c_str());
					exit(-1);
				}
			} else if (strncmp(buf, "# No supporting sequences were found", 36) == 0) {
				// no supporting sequences were found in the database used
				// query sequences itself will be used as the only supporting sequence
				Log("the supporting set file indicates that no supporting sequences were found in the database used.\n", true);
				this->AddQueryAsSupportingSequence();
				return;
			}
			continue;
		}

		id = strtok(buf, "\t\n");

		if (id == NULL) {
			Log("Parsing Error\n", true);
			exit(-1);
		}

		Sequence seq;

		seq.id_ = id;

		strptr = strtok(NULL, "\t\n");
		seq.cluster_id_ = atoi(strptr);

		if (seq.cluster_id_ > max_cluster_id) {
			max_cluster_id = seq.cluster_id_;
		}

		strptr = strtok(NULL, "\t\n");
		seq.e_value_ = atof(strptr);

		strptr = strtok(NULL, "\t\n");
		seq.bit_score_ = atof(strptr);

		subject_seqs_.push_back(seq);
		if (max_seq_len_ < seq.seq_.length()) {
			max_seq_len_ = seq.seq_.length();
		}
		if (num_clusters_ < seq.cluster_id_) {
			num_clusters_ = seq.cluster_id_;
		}
		num_seqs++;
	}

	this->num_clusters_ = max_cluster_id;
	this->CountNumOfSeqsInEachCluster();

	fclose(fp);

	this->SetSequences();

}


int SequenceDB::CreateFastaFileOfSubjectSequences(string out_file)
{
	char tmp_file_id[BUF_SIZE_MED];
	sprintf(tmp_file_id, "%s/ids", this->temp_dir_.c_str());
	FILE* fp_id_out = fopen(tmp_file_id, "w");
	char buf_id[BUF_SIZE_MED];

	for (vector<Sequence>::iterator it = subject_seqs_.begin(); it < subject_seqs_.end(); it++) {
		strcpy(buf_id, (*it).id_.data());
		char* p1 = strtok(buf_id, "|");
		if (p1 == NULL) {
			Log("Parsing Error\n", true);
			exit(-1);
		}
		char* p2 = strtok(NULL, "\t\n");
		if (p2 == NULL) {
			Log("Parsing Error\n", true);
			exit(-1);
		}
		fprintf(fp_id_out, "%s|%s\n", p1, p2);
	}

	fclose(fp_id_out);

	this->CreateFastaFileUsingBlastdbcmdWithTargetOnly(tmp_file_id, out_file);

	return 0;
}

// fill sequences in subject_seqs_ vector
int SequenceDB::SetSequences()
{
	string tmp_file_fasta = this->temp_dir_ + "/fasta";

	this->CreateFastaFileOfSubjectSequences(tmp_file_fasta);
	this->SetOnlySequencesFromFastaFile(tmp_file_fasta);

	return 0;
}

int SequenceDB::SetOnlySequencesFromFastaFile(string fasta_file)
{
	if (fasta_file.empty())
		return -1;

	FILE* fp = fopen(fasta_file.c_str(), "r");

	if (fp == NULL) {
		cerr << "cannot open: " << fasta_file << endl;
		return -1;
	}

	if (b_logging_) {
		Log("loading subject sequences from a FASTA file...\n", true);
	}

	int num_seq = 0;
	char* strptr;
	char buf[BUF_SIZE_LARGE];

	Sequence seq;
	while (!feof(fp)) {
		if (fgets(buf, BUF_SIZE_LARGE, fp) == NULL)
			break;

		if (buf[0] == '>') {
			if (!seq.id_.empty() && !this->IsTooShortComparedToQuery(&seq)) {
				if (subject_seqs_[num_seq].id_.compare(seq.id_) == 0) {
					subject_seqs_[num_seq].seq_.assign(seq.seq_);
					num_seq++;
					if (max_seq_len_ < seq.seq_.length()) {
						max_seq_len_ = seq.seq_.length();
					}
				} else {
					fprintf(stderr, "IDs are not matched (supporting sequence set file:%s, fasta file:%s)\n", subject_seqs_[num_seq].id_.data(), seq.id_.data());
					exit(-1);
				}
			}
			// def line
			string def(buf, 100);
			size_t pos = def.find_first_of('\n', 0);
			if (pos != string::npos) {
				seq.def_ = def.substr(0, pos);
			} else {
				seq.def_ = def + " ...";
			}

			string id(buf, 10000);
			size_t pos1 = id.find_first_not_of(" ", 1);
			size_t pos2 = id.find_first_of(" \n", pos1);
			if (pos2 != string::npos) {
				seq.id_ = id.substr(pos1, pos2-pos1);
			} else {
				seq.id_ = id.substr(pos1, 10000-1-pos1);
			}

			while (strchr(buf, '\n') == NULL) {
				fgets(buf, sizeof(buf), fp);
			}

			//strptr = strtok(buf, " \n");
			//seq.id_.assign(++strptr);
			seq.seq_.clear();
		} else {
			strptr = strtok(buf, " \n");
			if (strptr == NULL)
				break;
			seq.seq_.append(strptr);
		}

		if (num_seq > 100000) {
			cerr << fasta_file << " has too many sequences. Please check the file." << endl;
			exit(-1);
		}
	}

//	if (num_seq != 0 && !this->IsTooShortComparedToQuery(&seq)) {
	if (num_seq != (int)this->subject_seqs_.size() && !this->IsTooShortComparedToQuery(&seq)) {
		if (subject_seqs_[num_seq].id_.compare(seq.id_) == 0) {
			subject_seqs_[num_seq].seq_ = seq.seq_;
			num_seq++;
			if (max_seq_len_ < seq.seq_.length()) {
				max_seq_len_ = seq.seq_.length();
			}
		} else {
			fprintf(stderr, "IDs are not matched (supporting sequence set file:%s, fasta file:%s)\n", subject_seqs_[num_seq].id_.data(), seq.id_.data());
			exit(-1);
		}
	}

	fclose(fp);

	return num_seq;
}


void SequenceDB::PrintDBStat(FILE* out)
{
	fprintf(out, "\n## Query and supporting sequences information ##\n");
	fprintf(out, "# Query filename:\t%s\n", query_seq_.id_.data());
	fprintf(out, "# Query protein length:\t%d\n", (int)query_seq_.seq_.length());
	fprintf(out, "# Number of clusters:\t%d\n", num_clusters_);
	fprintf(out, "# Number of supporting sequences:\t%d\n", (int)subject_seqs_.size());
	fprintf(out, "# Subject sequences\n");
	int count = 0;

	fprintf(out, "# Seq_num\tSequence_ID\tCluster_ID\tE-value\tBit_Score\tSequence_length\tDef_line\n");
	for (vector<Sequence>::iterator it = subject_seqs_.begin(); it < subject_seqs_.end(); it++) {
		fprintf(out, "%d:", ++count);
		fprintf(out, "\t%s", (*it).id_.c_str());
		fprintf(out, "\t%d", (*it).cluster_id_);
		fprintf(out, "\t%.0e", (*it).e_value_);
		fprintf(out, "\t%.1f", (*it).bit_score_);
		fprintf(out, "\t%d", (int)(*it).seq_.length());
		fprintf(out, "\t%s", (*it).def_.c_str());
		fprintf(out, "\n");
	}
	fprintf(out, "\n");
}

void SequenceDB::PrintNumOfSubjectSequences(FILE* out)
{
	fprintf(out, "# Number of clusters:\t%d\n", num_clusters_);
	fprintf(out, "# Number of supporting sequences used:\t%d\n", (int)subject_seqs_.size());
}


//
// return the number of clusters selected
//
int SequenceDB::ClusterAndSelectFromFastaFile(string fasta_file)
{
	char output_prefix[BUF_SIZE_MED];
	sprintf(output_prefix, "%s/cdhit", temp_dir_.c_str());
	char cdhit_command[BUF_SIZE_MED];
	int ret;

	bool b_need_bak_option = false;

	// CD-HIT version check
	//
	char help_filename[BUF_SIZE_MED];
	sprintf(help_filename, "%s.help", output_prefix);
	sprintf(cdhit_command, "%s -h > %s",
			p_options_->cdhit_command_.c_str(),
			help_filename);
	ret = system(cdhit_command);

	FILE *p_help = fopen(help_filename, "r");
	if (p_help == NULL) {
		fprintf(stderr, "cd-hit version checking failed\n");
		exit(-1);
	}
	char buf[BUF_SIZE_MED];
	while (!feof(p_help)) {
		if (fgets(buf, BUF_SIZE_MED, p_help) == NULL)
			break;

		char* p;
		if ((p = strstr(buf, "version")) != NULL) {
			p += strlen("version");
			p = strtok(p, " ");
			char* p_v1 = strtok(p, ".");
			char* p_v2 = strtok(NULL, ".");
			char* p_v3 = strtok(NULL, ".");

			int v1 = (p_v1 == NULL) ? 0:atoi(p_v1);
			int v2 = (p_v2 == NULL) ? 0:atoi(p_v2);
			int v3 = (p_v3 == NULL) ? 0:atoi(p_v3);

			if (v1 > 4 || (v1 == 4 && v2 > 5) || (v1 == 4 && v2 == 5 && v3 >= 8)) {
				b_need_bak_option = true;
			}
			break;
		}
	}
	fclose(p_help);

	// run CD-HIT
	//
	if (b_need_bak_option) {
		sprintf(cdhit_command, "%s -i %s -o %s.cluster -c %.2f -s %.1f -n 5 -l %d -bak 1 > %s.log",
			p_options_->cdhit_command_.c_str(),
			fasta_file.c_str(),
			output_prefix,
			p_options_->cluster_threshold_,
			(double)this->p_options_->percent_length_rep_/100.0,
			(int)(p_options_->percent_length_query_*this->query_seq_.seq_.length()/100),
			output_prefix);
	} else {
		sprintf(cdhit_command, "%s -i %s -o %s.cluster -c %.2f -s %.1f -n 5 -l %d > %s.log",
			p_options_->cdhit_command_.c_str(),
			fasta_file.c_str(),
			output_prefix,
			p_options_->cluster_threshold_,
			(double)this->p_options_->percent_length_rep_/100.0,
			(int)(p_options_->percent_length_query_*this->query_seq_.seq_.length()/100),
			output_prefix);
	}

//	if (this->b_logging_) {
		Log("clustering subject sequences...\n", true);
//	}


	ret = system(cdhit_command);
	if (ret != 0) {
		// check whether the total number of sequences is zero, which causes an error in cd-hit
		char log_filename[BUF_SIZE_MED];
		sprintf(log_filename, "%s.log", output_prefix);
		FILE* fp_log = fopen(log_filename, "r");
		char buf[BUF_SIZE_MED];
		while (!feof(fp_log)) {
			if (fgets(buf, BUF_SIZE_MED, fp_log) == NULL)
				break;
			if (strncmp(buf, "total seq:", 10) == 0) {
				if (atoi(&buf[10]) == 0) {
					subject_seqs_.clear();
					return 0;
				} else {
					break;
				}
			}
		}
		fprintf(stderr, "cd-hit failed (%d)\n", ret);
		exit(-1);
	}

	if (this->b_logging_) {
		Log("selecting clusters...\n", true);
	}

	// parse CD-HIT result
	//
	char cluster_filename[BUF_SIZE_MED];
	sprintf(cluster_filename, "%s.cluster.bak.clstr", output_prefix);

	FILE* fp_cluster = fopen(cluster_filename, "r");
	if (!fp_cluster) {
		Log("file open error!!\n", true);
		cerr << cluster_filename << endl;
		return -1;
	}

	int* cluster_id = new int[subject_seqs_.size()];
	memset(cluster_id, 0x00, subject_seqs_.size()*sizeof(int));

	int num_cluster = 0;
	char* strptr;
	int i = 0;
	int id;	// cdhit cluster id
	int compare_len;
	char* p_dots;

	while (!feof(fp_cluster)) {
		if (fgets(buf, BUF_SIZE_MED, fp_cluster) == NULL)
			break;

		id = atoi(buf);

		strptr = strstr(buf, ">");
		p_dots = strstr(strptr, "...");



		compare_len = ((int)subject_seqs_[i].id_.length() < p_dots-strptr-1) ? (subject_seqs_[i].id_.length()) : p_dots-strptr-1;

		while (subject_seqs_[i].id_.compare(0, compare_len, strptr+1, compare_len) != 0) {
			subject_seqs_[i].cluster_id_ = -1;
			i++;
			compare_len = ((int)subject_seqs_[i].id_.length() < p_dots-strptr-1) ? (subject_seqs_[i].id_.length()) : p_dots-strptr-1;
		}

		if (cluster_id[id] == 0) {
			num_cluster++;
			cluster_id[id] = num_cluster;
		}

		subject_seqs_[i].cluster_id_ = cluster_id[id];
		i++;
	}

	fclose(fp_cluster);

	delete[] cluster_id;

	if (num_cluster > p_options_->num_cluster_) {
		num_clusters_ = p_options_->num_cluster_;
	} else {
		num_clusters_ = num_cluster;
	}
	// remove clusters other than top N clusters
	for (vector<Sequence>::iterator it=subject_seqs_.begin(); it<subject_seqs_.end();) {
		if ((*it).cluster_id_ > num_clusters_ || (*it).cluster_id_ == -1) {
			it = subject_seqs_.erase(it);
		} else {
			it++;
		}
	}

	UpdateMaxSeqLength();

	this->CountNumOfSeqsInEachCluster();

	if (this->b_logging_) {
		char temp[BUF_SIZE_MED];
		sprintf(temp, "%d subject sequences in %d clusters were selected for supporting sequences.\n", (int)subject_seqs_.size(), num_clusters_);
		Log(temp, true);
	}

	return num_clusters_;
}

int SequenceDB::AddQueryAsSupportingSequence()
{
	Log("use the query itself as a supporting sequence\n", true);
	// add query sequence itself to subject sequence set
	int num_seqs = SetSequencesFromFastaFile(p_options_->query_file_name_);

	if (num_seqs == 1) {
		subject_seqs_[0].cluster_id_ = 1;
		num_clusters_ = 1;
		num_seqs_in_cluster_ = new int[2];
		num_seqs_in_cluster_[1] = 1;

		string q("[Query]");
		subject_seqs_[0].id_ = q.append(subject_seqs_[0].id_);
	} else {
		printf("Something wrong in query sequence file!\n");
		exit(-1);
	}

	this->b_supporting_itself = true;

	return 0;
}

void SequenceDB::CountNumOfSeqsInEachCluster()
{
	// count number of seqs in each cluster
	if (this->num_seqs_in_cluster_ != NULL) {
		delete[] num_seqs_in_cluster_;
	}
	num_seqs_in_cluster_ = new int[num_clusters_+1];
	for (int i=0; i<=num_clusters_; i++) {
		num_seqs_in_cluster_[i] = 0;
	}
	for (vector<Sequence>::iterator it=subject_seqs_.begin(); it<subject_seqs_.end(); it++) {
		num_seqs_in_cluster_[(*it).cluster_id_]++;
	}
}

int SequenceDB::UpdateMaxSeqLength()
{
	max_seq_len_ = 0;
	for (vector<Sequence>::iterator it=subject_seqs_.begin(); it<subject_seqs_.end(); it++) {
		if (max_seq_len_ < (*it).seq_.length()) {
			max_seq_len_ = (*it).seq_.length();
		}
	}

	return max_seq_len_;
}



void SequenceDB::ComputeDeltaScores()
{
//	if (this->b_logging_) {
		Log("computing delta alignment scores...\n", true);
//	}

	sort(variants_.begin(), variants_.end(), Variant::Compare);

	for (vector<Variant>::iterator it_var=variants_.begin(); it_var<variants_.end(); it_var++) {
		(*it_var).sum_weights_ = 0;
		(*it_var).delta_score_ = 0;
	}

	for (vector<Sequence>::iterator it_seq=subject_seqs_.begin(); it_seq<subject_seqs_.end(); it_seq++) {
		int ref_score;

		this->p_tables_->SetTable(query_seq_.seq_, (*it_seq).seq_, &(p_options_->score_matrix_));

		p_tables_->FillForwardTable();
		ref_score = p_tables_->FillBackwardTable();

		//ref_score = p_tables_->GetRefAlignmentScore();
		for (vector<Variant>::iterator it_var=variants_.begin(); it_var<variants_.end(); it_var++) {
			//ref_score = p_tables_->FillBackwardTable();
			int var_score = p_tables_->GetVarAlignmentScore((*it_var).pos_flank_left_, (*it_var).pos_flank_right_, (*it_var).str_added_);

			double weight;
			weight = 1.0/this->num_seqs_in_cluster_[(*it_seq).cluster_id_];

			(*it_var).delta_score_ += (double)(var_score-ref_score)*weight;
			(*it_var).sum_weights_ += weight;
		}
	}

	for (vector<Variant>::iterator it_var=variants_.begin(); it_var<variants_.end(); it_var++) {
		(*it_var).delta_score_ /= (*it_var).sum_weights_;
	}

	return;
}


void SequenceDB::PrintDeltaScores(FILE* out)
{
	if (b_logging_) {
		Log("printing PROVEAN scores...\n", true);
	}
	fprintf(out, "## PROVEAN scores ##\n");
	fprintf(out, "# VARIATION\tSCORE\n");
	sort(variants_.begin(), variants_.end(), Variant::CompareOrder);
	for (vector<Variant>::iterator it_var=variants_.begin(); it_var<variants_.end(); it_var++) {
		fprintf(out, "%s\t%.3f\n", (*it_var).input_.c_str(), (*it_var).delta_score_);
	}
}

void SequenceDB::PrintDeltaScores(FILE* out, string prefix)
{
	sort(variants_.begin(), variants_.end(), Variant::CompareOrder);
	for (vector<Variant>::iterator it_var=variants_.begin(); it_var<variants_.end(); it_var++) {
		fprintf(out, "[%s]\t%s\t%.3f\n", prefix.c_str(), (*it_var).input_.c_str(), (*it_var).delta_score_);
	}
}

void SequenceDB::AddAllSAPs()
{
	char AAs[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
	char hgvs[1024];

	for (unsigned int i=1; i<=this->query_seq_.seq_.length(); i++) {
		for (int k=0; k<(int)(sizeof(AAs)/sizeof(char)); k++) {
			sprintf(hgvs, "%c%d%c", query_seq_.seq_[i-1], i, AAs[k]);
			this->AddVariant_HGVS(hgvs);
		}
	}
}

void SequenceDB::PrintAllSAPs()
{
	char AAs[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
	int aa2index[] = {
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA, 0,NA, 1, 2, 3, 4, 5, 6, 7,NA, 8, 9,10,11,NA,
		12,13,14,15,16,NA,17,18,NA,19,NA,NA,NA,NA,NA,NA,
		NA, 0,NA, 1, 2, 3, 4, 5, 6, 7,NA, 8, 9,10,11,NA,
		12,13,14,15,16,NA,17,18,NA,19,NA,NA,NA,NA,NA,NA,

		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
		NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA };

	double** scores = new double*[query_seq_.seq_.length()];
	for (unsigned int i=0; i<query_seq_.seq_.length(); i++) {
		scores[i] = new double[20];
		for (int j=0; j<20; j++) {
			scores[i][j] = -1000;
		}
	}

	for (vector<Variant>::iterator it_var=variants_.begin(); it_var<variants_.end(); it_var++) {
		if ((*it_var).input_.find("del") == string::npos && (*it_var).input_.find("ins") == string::npos && (*it_var).input_.find("dup") == string::npos) {
			char var_aa;
			//char ref_aa, var_aa;
			//ref_aa = (*it_var).input_[0];
			int var_pos = (*it_var).input_.find_first_not_of("0123456789", 1);
			var_aa = (*it_var).input_[var_pos];
			int pos = atoi((*it_var).input_.substr(1, var_pos-1).c_str());

			scores[pos-1][aa2index[(int)var_aa]] = (*it_var).delta_score_;
		}
	}

	// heading
	printf("# ENST AAPOS AA");
	for (int k=0; k<(int)(sizeof(AAs)/sizeof(char)); k++) {
		printf(" %c", AAs[k]);
	}
	printf("\n");

	for (unsigned int i=1; i<=this->query_seq_.seq_.length(); i++) {
		printf("[allsap score]\t%s %d %c", query_seq_.id_.c_str(), i, query_seq_.seq_[i-1]);
		for (int k=0; k<20; k++) {
			printf(" %.3f", scores[i-1][k]);
		}
		printf("\n");
	}
}

int SequenceDB::AddVariant_HGVS(string hgvs_in)
{
	Variant v;

	size_t pos1 = hgvs_in.find_first_not_of(" ");
	size_t pos2 = hgvs_in.find_first_of(" \t\f\v\r\n", pos1);
	if (pos2 != string::npos) {
		v.input_ = hgvs_in.substr(pos1, pos2-pos1);
	} else {
		v.input_ = hgvs_in;
	}

	string hgvs = v.input_;

	if (hgvs.find("delins") != string::npos) {
		// indels
		char del_start_aa = hgvs[0];
		int idx = hgvs.find_first_not_of("0123456789", 1);
		int del_start_pos = atoi(hgvs.substr(1, idx-1).c_str());
		int del_end_pos;
		char del_end_aa;
		if (hgvs[idx] == '_') {
			// more than 1 AA deletion
			del_end_aa = hgvs[idx+1];
			int idx2 = hgvs.find_first_not_of("0123456789", idx+2);
			del_end_pos = atoi(hgvs.substr(idx+2, idx2-idx-2).c_str());
		} else if (hgvs[idx] == 'd') {
			if (hgvs[idx+1] == 'e' && hgvs[idx+2] == 'l') {
				// 1 AA deletion
				del_end_pos = del_start_pos;
				del_end_aa = del_start_aa;
			} else {
				fprintf(stderr, "%s: invalid format\n", hgvs.c_str());
				return -1;
			}
		} else {
			// invalid format
			fprintf(stderr, "%s: invalid format\n", hgvs.c_str());
			return -1;
		}

		string inserted_aas = hgvs.substr(hgvs.find("delins")+6);

		if (inserted_aas.empty()) {
			fprintf(stderr, "%s: no inserted amino acids provided\n", hgvs.c_str());
			return -1;
		}

		if (inserted_aas.find_first_not_of("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy") != string::npos) {
			// invalid AAs
			fprintf(stderr, "%s: invalid AAs provided for insertion\n", hgvs.c_str());
			return -1;
		}

		// checking AA
		if ((int)query_seq_.seq_.length() < del_start_pos || (int)query_seq_.seq_.length() < del_end_pos || del_start_pos > del_end_pos) {
			fprintf(stderr, "%s: invalid position, (%d,%d) (length of query sequence: %d)\n", hgvs.c_str(), del_start_pos, del_end_pos, (int)query_seq_.seq_.length());
			return -1;
		}
		if (query_seq_.seq_[del_start_pos-1] != del_start_aa) {
			fprintf(stderr, "%s: reference AA does not match: %c in query sequence, but %c is provided\n", hgvs.c_str(), query_seq_.seq_[del_start_pos-1], del_start_aa);
			return -1;
		}
		if (query_seq_.seq_[del_end_pos-1] != del_end_aa) {
			fprintf(stderr, "%s: reference AA does not match: %c in query sequence, but %c is provided\n", hgvs.c_str(), query_seq_.seq_[del_end_pos-1], del_end_aa);
			return -1;
		}
		// END checking AA

		v.pos_flank_left_ = del_start_pos -1;
		v.pos_flank_right_ = del_end_pos +1;
		v.str_added_ = inserted_aas;

	} else if (hgvs.find("del") != string::npos) {
		// deletions
		char del_start_aa = hgvs[0];
		int idx = hgvs.find_first_not_of("0123456789", 1);
		int del_start_pos = atoi(hgvs.substr(1, idx-1).c_str());
		int del_end_pos;
		char del_end_aa;
		if (hgvs[idx] == '_') {
			// more than 1 AA deletion
			del_end_aa = hgvs[idx+1];
			int idx2 = hgvs.find_first_not_of("0123456789", idx+2);
			del_end_pos = atoi(hgvs.substr(idx+2, idx2-idx-2).c_str());
		} else if (hgvs[idx] == 'd'){
			if (hgvs[idx+1] == 'e' && hgvs[idx+2] == 'l') {
				// 1 AA deletion
				del_end_pos = del_start_pos;
				del_end_aa = del_start_aa;
			} else {
				fprintf(stderr, "%s: invalid format\n", hgvs.c_str());
				return -1;
			}
		} else {
			// invalid format
			fprintf(stderr, "%s: invalid format\n", hgvs.c_str());
			return -1;
		}

		// checking AA
		if ((int)query_seq_.seq_.length() < del_start_pos || (int)query_seq_.seq_.length() < del_end_pos || del_start_pos > del_end_pos) {
			fprintf(stderr, "%s: invalid position:,(%d,%d) (length of query sequence: %d)\n", hgvs.c_str(), del_start_pos, del_end_pos, (int)query_seq_.seq_.length());
			return -1;
		}
		if (query_seq_.seq_[del_start_pos-1] != del_start_aa) {
			fprintf(stderr, "%s: reference AA does not match: %c in query sequence, but %c is provided\n", hgvs.c_str(), query_seq_.seq_[del_start_pos-1], del_start_aa);
			return -1;
		}
		if (query_seq_.seq_[del_end_pos-1] != del_end_aa) {
			fprintf(stderr, "%s: reference AA does not match: %c in query sequence, but %c is provided\n", hgvs.c_str(), query_seq_.seq_[del_end_pos-1], del_end_aa);
			return -1;
		}
		// END checking AA

		v.pos_flank_left_ = del_start_pos -1;
		v.pos_flank_right_ = del_end_pos +1;
		v.str_added_ = "";

	} else if (hgvs.find("ins") != string::npos) {
		// insertions
		char ins_start_aa = hgvs[0];
		int idx = hgvs.find_first_not_of("0123456789", 1);
		int ins_start_pos = atoi(hgvs.substr(1, idx-1).c_str());
		if (hgvs[idx] != '_' && hgvs[idx] != 'i') {
			cerr << "Invalid format: " << hgvs << endl;
			return -1;
		}

		char ins_end_aa;
		int ins_end_pos;
		if (hgvs[idx] == '_') {
			ins_end_aa = hgvs[idx+1];
			int idx2 = hgvs.find_first_not_of("0123456789", idx+2);
			ins_end_pos = atoi(hgvs.substr(idx+2, idx2-idx-2).c_str());
			if (ins_start_pos+1 != ins_end_pos) {
				cerr << "Invalid format: " << hgvs << endl;
				return -1;
			}
		} else {
			ins_end_aa = ins_start_aa;
			ins_end_pos = ins_start_pos;
		}

		string inserted_aas = hgvs.substr(hgvs.find("ins")+3);
		if (inserted_aas.empty()) {
			fprintf(stderr, "%s: no inserted amino acids provided\n", hgvs.c_str());
			return -1;
		}

		if (inserted_aas.find_first_not_of("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy") != string::npos) {
			// invalid AAs
			fprintf(stderr, "%s: invalid AAs provided for insertion\n", hgvs.c_str());
			return -1;
		}

		// checking AA
		if ((int)query_seq_.seq_.length() < ins_start_pos || (int)query_seq_.seq_.length() < ins_end_pos) {
			fprintf(stderr, "Invalid position: (%d,%d) (length of query sequence: %d)\n", ins_start_pos, ins_end_pos, (int)query_seq_.seq_.length());
			return -1;
		}
		if (query_seq_.seq_[ins_start_pos-1] != ins_start_aa) {
			fprintf(stderr, "Reference AA does not match: %c in Sequence, but %c is provided\n", query_seq_.seq_[ins_start_pos-1], ins_start_aa);
			return -1;
		}
		if (query_seq_.seq_[ins_end_pos-1] != ins_end_aa) {
			fprintf(stderr, "Reference AA does not match: %c in Sequence, but %c is provided\n", query_seq_.seq_[ins_end_pos-1], ins_end_aa);
			return -1;
		}
		// END checking AA

		v.pos_flank_left_ = ins_start_pos;
		v.pos_flank_right_ = ins_end_pos;
		v.str_added_ = inserted_aas;

	} else if (hgvs.find("dup") != string::npos) {
		// duplications (same as insertions)
		char ins_start_aa = hgvs[0];
		int idx = hgvs.find_first_not_of("0123456789", 1);
		int ins_start_pos = atoi(hgvs.substr(1, idx-1).c_str());
		char ins_end_aa;
		int ins_end_pos;
		if (hgvs[idx] == 'd') {
			// 1 AA insertion
			ins_end_aa = ins_start_aa;
			ins_end_pos = ins_start_pos;
		} else if (hgvs[idx] == '_') {
			// more than 1 AAs insertion
			ins_end_aa = hgvs[idx+1];
			int idx2 = hgvs.find_first_not_of("0123456789", idx+2);
			ins_end_pos = atoi(hgvs.substr(idx+2, idx2-idx-2).c_str());
			if (ins_start_pos > ins_end_pos) {
				cerr << "Invalid format: " << hgvs << endl;
				return -1;
			}
		} else {
			cerr << "Invalid format: " << hgvs << endl;
			return -1;
		}

		// checking AA
		if ((int)query_seq_.seq_.length() < ins_start_pos || (int)query_seq_.seq_.length() < ins_end_pos) {
			fprintf(stderr, "Invalid position: (%d,%d) (length of query sequence: %d)\n", ins_start_pos, ins_end_pos, (int)query_seq_.seq_.length());
			return -1;
		}
		if (query_seq_.seq_[ins_start_pos-1] != ins_start_aa) {
			fprintf(stderr, "Reference AA does not match: %c in Sequence, but %c is provided\n", query_seq_.seq_[ins_start_pos-1], ins_start_aa);
			return -1;
		}
		if (query_seq_.seq_[ins_end_pos-1] != ins_end_aa) {
			fprintf(stderr, "Reference AA does not match: %c in Sequence, but %c is provided\n", query_seq_.seq_[ins_end_pos-1], ins_end_aa);
			return -1;
		}
		// END checking AA

		string inserted_aas = query_seq_.seq_.substr(ins_start_pos-1, ins_end_pos-ins_start_pos+1);

		ins_start_pos = ins_end_pos;
		ins_end_pos = ins_end_pos+1;

		v.pos_flank_left_ = ins_start_pos;
		v.pos_flank_right_ = ins_end_pos;
		v.str_added_ = inserted_aas;

	} else {
		// substitutions
		char ref_aa;

		ref_aa = hgvs[0];
		size_t var_pos = hgvs.find_first_not_of("0123456789", 1);
		if (var_pos == string::npos) {
			fprintf(stderr, "%s: invalid format\n", hgvs.c_str());
			return -1;
		}
		int position = atoi(hgvs.substr(1, var_pos-1).c_str());

		string var_aas = hgvs.substr(var_pos);

		if (var_aas.length() != 1) {
			fprintf(stderr, "%s: invalid format\n", hgvs.c_str());
			return -1;
		}

		if (var_aas.find_first_not_of("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy") != string::npos) {
			// invalid AAs
			fprintf(stderr, "%s: invalid AAs provided\n", hgvs.c_str());
			return -1;
		}

		// checking AA
		if ((int)query_seq_.seq_.length() < position) {
			fprintf(stderr, "Invalid position: %d (length of query sequence: %d)\n", position, (int)query_seq_.seq_.length());
			return -1;
		}
		if (query_seq_.seq_[position-1] != ref_aa) {
			fprintf(stderr, "Reference AA does not match: %c in Sequence, but %c is provided\n", query_seq_.seq_[position-1], ref_aa);
			return -1;
		}
		// END checking AA

		v.pos_flank_left_ = position -1;
		v.pos_flank_right_ = position +1;
		v.str_added_ = var_aas;

	}

	this->num_variants_++;
	v.order_ = this->num_variants_;
	this->variants_.push_back(v);

	return 0;
}


// add variant in comma-separated (or space-separated) value format
// <position>,<reference AAs>,<variant AAs>
// <position> <reference AAs> <variant AAs>
int SequenceDB::AddVariant_CSV(string csv_in)
{
	Variant v;
	size_t pos = csv_in.find_first_of("\t\f\v\r\n");
	if (pos != string::npos) {
		v.input_ = csv_in.substr(0, pos);
	} else {
		v.input_ = csv_in;
	}

	char* var_str = new char [csv_in.length()+1];
	strcpy(var_str, csv_in.c_str());

	char* strptr;

	if ((strptr = strtok(var_str, " ,")) == NULL) {
		fprintf(stderr, "%s: invalid format\n", v.input_.c_str());
		delete [] var_str;
		return -1;
	}

	int AA_pos = atoi(strptr);

	if (AA_pos == 0) {
		fprintf(stderr, "%s: invalid format\n", v.input_.c_str());
		delete [] var_str;
		return -1;
	}

	// checking position
	if ((int)query_seq_.seq_.length() < AA_pos) {
		fprintf(stderr, "%s: invalid position, %d (length of query sequence: %d)\n", v.input_.c_str(), AA_pos, (int)query_seq_.seq_.length());
		delete [] var_str;
		return -1;
	}


	if ((strptr = strtok(NULL, " ,")) == NULL) {
		fprintf(stderr, "%s: invalid format\n", v.input_.c_str());
		delete [] var_str;
		return -1;
	}

	char* ref_AA = strptr;

	// will add checking if reference AAs are from 20 AAs
	//


	if ((strptr = strtok(NULL, " ,\r\t\n")) == NULL) {
		fprintf(stderr, "%s: invalid format\n", v.input_.c_str());
		delete [] var_str;
		return -1;
	}

	char* var_AA = strptr;

	// will add checking if variant AAs are from 20 AAs
	//

	if (strchr(var_AA,'.') != NULL && strlen(var_AA) > 1) {
		fprintf(stderr, "%s: invalid format\n", v.input_.c_str());
		delete [] var_str;
		return -1;
	}

	// checking reference AAs
	if (query_seq_.seq_.compare(AA_pos-1, strlen(ref_AA), ref_AA, strlen(ref_AA)) != 0) {
		fprintf(stderr, "%s: reference AAs do not match: %s in query sequence, but %s is provided\n", v.input_.c_str(), query_seq_.seq_.substr(AA_pos-1, strlen(ref_AA)).c_str(), ref_AA);
		delete [] var_str;
		return -1;
	}

	if (var_AA[0] == '.') {
		var_AA[0] = '\0';
	}

	while (ref_AA[0] && var_AA[0] && (ref_AA[0] == var_AA[0])) {
		AA_pos++;
		ref_AA++;
		var_AA++;
	}

	while (strlen(ref_AA)>0 && strlen(var_AA)>0 && (ref_AA[strlen(ref_AA)-1] == var_AA[strlen(var_AA)-1])) {
		ref_AA[strlen(ref_AA)-1] = '\0';
		var_AA[strlen(var_AA)-1] = '\0';
	}


	v.pos_flank_left_ = AA_pos-1;
	v.pos_flank_right_ = AA_pos+strlen(ref_AA);
	v.str_added_ = var_AA;

	if (v.str_added_.find_first_not_of("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy") != string::npos) {
		// invalid AAs
		fprintf(stderr, "%s: invalid AAs provided\n", v.input_.c_str());
		return -1;
	}

	this->num_variants_++;
	v.order_ = this->num_variants_;
	this->variants_.push_back(v);

	delete [] var_str;

	return 0;
}

int SequenceDB::AddVariantsFromFile(string filename)
{
	if (!p_options_->flag_quiet_) {
		Log("loading variations...\n", true);
	}

	FILE* fp = fopen(filename.c_str(), "r");

	if (fp == NULL) {
		cerr << "cannot open: " << filename << endl;
		return -1;
	}

	char buf[BUF_SIZE_MED];

	int count = 0;

	while (!feof(fp)) {
		if (fgets(buf, BUF_SIZE_MED, fp) == NULL)
			break;

		if (strchr(buf, ',') != NULL) {
			this->AddVariant_CSV(buf);
		} else {
			int i = 0;
			while (buf[i]) {
				if (! isspace(buf[i])) {
					this->AddVariant_HGVS(buf);
					break;
				}
				i++;
			}
		}

		count++;
	}

	fclose(fp);

	return count;
}

int SequenceDB::CreateTempDir()
{
	char temp_dir[BUF_SIZE_MED];
	if (p_options_->tmp_dir_given_ == 1) {
		sprintf(temp_dir, "%s/proveanXXXXXX", p_options_->tmp_dir_.data());
	} else {
		sprintf(temp_dir, "%s/proveanXXXXXX", P_tmpdir);
	}

	if (mkdtemp(temp_dir) == NULL) {
		fprintf(stderr, "Error in creating temporary directory\n");
		return -1;
	}

	temp_dir_ = temp_dir;

	return 0;
}

void SequenceDB::RunBlast(string query_file, string blast_db_file, string blast_out)
{
	int max_num = 100000;
	char blast_command[BUF_SIZE_MED];
	int return_code;
	char blast_help[BUF_SIZE_MED];

	sprintf(blast_help, "%s/psiblast.help", temp_dir_.c_str());

	// version check
	sprintf(blast_command, "%s -version > %s", p_options_->psiblast_command_.c_str(), blast_help);

	return_code = system(blast_command);

	if (return_code != 0) {
		fprintf(stderr, "\"%s\" failed (exit:%d)\n", blast_command, return_code);
		exit(-1);
	}

	bool b_use_max_target = false;
	FILE *p_help = fopen(blast_help, "r");
	if (p_help == NULL) {
		fprintf(stderr, "psiblast version checking failed\n");
		exit(-1);
	}
	char buf[BUF_SIZE_MED];
	while (!feof(p_help)) {
		if (fgets(buf, BUF_SIZE_MED, p_help) == NULL)
			break;

		char* p;
		if ((p = strstr(buf, "psiblast:")) != NULL) {
			p += strlen("psiblast:");
			p = strtok(p, " ");
			char* p_v1 = strtok(p, ".");
			char* p_v2 = strtok(NULL, ".");
			char* p_v3 = strtok(NULL, ".");

			int v1 = (p_v1 == NULL) ? 0:atoi(p_v1);
			int v2 = (p_v2 == NULL) ? 0:atoi(p_v2);
			int v3 = (p_v3 == NULL) ? 0:atoi(p_v3);

			if (v1 > 2 || (v1 == 2 && v2 > 2) || (v1 == 2 && v2 == 2 && v3 >= 26)) {
				b_use_max_target = true;
			}
			break;
		}
	}
	fclose(p_help);

	if (b_use_max_target) {
		sprintf(blast_command, "%s -query %s -db %s -show_gis -max_target_seqs %d -num_iterations %d -outfmt \"6 sseqid sgi evalue bitscore\" -evalue %f -out %s -num_threads %d",
				p_options_->psiblast_command_.c_str(),
				query_file.c_str(),
				blast_db_file.c_str(),
				max_num,
				1,
				0.1,
				blast_out.c_str(),
				p_options_->num_threads_
				);
	} else {
		sprintf(blast_command, "%s -query %s -db %s -show_gis -num_alignments %d -num_descriptions %d -num_iterations %d -outfmt \"6 sseqid sgi evalue bitscore\" -evalue %f -out %s -num_threads %d",
				p_options_->psiblast_command_.c_str(),
				query_file.c_str(),
				blast_db_file.c_str(),
				max_num,
				max_num,
				1,
				0.1,
				blast_out.c_str(),
				p_options_->num_threads_
				);

	}
	Log("searching related sequences...\n", true);

//	if (this->b_logging_) {
//		Log(blast_command, true);
//	}

	return_code = system(blast_command);

	if (return_code != 0) {
		fprintf(stderr, "psiblast failed (exit:%d)\n", return_code);
		exit(-1);
	}

}

int SequenceDB::CreateFastaFileFromBlastOut(string blast_file, string fasta_file)
{
	if (blast_file.empty())
		return -1;

	char tmp_file_id[BUF_SIZE_MED];
	FILE* fp = fopen(blast_file.c_str(), "r");

	if (fp == NULL) {
		cerr << "cannot open: " << blast_file << endl;
		return -1;
	}

	char buf[BUF_SIZE_MED];

	char* db_id;
	char* gi;
	int num_seqs = 0;

	sprintf(tmp_file_id, "%s/ids", this->temp_dir_.c_str());
	FILE* fp_id_out = fopen(tmp_file_id, "w");

	char last_id[BUF_SIZE_MED] = "";

	// get subject seqs
	while (!feof(fp)) {
		if (fgets(buf, BUF_SIZE_MED, fp) == NULL)
			break;

		db_id = strtok(buf, "\t\n");
		gi = strtok(NULL, "\t\n");
		char id[BUF_SIZE_MED]="";

		if (db_id == NULL || gi == NULL) {
			Log("Parsing Error\n", true);
			return -1;
		}

		strcpy(id, db_id);

		if (db_id != NULL && gi != NULL) {
			if (strcmp(last_id, id) == 0) {
				continue;
			}

			fprintf(fp_id_out, "%s\n", id);
			strcpy(last_id, id);
			num_seqs++;
		} else {
			Log("Wrong entry in blast out\n", true);
			break;
		}
	}

	fclose(fp);
	fclose(fp_id_out);

	this->CreateFastaFileUsingBlastdbcmd(tmp_file_id, fasta_file);

	return num_seqs;
}

int SequenceDB::CreateFastaFileUsingBlastdbcmd(string id_file, string fasta_file)
{
	char blastdbcmd[BUF_SIZE_MED];
	sprintf(blastdbcmd, "%s -db %s -dbtype 'prot' -entry_batch %s > %s",
			this->p_options_->blastdbcmd_command_.c_str(),
			this->p_options_->blast_db_file_name_.c_str(),
			id_file.c_str(),
			fasta_file.c_str()
			);
	int return_code = system(blastdbcmd);

	if (return_code != 0) {
		fprintf(stderr, "blastdbcmd failed (exit:%d)\n", return_code);
		exit(-1);
	}

	return 0;
}

int SequenceDB::CreateFastaFileUsingBlastdbcmdWithTargetOnly(string id_file, string fasta_file)
{
	char blastdbcmd[BUF_SIZE_MED];
	sprintf(blastdbcmd, "%s -db %s -dbtype 'prot' -entry_batch %s -target_only > %s",
			this->p_options_->blastdbcmd_command_.c_str(),
			this->p_options_->blast_db_file_name_.c_str(),
			id_file.c_str(),
			fasta_file.c_str()
			);
	int return_code = system(blastdbcmd);

	if (return_code != 0) {
		fprintf(stderr, "blastdbcmd failed (exit:%d)\n", return_code);
		exit(-1);
	}

	return 0;
}
