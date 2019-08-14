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
// Name        : SequenceDB.h
// Author      : Yongwook Choi
//============================================================================

#ifndef SEQUENCEDB_H_
#define SEQUENCEDB_H_

#include <vector>
#include "Common.h"
#include "Sequence.h"
#include "ScoreMatrix.h"
#include "Variant.h"
#include "Tables.h"
#include "Options.h"

using namespace std;

class SequenceDB {
private:
	void RunBlast(string query_file, string blast_db_file, string blast_out);
	bool b_logging_;

	unsigned int max_seq_len_;
	Sequence query_seq_;
	vector<Sequence> subject_seqs_;
	vector<Variant> variants_;
	Tables* p_tables_;
	int num_clusters_;
	int* num_seqs_in_cluster_;
	int num_variants_;
	string temp_dir_;
	Options* p_options_;

	bool b_supporting_itself;

	void CountNumOfSeqsInEachCluster();

	bool IsTooShortComparedToQuery(Sequence* seq);
	int SetOnlySequencesFromFastaFile(string fasta_file);

	int SetSequences();
	int SetSequenceInfoFromBlastOut(string blast_file, int max_num_seqs);
	int CreateFastaFileUsingBlastdbcmd(string id_file, string fasta_file);
	int CreateFastaFileUsingBlastdbcmdWithTargetOnly(string id_file, string fasta_file);
	int CreateFastaFileOfSubjectSequences(string out_file);
	int AddQueryAsSupportingSequence();
	void SaveSupportingSet(string out_file);
	void LoadSupportingSet(string in_file);
	int ClusterAndSelectFromFastaFile(string fasta_file);
	int UpdateMaxSeqLength();

	int SetSequencesFromBlastOut(string blast_file, string blast_db_file, int max_num_seqs);
	int SetSequencesFromFastaFile(string fasta_file);
	int AddVariant_HGVS(string hgvs);
	int AddVariant_CSV(string csv_in);


public:
	SequenceDB();
	virtual ~SequenceDB();

	void SetOptions(Options* p_options);
	bool SetQuerySequenceFromFastaFile(const char* id, string fasta_file);

	int CreateFastaFileFromBlastOut(string blast_file, string fasta_file);
	void SetSubjectSequences();
	void SetLoggingOption(bool b);
	int CreateTempDir();

	int GetNumberOfVariants();
	void ComputeDeltaScores();
	void AddAllSAPs();
	int AddVariantsFromFile(string filename);
	void PrintAllSAPs();
	void PrintNumOfSubjectSequences(FILE* out);
	void PrintDeltaScores(FILE* out);
	void PrintDBStat(FILE* out);
	void PrintDeltaScores(FILE* out, string prefix);

};

#endif /* SEQUENCEDB_H_ */
