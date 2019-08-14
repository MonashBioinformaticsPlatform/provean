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
// Name        : provean.cpp
// Author      : Yongwook Choi
//============================================================================

#include <string>
#include <iostream>
#include <algorithm>

#include "Common.h"
#include "Options.h"
#include "Sequence.h"
#include "SequenceDB.h"

using namespace std;

int main(int argc, char** argv)
{

	Options options;
	options.SetOptions(argc, argv);

	fprintf(stdout, "## PROVEAN v1.1 output ##\n");

	if (!options.flag_quiet_) {
		fprintf(stdout, "## Input Arguments ##\n");
		for (int i=0; i<argc; i++) {
			fprintf(stdout, "%s ", argv[i]);
		}
		fprintf(stdout, "\n");
	}


	options.PrintOptions(stdout);

	SequenceDB db;

	if (options.flag_quiet_) {
		db.SetLoggingOption(false);
	} else {
		db.SetLoggingOption(true);
	}

	db.SetOptions(&options);

	if (db.CreateTempDir() != 0) {
		exit(-1);
	}

	db.SetQuerySequenceFromFastaFile(options.query_file_name_.c_str(), options.query_file_name_);

	if (options.variation_.compare("") != 0) {
		db.AddVariantsFromFile(options.variation_);
	}

	if (options.flag_all_sap_) {
		if (!options.flag_quiet_) {
			Log("adding all single amino acid polymorphisms as variations...\n", true);
		}
		db.AddAllSAPs();
	}

	if (db.GetNumberOfVariants() == 0) {
		fprintf(stdout, "No variations entered\n");
		return 0;
	}


	db.SetSubjectSequences();

	db.PrintNumOfSubjectSequences(stdout);

	db.ComputeDeltaScores();

	if (options.flag_all_sap_) {
		db.PrintAllSAPs();
	} else {
		db.PrintDeltaScores(stdout);
	}

	return 0;
}

