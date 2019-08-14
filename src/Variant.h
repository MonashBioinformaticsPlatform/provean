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
// Name        : Variant.h
// Author      : Yongwook Choi
//============================================================================

#ifndef VARIANT_H_
#define VARIANT_H_

#include <string>
#include "Common.h"

using namespace std;

class Variant {
public:
	Variant();
	virtual ~Variant();

	string input_;
	int pos_flank_left_;
	int pos_flank_right_;
	string str_added_;
	double sum_weights_;
	int order_;
	double delta_score_;

	static bool Compare(Variant v1, Variant v2);
	static bool CompareOrder(Variant v1, Variant v2);

};

#endif /* VARIANT_H_ */
