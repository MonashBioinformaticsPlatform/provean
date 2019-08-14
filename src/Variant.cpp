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
// Name        : Variant.cpp
// Author      : Yongwook Choi
//============================================================================

#include "Variant.h"

Variant::Variant() {
	delta_score_ = 0.0;
	sum_weights_ = 0.0;
}

Variant::~Variant() {
	// TODO Auto-generated destructor stub
}

bool Variant::Compare(Variant v1, Variant v2)
{
	if (v1.pos_flank_left_ == v2.pos_flank_left_) {
		return (v1.pos_flank_right_ < v2.pos_flank_right_);
	} else {
		return (v1.pos_flank_left_ < v2.pos_flank_left_);
	}
}

bool Variant::CompareOrder(Variant v1, Variant v2)
{
	return (v1.order_ < v2.order_);
}
