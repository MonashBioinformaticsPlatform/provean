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
// Name        : Common.h
// Author      : Yongwook Choi
//============================================================================

#ifndef COMMON_H_
#define COMMON_H_

#define BUF_SIZE_LARGE 0x100000
#define BUF_SIZE_MED 0x10000
#define BUF_SIZE_SMALL 0x1000

#define VARIATION_SUB 0
#define VARIATION_DEL 1
#define VARIATION_INS 2
#define VARIATION_INDEL 3

#define VERSION "1.1"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <cstdlib>


using namespace std;

int Max3(int a, int b, int c);

void Log(FILE* out, string message, bool with_time);
void Log(string message, bool with_time);

#endif /* COMMON_H_ */
