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
// Name        : Common.cpp
// Author      : Yongwook Choi
//============================================================================

#include <time.h>
#include "Common.h"

int Max3(int a, int b, int c)
{
	if (a > b) {
		if (a > c)
			return a;
		else
			return c;
	} else {
		if (b > c)
			return b;
		else
			return c;
	}
}

void Log(FILE* out, string message, bool with_time)
{
	if (with_time) {
		time_t rawtime;
		struct tm* timeinfo;

		time(&rawtime);
		timeinfo = localtime(&rawtime);

		fprintf(out, "[%02d:%02d:%02d] %s", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, message.data());
	} else {
		fprintf(out, "%s", message.data());
	}
	fflush(out);
}

void Log(string message, bool with_time)
{
	Log(stdout, message, with_time);
}
