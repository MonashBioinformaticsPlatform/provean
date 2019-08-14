AC_INIT([provean], [1.1.5], [ychoi@jcvi.org])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Makefile
src/Makefile
scripts/Makefile
examples/Makefile
])

# check psiblast
AC_PATH_PROG([PSIBLAST], [psiblast], [no])
if test "$PSIBLAST" = "no"
then
	AC_MSG_ERROR(NCBI psiblast missing from path. Specify full path. (e.g. ./configure PSIBLAST=/path/to/psiblast/exe))
fi

# check cdhit
AC_PATH_PROG([CDHIT], [cd-hit], [no])
if test "$CDHIT" = "no"
then
	AC_MSG_ERROR(cd-hit missing from path. Specify full path. (e.g. ./configure CDHIT=/path/to/cd-hit/exe))
fi

# check blastdbcmd
AC_PATH_PROG([BLASTDBCMD], [blastdbcmd], [no])
if test "$BLASTDBCMD" = "no"
then
	AC_MSG_ERROR(blastdbcmd missing from path. Specify full path. (e.g. ./configure BLASTDBCMD=/path/to/blastdbcmd/exe))
fi

AC_SUBST([BLAST_DB])

AC_OUTPUT

