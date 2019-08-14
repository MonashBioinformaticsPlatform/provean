bin_SCRIPTS = provean.sh
CLEANFILES = provean.sh
EXTRA_DIST = provean.sh.in

do_subst = sed -e 's,\[@\]PSIBLAST\[@\],$(PSIBLAST),g' \
				-e 's,\[@\]CDHIT\[@\],$(CDHIT),g' \
				-e 's,\[@\]BLASTDB\[@\],$(BLAST_DB),g' \
				-e 's,\[@\]BLASTDBCMD\[@\],$(BLASTDBCMD),g' \
				-e 's,\[@\]VERSION\[@\],$(VERSION),g'

provean.sh : provean.sh.in Makefile
	$(do_subst) < $(srcdir)/provean.sh.in > provean.sh
	chmod +x provean.sh


