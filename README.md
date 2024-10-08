##################
# IMPORTANT

> This is just a fork of the original to allow us to fix a issue we found.
> We cannot provide support for this tool.
> -- MonashBioinformaticsPlatform

##################

```
PROVEAN v.1.1.5 (May 7, 2014)


0. QUICK INSTALLATION

	$ tar zxvf provean-1.1.5.tar.gz
	$ cd provean-1.1.5
	$ ./configure
	$ make
	$ make install


I. PREREQUISITES

	PROVEAN requires the following software and database.

  1. NCBI BLAST 2.2.28+ (or more recent)

	This is available at the NCBI ftp site.
	ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

  2. CD-HIT 3.1.2 (or more recent, but currently v4.6 and v4.6.1 are not
  	recommended since those versions have a reported problem,
	https://code.google.com/p/cdhit/issues/detail?id=18)

	This is available at the CD-HIT website.
	http://weizhong-lab.ucsd.edu/cd-hit/download.php
 
  3. NCBI nr (non-redundant) protein database 

	Only BLAST pre-formatted databases are needed.The NCBI nr databases
	(released August 2011) used for our publication	is available at the
	JCVI ftp site. Download and unpack nr.*.tar.gz files.
	ftp://ftp.jcvi.org/pub/data/provean/
	
	Current version of nr database is available at the NCBI ftp site.
	ftp://ftp.ncbi.nih.gov/blast/db/
	

II. INSTALLATION INSTRUCTIONS

  1. Download and/or install prerequisites described above.

  2. Unpack the distribution:
	
	$ tar zxvf provean-1.1.5.tar.gz

  3. Change to the source directory:

	$ cd provean-1.1.5

  4. Run configure:

	Specifying the location of NCBI nr databases (including an alias file
	name without the file extension, e.g. nr for nr.pal)
	
	$ ./configure BLAST_DB=/path/to/blast/database/nr

	or if you do not have nr databases yet (In this case, you can set the
	location later manually.)

	$ ./configure

	By default this will place the PROVEAN binary files in /usr/local/bin
	and the associated files in /usr/local/share/data/provean or
	/usr/local/share/doc/provean. If you want to place the executables in
	a different location, or you do not have write permissions (i.e. are
	not the root or superuser) to these directories, then you may specify
	a different location using:

	$ ./configure --prefix=/other/path/

	If you receive an error during configuration that psiblast, blastdbcmd,
	or cd-hit cannot be found, and you have indeed installed it, then place
	it in your PATH variable, or specify its location (including the binary
	file names):

	e.g.)
	$ ./configure PSIBLAST=/path/to/psiblast
	$ ./configure CDHIT=/path/to/cd-hit
	$ ./configure PSIBLAST=/path/to/psiblast BLASTDBCMD=/path/to/blastdbcmd CDHIT=/path/to/cd-hit
	
  5. Compile:

	$ make

  6. Install the software and data:

	$ make install


III. SETTING UP DATABASE & SOFTWARE LOCATION
	
	If you did not specify the location of NCBI nr database during
	installation, then you may edit provean.sh file in /usr/local/bin or
	under your specified directory with --prefix option. Set BLAST_DB
	variable accordingly.

	You can also change the path to psiblast, blastdbcmd, or cd-hit in a
	similar way if you want to use different version.

	e.g.)
	BLAST_DB="/path/to/blast/database/nr"
	PSIBLAST="/path/to/psiblast"
	BLASTDBCMD="/path/to/blastdbcmd"
	CD_HIT="/path/to/cd-hit"


IV. RUNNING PROVEAN

  1. To see PROVEAN usage instruction, execute provean.sh with -h option:

	$ provean.sh -h
	PROVEAN v1.1.5

	USAGE:
	 provean.sh [Options]

	Example:
	 # Given a query sequence in aaa.fasta file, 
	 # compute scores for variations in bbb.var file 
	   provean.sh -q aaa.fasta -v bbb.var

	Required arguments:
	 -q <string>, --query <string>
	   Query protein sequence filename in fasta format
	 -v <string>, --variation <string>
	   Variation filename containing a list of variations:
	     one entry per line in HGVS notation,
	     e.g.: G105C, F508del, Q49dup, Q49_P50insC, Q49_R52delinsLI

	Optional arguments:
	 --save_supporting_set <string>
	   Saves supporting sequence set infomation into a given filename
	 --supporting_set <string>
	   Supporting sequence set filename saved with
	   '--save_supporting_set' option above
	   (This will save time for BLAST search and clustering.)
	 --tmp_dir <string>
	   Temporary directory used to store temporary files
	 --num_threads <integer>
	   Number of threads (CPUs) to use in BLAST search
	 -V, --verbose
	   Verbosely shows the information about procedure
	 -h, --help
	   Gives this help message

  2. Run PROVEAN with test examples:

	Test examples are provided in share/data/provean/examples directory.
	Change to the directory and execute provean.sh with the options
	shown below. You should get a similar result, but the scores could be
	different since they depend on the protein database used.

	$ cd /usr/local/share/data/provean/examples
	$ provean.sh -q P04637.fasta -v P04637.var --save_supporting_set P04637.sss
	## PROVEAN v1.1 output ##
	# Query sequence file:  P04637.fasta
	# Variation file:       P04637.var
	# Protein database:     /usr/local/projects/SIFT/ychoi/provean_genome/provean_result/nr_db/nr
	[13:50:20] searching related sequences...
	[14:02:56] clustering subject sequences...
	[14:02:57] supporting sequence set was saved at P04637.sss
	[14:02:57] supporting sequences were saved in FASTA format at P04637.sss.fasta
	# Number of clusters:   30
	# Number of supporting sequences used:  413
	[14:02:57] computing delta alignment scores...
	## PROVEAN scores ##
	# VARIATION     SCORE
	P72R    -0.461
	G105C   -8.119
	K370del -2.201
	H178_H179insPHP -10.945
	L22_W23delinsQS -10.392

	In the above case, you provided a filename with '--save_supporting_set'
	option so that the infomation on supporting sequence set is stored
	into a file. You can provide this file with '--supporting_set' option
	so that PROVEAN skips the BLAST search and clustering procedures to
	save time as below.

	$ provean.sh -q P04637.fasta -v P04637.var --supporting_set P04637.sss
	## PROVEAN v1.1 output ##
	# Query sequence file:  P04637.fasta
	# Variation file:       P04637.var
	# Protein database:     /usr/local/projects/SIFT/ychoi/provean_genome/provean_result/nr_db/nr
	# Number of clusters:   30
	# Number of supporting sequences used:  413
	[14:05:40] computing delta alignment scores...
	## PROVEAN scores ##
	# VARIATION     SCORE
	P72R    -0.461
	G105C   -8.119
	K370del -2.201
	H178_H179insPHP -10.945
	L22_W23delinsQS -10.392

  3. Interpreting PROVEAN scores:

	We suggest using a cutoff of -2.5 for the PROVEAN score when using the
	NCBI nr protein database released in August 2011. That is, consider a score
	higher than -2.5 to be neutral (tolerated) and that lower than or equal to
	-2.5 to be deleterious (damaging). The PROVEAN scores and optimal cutoff
	may slightly vary with different versions of nr database because the scores
	are computed based on the homologs in the DB. More detailed information on
	PROVEAN scores can be found at http://provean.jcvi.org/about.php



Yongwook Choi
ychoi@jcvi.org
```
