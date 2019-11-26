README file for Codon Evolution Simulation (CES)

*************************************************************************************

VERSION: CES 1.0

RELEASE DATE: Sept 26, 2009


DESCRIPTION:
	The code simulates the evolution of synonymous codon usage under selection against
	translational nonsense errors.

	The model is described in:
	"Gilchrist, M.A., P. Shah, and R. Zaretzki (2009 or 2010 depending on publication date). 
	Measuring and detecting molecular adaptation in codon usage against nonsense errors 
	during protein translation. Genetics"

	****When using this work please cite the above paper****


	NOTE: folders for output files must exist prior to running CES or else the
	program will crash. 
		


LICENSE: 
	Licensed under GNU General Public License. See LICENSE file for details.

BINARIES: *NIX and OSX:
	Precompiled binary for *NIX machines is included in bin/.  
	These binaries were compiled on a linux machine with GCC 4.3.3
	The binaries should run on all i386 and x86_64 machines running linux 
	or OSX with the source code for the GNU Scientific Libraries (GSL) installed on them.  
	If problems are encountered we encourage you to try recompiling the code for your local 
	machine before contacting the authors.

COMPILERS:
	The code can be recompiled from source and optimized for the hardware of the 
	machine it will be run on.
	
	The source code has been successfully compiled with the following compilers.
	Mac:		gcc, g++
	Linux:		g++, gcc -lm
	
	In addition, the development libraries of the GNU Scientific Libraries must be installed.

BUILD:	Builds on Ubuntu 9.04 machine with GSL libraries using command: 

	       g++   ./CES.c  -g -lgsl -lgslcblas -lm  -o ../bin/CES

	For use with openmp libraries which allow the use of multiple cpus that share memory
	include the argument: -fopenmp when  compiling

	Example shell scripts for building executables can be found in ./source


SYNOPSIS:
	./CES [options] -F <sequence_file>

OPTIONS:

	-F <FILE>	Location of the fasta file to be used for simulation. Note the program
			will likely crash if genes with internal stop codons are used.
			phi values must be included on the same line as the ORF name
			e.g.
			>YAL001C	phi=	0.0088862 
			Note the use of TABS between the ORF name and "phi=" and the phi value
			[DEFAULT]  -F fasta/S.cerevisiae.S288c.with.phi.fasta

	-A1 <REAL>	Specify ribosome initation cost in ATPs,  a1 in SEMPPR
			[DEFAULT] a1=4

	-A2 <REAL>	Specify ribosome elongation cost in ATPs,  a2 in SEMPPR
			[DEFAULT] a2=4
	
	-AT <REAL>	Specify genome AT bias (0.0-1.0).
			[DEFAULT] at=0.5		

	-B <REAL>	Specify the B parameter value
			[DEFAULT] B=0.0025

	-D <BOOLEAN>	Indicate whether to print out delta_eta files for each evolutionary step.
			Print outs can get quite large
			[DEFAULT] D=0.  

	-I <INT>	Specify whether CES should relax selection on the last <INT> 
			amino acids of each sequence.  Used because it is thought that
			most genes can lose a few aa at the end, but still be 
			functional.

	-M <INT>	Indicates simulation time. Simulation time is scaled by the mutation rate mu,  
			where mu = 1E-9. 
			If the argument X is > 0, then we expect to see X * mu substitutions/codon.
			If X < 0, then we run for -X/mu time or, in other words, we expect to see -X
			substitutions if the genes were evolving neutrally.
			[DEFAULT] = -20 

	-Ne <REAL>	Effective population size. 
			[Default] = 1.36E7
			
		
	-O <STRING>	Specifies the folder for output files as well as prefix for specific
			output and log files.
			[DEFAULT] "output/out"


	-T <STRING>	Specify the location of the tRNA abundance/codon 
			translation file.
			[DEFAULT] "tRNA_files/S.cerevisiae.tRNA.tsv"

	-V		Ratio of transition to transversion mutation rates.
			[DEFAULT] 1.0 



FILE FORMATS:
	INPUT FILES:
		Sequence file:
			* The file should contain genes in the standard FASTA format.
			* The analysis will terminate at the first non-sense codon or
			  if the codon contains any character apart from "A", "T", "G" and "C".
			* All nucleotides should be in uppercase characters.
			* The sequence should contain only an ORF WITHOUT the UTRs and intronic
			  regions.

		tRNA abundance/codon translation file:
			* The codon for a given amino acid should be in a
			  consecutive order
			* The amino acids need not be in any particular order

			[AA]	[CODON]	[ABUNDANCE/RATE]

			example:
			A	GCG	7.82
			A	GCT	18.2
			C	TGC	7.47
			.	.	.
			.	.	.
			.	.	.
			W	TGG	11.7
			Y	TAC	17.0
			Y	TAT	10.4

	OUTPUT FILES:

		*.ces.fasta
			A fasta file with the evolved sequences


UPDATES:
	Updates for this code can be found at the following website: 
	www.tiem.utk.edu/~mikeg/SupplementaryMaterials/CES/ces.html

BUGS:
	In case of any bugs or trouble with the code, send a mail to mikeg@utk.edu
