
```
RepeatMasker
Developed by Arian Smit and Robert Hubley
Please refer to: Smit, AFA, Hubley, R. & Green, P "RepeatMasker" at
http://www.repeatmasker.org

NOTE: GitHub 'master' branch is a development branch and not a release.  
      Some large data files such as the Taxonomy database, TE metadata
      and default libraries are not in the repository.  Please download 
      an official release from:

        http://www.repeatmasker.org/RMDownload.html

```

RepeatMasker
------------

RepeatMasker is a program that screens DNA sequences for interspersed
repeats and low complexity DNA sequences. The output of the program is
a detailed annotation of the repeats that are present in the query
sequence as well as a modified version of the query sequence in which
all the annotated repeats have been masked (default: replaced by
Ns). Sequence comparisons in RepeatMasker are performed by the program
cross_match, an efficient implementation of the Smith-Waterman-Gotoh
algorithm developed by Phil Green, or by WU-Blast developed by Warren
Gish.

See "INSTALL" for instructions on how to install RepeatMasker.
See "repeatmaker.help" for a detailed program manual.

Libraries Overview
------------------

Updates of the RepeatMasker program are distributed with a copy of the
Dfam database ( www.dfam.org ) and a copy of the new Dfam_consensus 
database ( site coming soon ).  Dfam and Dfam_consensus are small but
growing "open" databases of Transposable Element profile hidden markov
models and consensus sequences respectively.  

RepeatMasker is also compatible with the RepBase database managed by 
the Genetic Information Research Institute and requires a license to 
use.  We maintain the "Repbase RepeatMasker Edition" libraries as 
co-editor of RepBase Update and aim to keep them in sync with the 
RepBase Update libraries.  However, at any one time there are 
differences.  Entries can differ somewhat in sequence, generally not 
by more than a few percent.  The nomenclature is by and large identical.
Discrepant RepBase and RepeatMasker names for identical sequences 
are indicated in the EMBL formatted version of the RepeatMasker database.  

An inevitable origin of differences is RepeatMasker's extensive
post-alignment processing (improvement) of the repeat annotation. For
one of many examples, internal sequences of LTR elements can be named
after the LTRs, even if there is no specific entry for that element in
the databases. The new alignment file format ( 4.0 and higher ) displays
the actual consensus sequence identifier and a ID to cross-reference to
the final annotation in the *.out file.

For more information about RepBase please go to http://www.girinst.org 
and look for the "RepBase RepeatMasker Edition".

The script queryRepeatDatabase.pl in the utilities directory
( RepeatMasker/util/ ) will display how many repeats are in the library
for a given species. For poorly-covered species, you can create your
own libraries and use these with RepeatMasker. Alternatively, you can
mask your sequence using our transposon-protein database using the
RepeatProtienMask utility included with the RepeatMasker package.  
This utility uses the repeat protein database which is also included
in the RepeatMasker package.

RepeatMasker "open-4.0" and later versions are distributed under the
Open Source License.  Please read LICENSE file for more information.

