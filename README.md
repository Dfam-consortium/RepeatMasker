
```
RepeatMasker
Developed by Arian Smit and Robert Hubley
Please refer to: Smit, AFA, Hubley, R. & Green, P "RepeatMasker" at
http://www.repeatmasker.org
```

RepeatMasker
------------

RepeatMasker is a program that screens DNA sequences for interspersed
repeats and low complexity DNA sequences. The output of the program is
a detailed annotation of the repeats that are present in the query
sequence as well as a modified version of the query sequence in which
all the annotated repeats have been masked (default: replaced by
Ns). Sequence comparisons in RepeatMasker are performed by one of
several available alignment programs:

  - RMBlast, a variant of NCBI blastn that supports substitution 
    matrices, complexity adjusted scoring and masklevel filtering.
  - crossmatch, an efficient implementation of the Smith-Waterman-Gotoh
    algorithm developed by Phil Green.
  - NHMMER, a profile Hidden Markov Model aligner written by Travis
    Wheeler and Sean Eddy.
  - ABBLAST, A blast variant developed by Warren Gish.

See "repeatmasker.help" for a detailed program manual.

RepeatMasker "open-4.0" and later versions are distributed under the
Open Source License.  Please read LICENSE for more information.

Libraries Overview
------------------

RepeatMasker works out-of-the-box with user-supplied libraries provided
via the `-lib` option: FASTA files for use with RMBlast, crossmatch, or
ABBLAST, and profile HMM files for use with NHMMER.

For automated, species/taxa-specific queries against the Dfam database,
RepeatMasker supports FamDB as an optional (but highly recommended)
dependency. FamDB manages Dfam library partitions and can generate
organism-specific consensus or HMM libraries on the fly. The FamDB
project and installation instructions are at:

  https://github.com/Dfam-consortium/famdb

The FamDB project also documents how to combine RepBase sequences with
Dfam. RepeatMasker is compatible with RepBase data, but merging RepBase
with FamDB is handled entirely through the FamDB installation process.


Installation
------------

### Prerequisites

  - A UNIX based operating system.

  - Perl 5.8.0 or higher.

  - Python 3.0 or higher.

  - TRF 4.09 or higher ( http://tandem.bu.edu/trf/trf.html )

  - A search engine — at least one of the following is required:

       RMBlast    :  http://www.repeatmasker.org/RMBlast.html
       crossmatch :  http://www.phrap.org
       NHMMER     :  https://hmmer.org  ( Dfam/FamDB required )
       ABBLAST    :  http://blast.advbiocomp.com/licensing/

  - FamDB (optional, but highly recommended) for species-specific
    queries against the Dfam TE database:

       https://github.com/Dfam-consortium/famdb

    Follow the FamDB installation instructions to download and install
    Dfam library partitions into the RepeatMasker `Libraries/famdb/`
    directory.

### Installing RepeatMasker

  1. Unpack the distribution in the desired location (e.g. `/usr/local/`).
     Do not extract into a directory that already contains a `RepeatMasker`
     subdirectory, as it will attempt to overwrite existing files. For
     example:

         % cp RepeatMasker-open-4-#-#.tar.gz /usr/local
         % cd /usr/local
         % gunzip RepeatMasker-open-4-#-#.tar.gz
         % tar xvf RepeatMasker-open-4-#-#.tar

  2. RepeatMasker is not distributed with a TE library. You can use it
     immediately with a custom library (`-lib mylib.fa`), or install FamDB
     and Dfam library partitions for automated species-specific annotation.
     See the FamDB releases page for downloadable partitions:

         https://github.com/Dfam-consortium/FamDB/releases

  3. Configure the distribution by running the configure script:

         % perl ./configure

     The configure script will prompt for the locations of the search
     engine(s) and any optional dependencies.


Library Cache Directories
-------------------------

Since version 3.0, RepeatMasker creates a cache of species-specific libraries 
extracted from FamDB to speed up repeated searches.  It uses the first 
writable directory in the following path:

  1. The `Libraries/` subdirectory of the RepeatMasker installation.
  2. The `.RepeatMaskerCache` subdirectory of the user's home directory.
  3. The temporary processing directory `RM_#` created alongside the
     sequence file and removed at the end of the run.

If the cache cannot be written to paths 1 or 2, libraries are rebuilt on
every run, which will slow down jobs on shorter sequences.


