
```
RepeatMasker
Developed by Arian Smit and Robert Hubley
Please refer to: Smit, AFA, Hubley, R. & Green, P "RepeatMasker" at
http://www.repeatmasker.org

IMPORTANT:
The github 'master' branch does not contain any TE libraries.  As of
4.1.7, you may use the release as-is for custom library searches 
(e.g "-lib mylib.fasta").  For searches against community developed
TE libraries, you will need to download and install one or more of the
Dfam database FamDB partition.  See the RepeatMasker website for more
information: http://www.repeatmasker.org/RepeatMasker/
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

See "INSTALL" for instructions on how to install RepeatMasker.
See "repeatmaker.help" for a detailed program manual.

Libraries Overview
------------------

Updates of the RepeatMasker program are distributed with a copy of the
Dfam database ( www.dfam.org ). Dfam is an "open" databases of 
Transposable Element seed alignments, profile Hidden Markov Models 
and consensus sequences.

RepeatMasker is also compatible with the RepBase database managed by 
the Genetic Information Research Institute and requires a license to 
use. Up until 2019 we maintained the "Repbase RepeatMasker Edition" 
libraries as co-editor of RepBase Update.  For newer versions of 
RepBase users will need to use the sequences in FASTA format with
RepeatMasker's "-lib" option.

RepeatMasker "open-4.0" and later versions are distributed under the
Open Source License.  Please read LICENSE file for more information.

