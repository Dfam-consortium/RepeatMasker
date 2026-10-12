Next
  - ProcessRepeats memory/processing time fixes for scaffold assemblies (high contig
    counts)
  - RepeatMasker reports the FamDB data directory in use ("Data Dir") next to
    the famdb.py location.
  - When famdb.py fails, RepeatMasker now shows what it wrote to stderr (e.g.
    a Python traceback) along with the command and exit status.  Taxonomy.pm
    used to treat a crashed famdb.py as a valid empty answer, so a broken
    FamDB install (e.g. h5py missing from the python on PATH) made -species
    runs silently fall back to generic, non-mammalian processing.  It now
    stops with famdb.py's output.
  - BUGFIX: Cached species libraries did not record which FamDB partitions
    were installed when they were built.  A -species run made before the
    partition for that taxon was installed cached only the families
    available then (e.g. those of higher-level clades) and kept using them
    after the partition was added.  Each species cache directory is now
    named for the applicable families installed when it was built (e.g.
    homo_sapiens.1a2b3c4d; the counts are kept in famdbSignature.txt), so
    a change in installed partitions selects a new cache instead of
    reusing or overwriting the old one.  Caches made by earlier versions
    are no longer used (they are rebuilt once) and may be deleted.
  - New -famdb_data_dir option ( FAMDB_DATA_DIR environment variable, or
    set once with configure ): the directory holding the FamDB *.h5 files.
    It is passed to famdb.py as -i by RepeatMasker, ProcessRepeats,
    DateRepeats, RepeatProteinMask and util/buildSummary.pl, so a Dfam
    installation can be chosen per run without editing famdb.conf.
  - New -cache_dir option ( REPEATMASKER_CACHE_DIR environment variable,
    or set once with configure ): where to keep the libraries RepeatMasker
    builds from FamDB.  It is tried first, before Libraries/ and
    ~/.RepeatMaskerCache, and is created if missing.  RepeatMasker now says
    when no cache location is writable and libraries are being built in
    the run's temporary directory, to be rebuilt by every run.
  - A cache directory can be shared by concurrent runs.  Libraries are
    built in a directory private to each run and moved into place when
    complete, so a run never removes or reads another run's partial build.
    Builds abandoned by killed runs are cleaned up once their process has
    exited (or after two days for builds from other hosts).
    RepeatProteinMask writes its RepeatPeps library the same way.
  - The -h option list now shows the environment variable for each
    setting that has one.
  - RepeatMasker now accepts the common-name aliases already used by
    ProcessRepeats and DateRepeats (e.g. mouse, rat, drosophila, anopheles,
    worm, mammal, diatom).  Most of these are ambiguous to FamDB and were
    rejected.  All alias targets were checked against Dfam 4.0.  "maize" is
    no longer an alias for the Zea genus; FamDB resolves it to Zea mays,
    which is what RepeatMasker already used.
  - BUGFIX: Taxonomy::isA() returned 0 rather than -1 for a species or clade
    not found in FamDB.  ProcessRepeats tested the result for truth, where
    -1 counts as true, so it now compares against 1 and stops when run
    directly with an unknown -species instead of quietly producing generic
    output.
  - BUGFIX: DateRepeats could not use aliases (Taxonomy::isSpecies()
    applied them after the FamDB lookup) or taxa whose names contain
    "<", ">", "," or parentheses (e.g. "Drosophila <flies,genus>"); the
    stripped form it returned is not a name FamDB knows.  It now returns
    FamDB's sanitized name (e.g. drosophila_flies_genus).
  - BUGFIX: RepeatMasker died with "Can't use string ... as a symbol ref"
    when appending a batch's log or stderr file to the run's (introduced
    with the runtime fix for high sequence counts above).
  - BUGFIX: Simple repeats and low complexity regions could be silently
    missing when RepeatMasker was run from a deep directory.  TRF 4.09 and
    4.10 overflow a buffer when given a path of about 260 characters or
    more and then abort or write nothing, and RepeatMasker passed TRF the
    full path of each batch file.  The run still exited 0, with no
    Simple_repeat or Low_complexity lines in the .out.  TRF is now given
    the file name relative to the directory it runs in.  A TRF run that
    produces no output for a non-empty batch is now treated as a failure
    and reported with TRF's error message; previously only a signal was
    caught, and the shell wrapping TRF hid that too.  (TRF's exit status
    is the number of sequences read, so it cannot be used to detect
    errors.)
  - BUGFIX: ProcessRepeats could split a forward-strand SVA whose VNTR is
    longer than the consensus into 2-4 .out lines sharing one ID, while the
    reverse complement of the same copy was reported as one line.  An SVA
    looks ahead at most 3 candidates when assembling its fragments, and
    only merges made by the SVA/LAVA rule stopped counting against that
    limit.  Forward-strand VNTR fragments often start just before the VNTR
    range and merge through the general "joined pair" rule instead, using
    up the look-ahead before the next fragment is reached.  Merges by every
    rule now free a look-ahead slot.  E.g. hg38 chr2:41746002-41748725
    (SVA_F) is now one line, 53-2716 + SVA_F 1 1375 (0), on both strands.
  - BUGFIX: util/rmToTrackHub.pl exited 0 even when bedToBigBed failed,
    leaving a corrupt or partial bigBed file and a hub that pointed to it.
    It now stops with bedToBigBed's exit status and the command that
    failed, and removes the partial file.  Failures of the sort and
    intersectBed steps it runs are reported the same way.
  - t/RepeatMasker.t now fails when an expected output file is missing;
    it used to pass.
  - BUGFIX: RepeatMasker::createSubElements() reported fragments made only of
    genome bases inserted relative to the consensus when it cut a hit at
    excised simple repeats.  Such a fragment aligns to no consensus bases.
    On the reverse strand it was given an empty consensus range (end one below
    start) and the next fragment's range was shifted by one base; on the
    forward strand it claimed one consensus base it did not align.  These
    fragments are now dropped.  A remaining case of the 4.2.4 fix to the same
    coordinate arithmetic.
  - BUGFIX: The -famdb_dir option and the FAMDB_DIR environment variable were
    ignored for species searches.  Taxonomy.pm read FAMDB_DIR from the
    configuration file when the module loaded, before overrides were applied,
    and ProcessRepeats, DateRepeats and util/buildSummary.pl never applied
    them.  Taxonomy now reads the setting when constructed, RepeatMasker
    passes the resolved location to ProcessRepeats, and the other scripts
    apply the environment override.
  - BUGFIX: Alignments next to excised repeats could have their query
    letters shifted one base from the reported range, with an excision
    marker ("N") as their first letter.  Two causes:
      - Abutting excisions each leave a marker, so two fragments of a hit
        can be separated by more than one.  createSubElements() always
        stepped past exactly one, leaving the rest at the front of the next
        fragment.  fragmentAnnotations() now records how many markers
        precede each fragment, and createSubElements() steps past exactly
        those, moving the consensus position by the consensus letters
        actually aligned to them.
      - addToRefinementCollection() compared excision coordinates with hit
        coordinates in the cut sequence, and skipped markers next to the
        start of the hit or to another excision, so refinement seeds (e.g.
        Alu subfamilies) could keep a marker.  It now places each marker in
        the cut sequence and removes every marker inside the hit.
    On a 1 Mb human test region 6 fragment and 5 refined alignments
    changed, all by one consensus base; masking was unchanged.
    t/fragmentCoords.t covers both cases.

4.2.5
  - ProcessRepeats speed improvements(25-45% faster on mammalian genomes).
  - Search engines now prepare their own subject databases.  Added
    prepareSubject(), isSubjectPrepared() and getSubjectArtifacts() to
    SearchEngineI, implemented for the blast and HMMER engines.
  - RepeatProteinMask looks for blastx and makeblastdb on demand, first in
    RMBLAST_DIR and then on PATH, and reports what to install if neither
    has them.  These are NCBI toolkit programs that rmblastn 3.x does not
    ship.
  - BUGFIX: ProcessRepeats scoreLINEPair() tested getQueryRemaining()
    where it should have tested getSubjRemaining() when deciding whether a
    reverse-orientation LINE pair was compatible with the model.
  - BUGFIX: ProcessRepeats isTransposonPairCompatible() only scored pairs
    on the "C" strand.  The "+" strand case was missing, along with the
    consensus length needed to write it, so forward-strand DNA transposon
    fragments that should have joined were reported as separate elements.
  - Native support for the rmblastn 3.x series.  The reimplemented
    rmblastn takes GNU-style long options and searches a FASTA file
    directly, with no makeblastdb step.  NCBIBlastSearchEngine::new()
    probes the binary and re-blesses the object into a new
    RMBlastSearchEngine class at major version 3 or higher, so callers
    ask for the same engine and get whichever implementation matches the
    installed binary.  rmblastn 2.x still works and gives the same
    results as before.

  
4.2.4
  - Support for new component-based FamDB and Dfam 4.0
  - Updated documentation
  - An updated configure script without screen clears supports easier 
    history tracking and debugging.
  - FamDB is now a fully separate and optional package for RepeatMasker.  For some
    time RepeatMasker has been able to run custom library searches without the need
    for a local copy of the Dfam database.  This step simply removes the need to
    update FamDB every time RepeatMasker is updated.
  - BUGFIX: Negative subject coordinates were produced when splitting alignments
    in RepeatMasker::createSubElements().  The coordinate arithmetic contained
    an off-by-one error and incorrectly assumed the X/N separator always consumed
    one subject position, causing coordinates to drift negative when the separator
    aligned to a gap in the subject.  Reported by Mark Diekhans.

4.2.3
  - BUGFIX: Since the 4.1.1 release the staged approach to searching
    TE families has not been applied species in the 'rodentia' order
    due to a flag being set incorrectly.  While this doesn't prohibit
    families from being annotated, it does decrease the sensitivity of
    the search.  This release fixes this issue.  This problem did not 
    impact any other orders.
  - Feature Request: Added support for BGZIP files by changing the 
    method used to process compressed files.  Instead of uncompressing
    input files in-place, the program now leaves the original 
    compressed file intact.

4.2.2
  - BUGFIX: In the 4.2.1 release a bug was introduced that impacted
    the handling of poly-A tails and the naming of LINE annotations.
    While the overall masking level remained the same, the change 
    caused some previously hidden and overlapping annotations to 
    appear with invalid query ranges ( query start > query end ). 
    (github #362)
  - Added support for '+' symbol in input filenames. (github #365)

4.2.1
  - BUGFIX: An infinite loop was possible in rare circumstances
    involoving interleaved satellite consensi while using the
    "-a" option.  This caused the alignment file to keep growing
    and eventually ProcessRepeats to crash.
  - BUGFIX: Fix a one-off issue in alignment range calculations
    induced by the clipping mechanism in RepeatMasker.
  - BUGFIX: The configure script was crashing while attempting to
    import Repbase in the previous version.
  - Updated to the latest famdb version (2.0.5)
  - Crossmatch alignment format parser hardening.  While not proper, 
    negative coordinates can appear in the alignment file due to 
    incorrect range calculations.  This change allows for the parser
    to handle these cases and not mis-align the data columns.
  - Updated the length of the Charlie7 family in RepeatAnnotationData.pm
  
4.2.0
  - BUGFIX: For mammalian genomes searched using the "-species" option,
    newer Dfam families without search stages, or with the "-uncurated"
    flag Dfam DR families were not being included in the search.  
  - BUGFIX: RepeatMasker exits with an error code if no repetitive
    sequences were found in the input file.
  - BUGFIX: famdb.py was not reporting the correct counts for families
    ancestral to the query term.  RepeatMasker uses this count and therefore
    was also reporting an incorrect value.
  - Remove duplicate families when attempting to merge Repbase 
    RepeatMaskerEdition into FamDB.
  - The configure utility was not using the -libdir path when creating
    the RepeatMasker.lib file. (github #344)

4.1.9
  - BUGFIX: When the RepeatMasker adjudicator (ProcessRepeats) is 
    faced with a huge number of overlapping and redundant alignments 
    the joining can go awry, creating a infinite loop.  The effect
    is that ProcessRepeats never completes.
  - BUGFIX: When using nhmmer/pHMMs and a non-species-level "-species"
    parameters (e.g. "-species mammals"), RepeatMasker was using an
    exteremly permissive per-family score cutoff.  This overwhelmed
    ProcessRepeats with many low-scoring alignments.
  - BUGFIX: ProcessRepeats dies around line 8190 with a "Division by
    zero" error.
  - BUGFIX: Fixed call to famdb for buildSummary.pl in the -species 
    option.  Updated documentation.
  - Updated repeatmasker.help.

4.1.8
  - Fix table format for non-mammals.  The "Retroelements" heading
    didn't include Penelope elements ( "PLE/" ) in it's tabulation.
  - ProcessRepeats - A major change in how annotations are processed
    to reduce PR's memory footprint.  Instead of loading all anntations
    in memory, they are now streamed and processed a sequnce at a time.
    For chromosomal assemblies this reduced the memory footprint to
    the size of the largest chromosome, instead of the whole genome.
  - Fixed a bug that caused general/is.lib to be left in a unfrozen
    state for rmblastn.
  - Updated to Dfam3.9

4.1.7-p1
  - This is a critical fix to the 4.1.7 release.  The -species option
    was not being interpreted correctly by RepeatMasker/ProcessRepeats.
    The most obvious change for mammals/primates was that the output
    format (*.tbl file) looked different than it had before.  Under the
    hood this also changed the way some species were handled in the 
    library searches. 

4.1.7
  - Handle the case were the default search engine is not configured
    and require the user to provide it using -e/-engine.
  - Use the global alignment bandwidth simulation for the refinement
    step.
  - Allow for minimal database install and warn when used with -species.
  - Fixed a bug on configure with setting the default search engine when
    only one search engine is configured.
  - RepeatProteinMask wasn't obtaining the correct search engine if
    configure didn't set it up correctly.
  - FamDB wasn't exporting RepeatMasker.lib with class-name suffixes.
    This had an impact on the RepeatModeler/RepeatClassifier package
    but not on RepeatMasker itself. Upgraded to FamDB 1.0.5
  - Fixed an issue with parsing RepBase.  Some records were not being
    parsed correctly for records that do not include the keyword "repbase"
    following the identifier.
  - New semantantics for NCBIBlastSearchEngine bandwidth interpretation.
  - Fixed a divide-by-zero error in ProcessRepeat's merging of 
    DNA tranpsoson overlapping fragments.
  - FamDB fixed inconsistency between curated/uncurated counts vs
    exports.

4.1.6
  - Upgraded to FamDB 1.0.2 to support Dfam 3.8 and the new partitioned
    database format.
  - Added Libraries/RMRB_spec_to_tax.json to project.  This maps the
    RepBase taxanames to current NCBI tax_ids and needs to be refreshed
    with each new Dfam release.
  - Added softmasking support to NCBIBlastSearchEngine.pm.
  - Added new '--uncurated' flag to handle single export Dfam format.     
    If this flag is used the CONS/HMM cached directories will be suffixed
    with "_wunc".
  - Fixed sunk error messages from famdb.py.  Now they will be displayed
    and cause RepeatMasker to quit. 
  - Additional library setup steps and error checking for configure utility.
  - CAF documentation in SearchResult.
  - calcDivergenceFromAlign clarified use of "-a" in documentation.

4.1.5
  - Updated codebase for Dfam 3.7 compatibility (famdb format 4.3)
  - Penelope classification change caused *.tbl file accounting
    to place them in the Unknown category.  Also fixed landscape
    generation tool.
  - Added a new utility to merge *.out *.align files generated by
    running RepeatMasker serially.
  - Repbase Metadata is seriously out-of-date, updated species names
    so that they match the current NCBI Taxonomy names.
  - Fixed an issue with the HMM parser.  It wasn't recognizing 
    negative values for Tau with models that do not have GA thresholds.

4.1.4
  - Added support for RMBlast 2.13.0
  - Release of the TE genome browser visualization (UCSC) and trackhub
    generation tool.
  - New CpGSites and unadjusted Kimura stats in the *.align file.
  - Fixed a bug that caused the read-only state of the input fasta file
    to propogate to the intermediate files and cause the program to exit.
  - Removed DateRepeats as it's based on old library formats - this functionality
    will return with the refactored version of RM in the works.

4.1.3-p2
  - A recent change in 4.1.3 to correct blank fragment ID fields can in rare
    cases causing the error message: 'Can't call method "setLeftLinkedHit"'.
  - The RepeatAnnotationData.pm file containing necessary information for
    recognizing equivalent fragments of DNA transposons was missing data.
  - The MULE-MuDR class was added to the *.tbl file for "-lib" searches.

4.1.3
  - A new utility for generating trackHubs for our new UCSC TE visualization
  - Fix a bug where killing RM while starting up can leave the cached libraries
    in an inconsistent state. 
  - Fixed a bug where in rare cases the joined fragment ID field is blank
  - Merged in changes to Dupmasker supporting multi-threaded use
  - Fixed legacy RepBase taxonomic labels
  - Added support for GFF v3 output and fixed the utility/rmOutToGFF3.pl

4.1.2-p1
  - Releases 4.1.1-4.1.2 contained a bug with the processing of Alu sequences 
    in primates. The step where an initial annotation is refined into a 
    particular Alu subfamily was not performed and the annotations remained 
    labeled with the initial capture sequence ( AluJb, AluSx, or AluY ). 
    This patch release fixes this one issue.
  
4.1.2
  - Fixed 21 protein family classifications in RepeatProteinLib.
  - Fixed a problem with the generation of the RepeatMasker.lib file for use by 
    RepeatModeler. In release 4.1.1 it did not add the classification info to 
    this auxilary file.
  - Fixed a "log(0)" error that can cause the program to fault in rare circumstances.
  - buildSummary now supports FamDB and has improved documentation.
  - Bugfixes and improvements to FamDB.

4.1.1
  - Dfam (starting with version 3.2) is now distributed in the FamDB file format 
    based on HDF5, which has improved support for large datasets compared to the 
    EMBL and HMM formats that were previously used. RepeatMasker therefore includes
    a copy of famdb.py, and depends on the python package h5py.

      o The 'configure' script and other parts of RepeatMasker have been updated to 
        accomodate these changes.
      o The utilities 'queryTaxonomyDatabase.pl' and 'queryRepeatDatabase.pl' are no 
        longer included, since that data is now included in FamDB. The 'famdb.py' 
        tool can be used to make many of the same queries as the removed utilies, 
        and even more.

4.1.0
  - Refactored the RepeatMasker configuration system.  
    The changes make it easier to automate installation and/or
    bundle RepeatMasker.

4.0.9-p2
  - Another improvement in TRF error handling

4.0.9-p1 : Tue Apr 16 08:51:22 PDT 2019
  - Input files containing multiple FASTA sequences
    caused RepeatMasker to error out with a message
    like:

      "WARNING: TRF returned an error (Return code = ### )
       TRF parameters: 2.7.7.80.10.50.10
       A search phase could not complete on this batch.
       The batch file will be re-run and if possible the
       program will resume.
       WARNING: Retrying batch ( 1 ) [ 255,, 195]..."

    This bug was introduced in an attempt to improve TRF
    error catching.  Unfortunatly the return codes are
    not documented for TRF and the assumption that 256
    is the only successful return code is wrong.  The
    "success" code appears to change depending on the
    number of sequences in the file.  The workaround is
    to fail only if there is a message in the error output
    file.

