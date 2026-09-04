#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) NCBIBlastSearchEngine.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of SearchEngineI for the
##      the NCBI Blast search engine.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2004 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
# bless(
#      'NCBIBlastSearchEngine' );
#
###############################################################################
# ChangeLog
#
#     $Log$
#
###############################################################################
# To Do:
#
#

=head1 NAME

NCBIBlastSearchEngine

=head1 SYNOPSIS

use NCBIBlastSearchEngine

Usage: 

  use SearchEngineI;
  use NCBIBlastSearchEngine;
  use SearchResultCollection;

  my $NCBIEngine = NCBIBlastSearchEngine->new( 
                    pathToEngine=>"/usr/local/ncbi/bin/rmblastn" );

  $NCBIEngine->setMatrix( "/users/bob/simple.matrix" );
  $NCBIEngine->setQuery( "/users/bob/query.fasta" );
  $NCBIEngine->setSubject( "/users/bob/subject.fasta" );
  my $searchResults = $NCBIEngine->search();

=head1 DESCRIPTION

  A concrete implementation of the abstract class / interface SearchEngineI
  which use the NCBI rmblastn sequence search engine.

=head1 INSTANCE METHODS

=cut 

package NCBIBlastSearchEngine;
use strict;
use SearchEngineI;
use SearchResultCollection;
use Data::Dumper;
use FileHandle;
use File::Basename;
use File::Spec;
use Carp;
# For debugging only
#use Time::HiRes qw(gettimeofday); 

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter SearchEngineI);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Version
#
my $VERSION = 0.1;
my $CLASS   = "NCBIBlastSearchEngine";

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  croak $CLASS
      . "::new: Missing path to search engine!\n\n"
      . "use \$searchEngine = $CLASS->new( pathToEngine=>\"/usr/local/"
      . "bin/rmblastn\")\n"
      if ( not defined $nameValuePairs{'pathToEngine'} );

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  $this->setPathToEngine( $nameValuePairs{'pathToEngine'} );

  #
  # Version dispatch.  The 3.x series takes the same settings but spells
  # them differently on the command line and has no separate database
  # step, so hand the caller a subclass that knows how to drive it.
  # Substituting here rather than at every call site keeps existing code
  # working: callers still ask for an NCBIBlastSearchEngine, and
  # isa("NCBIBlastSearchEngine") tests still match.
  #
  if ( $this->getMajorVersion() >= 3 && ref( $this ) eq $CLASS ) {
    require RMBlastSearchEngine;
    bless $this, "RMBlastSearchEngine";
  }

  # TODO: Figure out a better design
  $this->setUseDustSeg( 1 );

  # Allow import of values
  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak( $CLASS . "::set: Instance variable $name doesn't exist." . "" );
      }
      $this->$method( $value );
    }
  }

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setOverrideParameters()

  Use: my $value    = getOverrideParameters( );
  Use: my $oldValue = setOverrideParameters( $value );

  Get/Set the the override paramters.  These are used instead
  of all the SearchEngineI default parameters if set.

=cut

##-------------------------------------------------------------------------##
sub getOverrideParameters {
  my $this = shift;

  return $this->{'overrideParameters'};
}

sub setOverrideParameters {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'overrideParameters'};
  $this->{'overrideParameters'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 MaskLevelSequence()

  Use: my $value    = getMaskLevelSequence( );
  Use: my $oldValue = setMaskLevelSequence( $value );

  Get/Set the MaskLevelSequence paramter.  This is the
  sequence ( SearchResult::Query or SearchResult::Subject )
  which will be considered when applying the mask level.  The
  default is SearchResult::Query.

=cut

##-------------------------------------------------------------------------##
sub getMaskLevelSequence {
  my $this = shift;

  return $this->{'maskLevelSequence'};
}

sub setMaskLevelSequence {
  my $this  = shift;
  my $value = shift;

  croak $CLASS
      . "::setMaskLevelSequence: Invalid value ( $value ). "
      . "Should be either SearchResult::Query or "
      . "SearchResult::Subject\n"
      if (    $value != SearchResult::Query
           && $value != SearchResult::Subject );

  my $oldValue = $this->{'maskLevelSequence'};
  $this->{'maskLevelSequence'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 UseDustSeg()

  Use: my $value    = getUseDustSeg( );
  Use: my $oldValue = setUseDustSeg( $value );

  Turn on / off the Dust/Seg screening of words.

=cut

##-------------------------------------------------------------------------##
sub getUseDustSeg {
  my $this = shift;

  return $this->{'useDustSeg'};
}

sub setUseDustSeg {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'useDustSeg'};
  $this->{'useDustSeg'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 TempDir()

  Use: my $value    = getTempDir( );
  Use: my $oldValue = setTempDir( $value );

  Set the directory to use as a temp directory for a search.  
  The default is to use the directory which contains the query sequence.

=cut

##-------------------------------------------------------------------------##
sub getTempDir {
  my $this = shift;

  return $this->{'tempDir'};
}

sub setTempDir {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'tempDir'};
  $this->{'tempDir'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setAdditionalParameters()

  Use: my $value    = getAdditionalParameters( );
  Use: my $oldValue = setAdditionalParameters( $value );

  Get/Set the additional paramters.  These are used in addition
  to the existing parameter set. 

=cut

##-------------------------------------------------------------------------##
sub getAdditionalParameters {
  my $this = shift;

  return $this->{'additionalParameters'};
}

sub setAdditionalParameters {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'additionalParameters'};
  $this->{'additionalParameters'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPathToEngine()

  Use: my $value    = getPathToEngine( );
  Use: my $oldValue = setPathToEngine( $value );

  Get/Set the fully qualified path to the search engine
  binary file.

=cut

##-------------------------------------------------------------------------##
sub getPathToEngine {
  my $this = shift;

  return $this->{'pathToEngine'};
}

sub setPathToEngine {
  my $this  = shift;
  my $value = shift;

  croak $CLASS
      . "::setPathToEngine(): Missing parameter!  Must specify "
      . "a path to the RMBlastN program.\n"
      if ( $value =~ /^\s*$/ );

  croak $CLASS. "::setPathToEngine( $value ): Program does not exist!"
      if ( not -x $value || `which $value` );

  $this->_probeVersion( $value );

  my $oldValue = $this->{'pathToEngine'};
  $this->{'pathToEngine'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
##  Use: $this->_probeVersion( $pathToBinary );
##
##  Ask the binary what it is, and set the version and capability flags.
##  Broken out of setPathToEngine() so that a subclass driving a different
##  rmblastn series can extend it.
##-------------------------------------------------------------------------##
sub _probeVersion {
  my $this  = shift;
  my $value = shift;

  # NCBI toolkit binaries accept "-version".  The 3.x series uses GNU
  # style long options, so try both before giving up.
  my $result = `$value -version 2>&1`;
  if ( $result !~ /rmblast/i ) {
    $result = `$value --version 2>&1`;
  }

  croak $CLASS
      . "::setPathToEngine( $value ): Cannot determine "
      . "engine variant and version!\n"
      if ( $result !~ /rmblast[n]\s*:?\s*(\S.*)/ );

  $this->{'engineName'} = "rmblastn";
  $this->{'version'}    = $1;

  # 2.x reports "2.17.1+" following the NCBI convention; 3.x reports a
  # bare "3.0.4".  Do not require the trailing "+".  This test used to,
  # which meant a 3.x binary parsed as having no version, left
  # hasTabFormat unset, and fell through to the legacy report parser
  # without any error.
  croak $CLASS
      . "::setPathToEngine( $value ): Could not parse a version number "
      . "out of \""
      . $this->{'version'}
      . "\"!\n"
      if ( $this->{'version'} !~ /(\d+)\.(\d+)\.(\d+)/ );

  my $majorVer = $1;
  my $minorVer = $2;
  my $revision = $3;
  $this->{'majorVersion'} = $majorVer;
  $this->{'minorVersion'} = $minorVer;
  $this->{'revision'}     = $revision;

  if ( $majorVer == 2 ) {
    if ( $minorVer >= 13 ) {
      # Since the 2.13.0+ release of RMBlast we now have:
      #    - pre-computed Kimura divergence, Kimura CpG adjusted,
      #      transitions, transversions, cpg_sites, and the cross_match
      #      stats (perc_sub, perc_query_gap, perc_subj_gap).
      #    - Ability to thread on the query sequences rather than just
      #      the subject sequences.
      #    - The ability to output tab delimited format with all the above
      #      fields.
      $this->{'hasQueryThreading'} = 1;
      $this->{'hasTabFormat'}      = 1;
    }
    if ( $minorVer > 14 || ( $minorVer == 14 && $revision >= 1 ) ) {
      $this->{'hasDBSoftMasking'} = 1;
    }
  }
  elsif ( $majorVer == 3 ) {
    # The 3.x series is a reimplementation with GNU style options and no
    # separate database formatting step.  RMBlastSearchEngine drives it;
    # new() substitutes that class for this one.
    $this->{'hasQueryThreading'} = 1;
    $this->{'hasTabFormat'}      = 1;
    $this->{'hasDBSoftMasking'}  = 1;
  }
  else {
    # Refuse rather than guess.  An unrecognized series may differ from
    # both of the ones we know how to drive, and guessing wrong produces
    # bad annotations with no error to show for it.
    croak $CLASS
        . "::setPathToEngine( $value ): Unsupported rmblastn version "
        . $this->{'version'}
        . ".  This release can drive the 2.x and 3.x series.\n";
  }

  return;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getMajorVersion( );

Get the major version number of the engine binary.  2 for the NCBI
toolkit derived series, 3 for the reimplemented series.

=back

=cut

##-------------------------------------------------------------------------##
sub getMajorVersion {
  my $this = shift;

  return $this->{'majorVersion'};
}

##-------------------------------------------------------------------------##

=head2 get_setPathToDBFormatter()

  Use: my $value    = getPathToDBFormatter( );
  Use: my $oldValue = setPathToDBFormatter( $value );

  The program used to format a subject database for this engine.  Derived
  from the engine binary by default, since makeblastdb ships alongside
  rmblastn, so callers need not carry a second configuration value.

=cut

##-------------------------------------------------------------------------##
sub getPathToDBFormatter {
  my $this = shift;

  return $this->{'pathToDBFormatter'}
      if ( defined $this->{'pathToDBFormatter'} );

  return dirname( $this->getPathToEngine() ) . "/makeblastdb";
}

sub setPathToDBFormatter {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'pathToDBFormatter'};
  $this->{'pathToDBFormatter'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 getSubjectArtifacts()

  Use: my @files = getSubjectArtifacts( $path );

  The files makeblastdb produces for a nucleotide database.

=cut

##-------------------------------------------------------------------------##
sub getSubjectArtifacts {
  my $this = shift;
  my $path = shift;

  return () if ( !defined $path );

  # BLASTDB version 4 and 5 suffixes.  Not all are produced for every
  # database, so callers must tolerate absent members.
  return map { "$path.$_" }
      qw( nhr nin nsq ndb not ntf nto njs nog nos nod );
}

##-------------------------------------------------------------------------##

=head2 isSubjectPrepared()

  Use: my $bool = isSubjectPrepared( $path );

  True if $path names a formatted BLAST database.

=cut

##-------------------------------------------------------------------------##
sub isSubjectPrepared {
  my $this = shift;
  my $path = shift;

  return 0 if ( !defined $path );

  return ( -s "$path.nin" || -s "$path.nhr" || -s "$path.nsq" );
}

##-------------------------------------------------------------------------##

=head2 prepareSubject()

  Use: my $subjectPath = prepareSubject( $seqFile,
                                         [outputDir  => $dir],
                                         [dbName     => $name],
                                         [checkStale => 1],
                                         [force      => 1] );

  Run makeblastdb over $seqFile and return the database basename to hand
  to setSubject().  This is not always $seqFile: when outputDir is given
  the artifacts are written there, and the returned path points at that
  directory, which holds the index but not the sequence file.

  By default an existing database is left alone.  Pass checkStale to also
  rebuild when $seqFile is newer than its artifacts, or force to rebuild
  unconditionally.

=cut

##-------------------------------------------------------------------------##
sub prepareSubject {
  my $this    = shift;
  my $seqFile = shift;
  my %params  = @_;

  croak $CLASS
      . "::prepareSubject(): Sequence file ($seqFile) does not "
      . "exist or is empty!\n"
      if ( !-s $seqFile );

  my ( $vol, $dir, $file ) = File::Spec->splitpath( $seqFile );
  my $outputDir = $params{'outputDir'};
  $outputDir = ( $dir eq "" ? "." : $dir ) if ( !defined $outputDir );
  $outputDir =~ s/\/+$//;
  my $dbName = $params{'dbName'};
  $dbName = $file if ( !defined $dbName );
  my $dbPath = "$outputDir/$dbName";

  if ( !$params{'force'} && $this->isSubjectPrepared( $dbPath ) ) {
    my $stale = 0;
    $stale =
        $this->_artifactsAreStale( $seqFile,
                                   $this->getSubjectArtifacts( $dbPath ) )
        if ( $params{'checkStale'} );
    if ( !$stale ) {
      print $CLASS
          . "::prepareSubject(): $dbPath is already prepared, skipping.\n"
          if ( $this->getDEBUG() );
      return $dbPath;
    }
  }

  my $formatter = $this->getPathToDBFormatter();
  croak $CLASS
      . "::prepareSubject(): Cannot find the database formatting program "
      . "($formatter).  A 2.x rmblast installation must provide "
      . "makeblastdb alongside rmblastn.\n"
      if ( !-x $formatter );

  my $log = "$dbPath.makeblastdb.log";
  system( "$formatter -dbtype nucl -out $dbPath -in $seqFile > $log 2>&1" ) == 0
      or croak $CLASS
      . "::prepareSubject(): Error running $formatter on $seqFile.\n"
      . "See $log for details.\n";

  return $dbPath;
}

##-------------------------------------------------------------------------##
## Instance Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getForceLegacyParsing( );

=item Use: my $oldValue = setForceLegacyParsing( $value );

Get/Set the use of legacy parsing when used with rmblastn 2.13.0+ or 
newer.  This serves to assist with debugging.

  $value :  Integer >= 0

=back

=cut

##-------------------------------------------------------------------------##
sub getForceLegacyParsing {
  my $this = shift;

  return $this->{'useLegacyParser'};
}

sub setForceLegacyParsing {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'useLegacyParser'};
  $this->{'useLegcyParser'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getThreadByQuery( );

=item Use: my $oldValue = setThreadByQuery( $value );

Get/Set the use of threading over the query sequences.  This is
used when there is a much larger query set than database.  

  $value :  Integer >= 0

=back

=cut

##-------------------------------------------------------------------------##
sub getThreadByQuery {
  my $this = shift;

  return $this->{'useThreadByQuery'};
}

sub setThreadByQuery {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'useThreadByQuery'};
  $this->{'useThreadByQuery'} = $value;

  return $oldValue;
}



##-------------------------------------------------------------------------##

=head2 getParameters()

  Use: my  $ncbiBlastParamString  = getParameters( );

  Convert object parameters into NCBI rmblastn command line parameters.

=cut

##-------------------------------------------------------------------------##
sub getParameters {
  my $this = shift;

  # Test if engine is available
  my $engine = $this->getPathToEngine();
  if ( !defined $engine || !-f "$engine" ) {
    croak $CLASS
        . "::search: The path to the search engine is undefined or\n"
        . "is set incorrectly: $engine\n";
  }

  my $parameters = $this->_renderParameters( $this->_computeSearchParameters() );

  my $runParameters;
  if ( defined $this->{'overrideParameters'}
       && $this->{'overrideParameters'} ne "" )
  {
    $runParameters = $this->{'overrideParameters'};
  }
  else {
    $runParameters = $parameters . " ";
  }

  if ( defined $this->{'additionalParameters'}
       && $this->{'additionalParameters'} ne "" )
  {
    $runParameters .= " " . $this->{'additionalParameters'};
  }

  return ( "$engine $runParameters" );
}

##-------------------------------------------------------------------------##
##  Use: my $paramsRef = $this->_computeSearchParameters();
##
##  Translate the engine-neutral SearchEngineI settings into the concrete
##  values rmblastn needs, without committing to any particular spelling
##  of the options.  Rendering those values onto a command line is
##  _renderParameters()'s job.
##
##  The 2.x and 3.x series differ in option syntax but not in alignment
##  semantics, so this split lets them share one copy of the score and
##  x-drop translation below.  Changing the numbers here changes both
##  engines; changing option names changes only one.
##-------------------------------------------------------------------------##
sub _computeSearchParameters {
  my $this = shift;

  my %p = ( 'num_alignments' => 9999999 );
  my $value;

  if ( ( $value = $this->getSubject() ) ) {
    # Make sure the subject has been prepared for this engine.  This is
    # checked here rather than in setSubject() because callers are
    # permitted to set the subject before preparing it.
    if ( $this->isSubjectPrepared( $value ) ) {
      $p{'db'} = $value;
    }
    else {
      croak $CLASS
          . "::search: Error...subject database ($value) has not been "
          . "prepared.  Call prepareSubject() before searching.\n";
    }
  }
  else {
    croak $CLASS. "::search: Error subject undefined!\n";
  }

  if ( ( $value = $this->getQuery() ) ) {
    if ( -f $value ) {
      $p{'query'} = $value;
    }
    else {
      croak $CLASS. "::search: Error...query ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error query undefined!\n";
  }

  if ( defined( $value = $this->getGapInit() )
       && $value =~ /\d+/ )
  {
    $p{'gapopen'} = abs( $value - $this->getInsGapExt() );
  }
  else {
    $p{'gapopen'} = 12;
  }

  if ( defined( $value = $this->getInsGapExt() )
       && $value =~ /\d+/ )
  {
    $p{'gapextend'} = abs( $value );
  }
  else {
    $p{'gapextend'} = 2;
  }

  if ( ( $value = $this->getMaskLevel() ) ) {
    $p{'mask_level'} = $value if ( $value > 0 );
  }

  if (    ( $value = $this->getScoreMode() )
       && ( $value == SearchEngineI::basicScoreMode ) )
  {

    # Do nothing
  }
  else {
    $p{'complexity_adjust'} = 1;
  }

  if ( defined( $value = $this->getMinMatch() )
       && $value =~ /\d+/ )
  {
    $p{'word_size'} = $value;
  }
  else {
    $p{'word_size'} = 14;
  }

  # Translate SearchEngine minScore/Bandwidth
  # to RMBlast's -xdrop_ungap/-xdrop_gap_final/
  # -xdrop_gap parameters.  NOTE: There is much
  # legacy support encoded here.
  if ( defined( $value = $this->getMinScore() )
       && $value =~ /\d+/ )
  {
    my $minScore = $value;

    ## There are three approaches encoded here
    ## that are selectable by using +bandwidth,
    ## 0 bandwidth, and -bandwidth.  An 'undefined'
    ## bandwidth defaults to the +bandwidth method.
    ##
    ##  +bandwidth: use legacy MaskerAid tranlsations
    ##              for minScore to RMBlast parameters.
    ##              [default]
    ##  0 bandwidth: unused? special case
    ##  -bandwidth: use bandwidth magnitude, in addition
    ##              to minScore to simulate a bandwidth
    ##              in RMBlast. This uses the gap penalties
    ##              to work out a correct -xdrop_gap_final
    ##              parameter.
    ## TODO: Document diff between xdrops for NCBI Blast
    ##    xdrop_ungap: 
    ##    xdrop_gap:
    ##    xdrop_gap_final:
    ##

    # NOTE: This was assigned to the refinement case in 4.1.7
    # and led to tons of refinement alignments.  Reverting to
    # the previous method ( bandwidth = "-1" ) for Repeatmasker
    # refinement.
    # NOTE: This is not equivalent
    # to "undefined".  It must have a value of "0".
    if ( $this->getBandwidth() eq "0" ) {
      $p{'xdrop_ungap'}     = $minScore * 2;
      $p{'xdrop_gap_final'} = $minScore * 4;
      $p{'xdrop_gap'}       = int( $minScore / 2 );
    }
    elsif ( defined( $value = $this->getBandwidth() )
            && $value < 0 )
    {
      # In crossmatch the bandwidth parameter is the off-diagonal
      # distance tolerated.  So the full band is (cm_bandwidth * 2 + 1)
      # wide.  Here we use bandwidth to indicate the full width of the band
      # (legacy...should have used the same definition as cm) so a bandwidth
      # of -29 is equivalent to the crossmatch bandwidth of 14.
      $p{'xdrop_ungap'} = $minScore * 2;

      # Ins/Del extension penalties are the same for RMBlast ( only cm differentiates )
      # The tolerated gapped xdrop should tolerate a gap of bandwidth #.  So
      # Gap init penalty + (extension penalty * bandwidth )
      $p{'xdrop_gap_final'} =
          ( abs( $value ) * abs( $this->getInsGapExt() ) ) +
          abs( $this->getGapInit() );
      $p{'xdrop_gap'} = int( $minScore / 2 );
    }
    else {
      # These are inherited from MaskerAid.  It's a strange choice as
      # it creates a side effect on indel size.  Lower minscores
      # reduce the allowable indel length whereas higher minscores
      # allow really large indel sizes ( assuming they also reach
      # the score threshold ).
      $p{'xdrop_ungap'}     = $minScore * 2;
      $p{'xdrop_gap_final'} = $minScore;
      $p{'xdrop_gap'}       = int( $minScore / 2 );
    }
    $p{'min_raw_gapped_score'} = $minScore;
    $p{'dust'}                 = "no";
  }

  if ( $this->hasTabOutput() ) {
    $p{'outfmt_fields'} = [ $this->getTabOutputFields() ];
  }

  #
  # TODO: Is there some way we can check to see if this
  #       is a MT version of rmblastn?  Also a good way
  #       to know if we should call with threads turned
  #       on?
  if ( defined( $value = $this->getCores() ) ) {
    $p{'num_threads'} = $value;
  }
  else {
    $p{'num_threads'} = 4;
  }

  if ( defined( $value = $this->getThreadByQuery() ) && $value > 0 ) {
    $p{'mt_mode'} = 1;
  }

  if ( defined( $value = $this->getMatrix() ) ) {

    # Test if matrix exists
    if ( -f $value ) {

      # NCBIBLAST requires that the matrix filename parameter
      # be relative to a directory path specified in
      # environment variables.
      my @path = split( /[\\\/]/, $value );
      $p{'matrix'}     = pop @path;
      $p{'matrix_dir'} = join( "/", @path );

    }
    else {
      croak $CLASS. "::search: Error...matrix ($value) does not exist!\n";
    }
  }

  return ( \%p );
}

##-------------------------------------------------------------------------##
##  Use: my $bool = $this->hasTabOutput();
##
##  Whether this engine should be asked for tab delimited output.
##-------------------------------------------------------------------------##
sub hasTabOutput {
  my $this = shift;

  return ( exists $this->{'hasTabFormat'} && !$this->{'forceLegacyParser'} );
}

##-------------------------------------------------------------------------##
##  Use: my @fields = $this->getTabOutputFields();
##
##  The tab delimited output columns, in the order parseTabOutput()
##  expects them.  qseq and sseq are added only when the caller wants
##  alignments.
##-------------------------------------------------------------------------##
sub getTabOutputFields {
  my $this = shift;

  my @fields = qw( score perc_sub perc_query_gap perc_db_gap qseqid qstart
      qend qlen sstrand sseqid sstart send slen kdiv cpg_kdiv transi transv
      cpg_sites );
  push @fields, qw( qseq sseq ) if ( $this->getGenerateAlignments() );

  return @fields;
}

##-------------------------------------------------------------------------##
##  Use: my $string = $this->_renderParameters( $paramsRef );
##
##  Render computed parameter values onto an NCBI toolkit ( 2.x ) command
##  line.  A subclass with a different option syntax overrides this one
##  method.
##-------------------------------------------------------------------------##
sub _renderParameters {
  my $this = shift;
  my $p    = shift;

  my $parameters = " -num_alignments " . $p->{'num_alignments'};
  $parameters .= " -db " . $p->{'db'};
  $parameters .= " -query " . $p->{'query'};
  $parameters .= " -gapopen " . $p->{'gapopen'};
  $parameters .= " -gapextend " . $p->{'gapextend'};
  $parameters .= " -mask_level " . $p->{'mask_level'}
      if ( defined $p->{'mask_level'} );
  $parameters .= " -complexity_adjust " if ( $p->{'complexity_adjust'} );
  $parameters .= " -word_size " . $p->{'word_size'};

  if ( defined $p->{'xdrop_ungap'} ) {
    $parameters .=
          " -xdrop_ungap "
        . $p->{'xdrop_ungap'}
        . " -xdrop_gap_final "
        . $p->{'xdrop_gap_final'}
        . " -xdrop_gap "
        . $p->{'xdrop_gap'} . " ";
  }
  $parameters .= " -min_raw_gapped_score "
      . $p->{'min_raw_gapped_score'} . " -dust "
      . $p->{'dust'} . " "
      if ( defined $p->{'min_raw_gapped_score'} );

  $parameters .=
      " -outfmt=\"6 " . join( " ", @{ $p->{'outfmt_fields'} } ) . "\" "
      if ( defined $p->{'outfmt_fields'} );

  $parameters .= " -num_threads " . $p->{'num_threads'} . " ";
  $parameters .= " -mt_mode " . $p->{'mt_mode'} . " "
      if ( defined $p->{'mt_mode'} );
  $parameters .= " -matrix " . $p->{'matrix'} if ( defined $p->{'matrix'} );

  return $parameters;
}

##-------------------------------------------------------------------------##

=head2 search()

  Use: my ( $resultCode, $SearchResultCollectionI ) = search( );
 or                                                                            
  Use: my ( $resultCode, $SearchResultCollectionI )
                          = search( matrix=>"7p16g.matrix",
                                    ...
                                  );

  Run the search and return a SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub search {
  my $this           = shift;
  my %nameValuePairs = @_;

  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak( $CLASS . "::search: Instance variable $name doesn't exist." );
      }
      $this->$method( $value );
    }
  }

  # Form the command line
  my $cmdLine = $this->getParameters();

  my $matrixName;
  if ( defined( my $value = $this->getMatrix() ) ) {

    # Test if matrix exists
    if ( -f $value ) {

      # NCBIBLAST requires that the matrix filename parameter
      # be relative to a directory path specified in
      # environment variables.
      my @path = split( /[\\\/]/, $value );
      $matrixName = pop @path;

      # Set the environment
      $ENV{BLASTMAT} = join( "/", @path );
      print "Setting BLASTMAT to: " . $ENV{BLASTMAT} . "\n"
          if ( $this->getDEBUG() );
    }
    else {
      croak $CLASS. "::search: Error...matrix ($value) does not exist!\n";
    }
  }

  my $outputDirName;
  if ( defined $this->getTempDir() && -d $this->getTempDir() ) {
    $outputDirName = $this->getTempDir();
  }
  else {
    $outputDirName = dirname( $this->getQuery() );
  }

  # Invoke engine and handle errors
  my $POUTPUT = new FileHandle;
  my $errFile;
  my $currentTime;
  my $rand = 0;
  do {
    $currentTime = time();
    $rand = rand(5000);
    $errFile     = $outputDirName . "/ncResults-$currentTime-$$-$rand.err";
  } while ( -f $errFile );
  my $outFile = $outputDirName . "/ncResults-$currentTime-$$-$rand.out";
  my $pid;

  print $CLASS
      . "::search() Invoking search engine as: $cmdLine "
      . " 2>$errFile |\n"
      if ( $this->getDEBUG() );

  # DISABLE: For timing only
  #my $t0 = gettimeofday( ); 
  $pid = open( $POUTPUT, "$cmdLine 2>$errFile |" );

  my %parseParams = ();
  $parseParams{'debug'} = $this->getDEBUG() if ( $this->getDEBUG );
  $parseParams{'excludeAlignments'} = 1 if ( !$this->getGenerateAlignments() );
  $parseParams{'matrixName'}        = $matrixName;
  $parseParams{'format'} = "tab" if ( $this->{'hasTabFormat'} && !$this->{'forceLegacyParser'} );

  # Places to hold result codes and alignment data
  my $resultCode;
  my $searchResultsCollection;

  # Passing $POUTPUT directly to parseOutput and  saving to 
  # a file first and then passing the file to parseOutput work
  # out to be about the same in overall timing.  For now
  # I am going to avoid saving it to disk if we don't need to.
  if ( $this->getDEBUG() ) {
    # TODO DEBUGGING
    #if ( $this->getDEBUG ) {
    #  system(   "cp "
    #          . $this->getQuery()
    #          . " $outputDirName"
    #          . "/before-$currentTime-$$.fa" );
    #}

    open OUT, ">$outFile";
    while ( <$POUTPUT> ) {
      print OUT $_;
    }
    close OUT;
    close $POUTPUT;
    $resultCode = ( $? >> 8 );

    # DISABLE: For timing only ( normally disabled )
    #my $t1 = gettimeofday( ); 
    #my $elapsed = $t1 - $t0; 
    #print "NCBIBlast runtime: $elapsed secs\n";

    $parseParams{'searchOutput'} = $outFile;
    $searchResultsCollection = $this->parseOutput( %parseParams );
  }else {
    $parseParams{'searchOutput'} = $POUTPUT;
    $searchResultsCollection = $this->parseOutput( %parseParams );
    close $POUTPUT;
    $resultCode = ( $? >> 8 );
  }

  print "NCBIBlast returned a the following result code >$resultCode<\n"
      if ( $this->getDEBUG() );


  #
  # Postprocess the results
  #
  if ( defined $searchResultsCollection
       && $searchResultsCollection->size() > 0 )
  {

    if ( ! $this->{'hasTabFormat'} || $this->{'forceLegacyParser'} )
    {
      ## For some reason when complexity adjustment is turned off
      ## RMBLAST isn't respecting the min_raw_gapped_score. For
      ## now I am doing it as a perl postprocessing step below.
      my $minScore = $this->getMinScore();
      if ( defined $minScore ) {
        print $CLASS
            . "::search: "
            . $searchResultsCollection->size()
            . " hits before minScore filtering\n"
            if ( $this->getDEBUG() );
  
        for ( my $i = $searchResultsCollection->size() - 1 ; $i >= 0 ; $i-- ) {
          if ( $searchResultsCollection->get( $i )->getScore() < $minScore ) {
            $searchResultsCollection->remove( $i );
          }
        }
        print $CLASS
            . "::search: "
            . $searchResultsCollection->size()
            . " hits after minScore filtering\n"
            if ( $this->getDEBUG() );
      }
  
      for ( my $i = 0 ; $i < $searchResultsCollection->size() ; $i++ ) {
        my $result = $searchResultsCollection->get( $i );
  
        #
        # Calculate Kimura divergence using the CpG modification described
        # in SearchResult.pm
        #
        my ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
            $result->calcKimuraDivergence( divCpGMod => 1 );
        $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
      }
    }
    # The final result collection should be sorted by
    #  queryname and secondarily by query start position.
    $searchResultsCollection->sort(
      sub ($$) {
        $_[ 0 ]->getQueryName cmp $_[ 1 ]->getQueryName()
            || $_[ 0 ]->getQueryStart() <=> $_[ 1 ]->getQueryStart();
      }
    );


  }

  if ( $this->getDEBUG() ) {
    return ( $resultCode, $searchResultsCollection, $outFile, $errFile );
  }else {
    unlink $errFile;
    return ( $resultCode, $searchResultsCollection, "", $errFile );
  }
}

##-------------------------------------------------------------------------##

=head1 Class Methods

=cut

##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 parseOutput()

  Use: my $SearchResultCollection = NCBIBLASTSearchEngine::parseOutput(
                                     searchOutput => $filename|$FH,
                                     [format => 'tab'],
                                     [matrixName => $matrixName],
                                     [excludeAlignments => 1]  );

  Parse the result of a search and return a SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub parseOutput {
  # Callable two ways: as a class function, which is the historical
  # interface that external tools use, and as an instance method.  Only
  # the method form can dispatch to a subclass's parser, so search() uses
  # that.
  my $this;
  $this = shift if ( ref( $_[ 0 ] ) && UNIVERSAL::isa( $_[ 0 ], $CLASS ) );
  my %nameValueParams = @_;

  croak $CLASS. "::parseOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  if ( exists $nameValueParams{'format'} &&  $nameValueParams{'format'} eq "tab" )
  {
    return ( $this
             ? $this->parseTabOutput( %nameValueParams )
             : &parseTabOutput( %nameValueParams ) );
  }else {
    return ( $this
             ? $this->parseReportOutput( %nameValueParams )
             : &parseReportOutput( %nameValueParams ) );
  }
}


##-------------------------------------------------------------------------##

=head2 parseTabOutput()

  Use: my $SearchResultCollection = NCBIBLASTSearchEngine::parseTabOutput(
                                     searchOutput => $filename|$FH,
                                     [matrixName => $matrixName],
                                     [excludeAlignments => 1]  );

  Parse the result of a search in NCBI Blast tab delemited format
  (2.13.0+ and up) and return a SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub parseTabOutput {
  my $this;
  $this = shift if ( ref( $_[ 0 ] ) && UNIVERSAL::isa( $_[ 0 ], $CLASS ) );
  my %nameValueParams = @_;

  croak $CLASS. "::parseTabOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  my $NCBIFILE;
  if ( ref( $nameValueParams{'searchOutput'} ) !~ /GLOB|FileHandle|IO::File/ ) {
    print $CLASS
        . "::parseTabOutput() Opening file "
        . $nameValueParams{'searchOutput'} . "\n"
        if ( $nameValueParams{'debug'} );
    open $NCBIFILE, $nameValueParams{'searchOutput'}
        or croak $CLASS
        . "::parseTabOutput: Unable to open "
        . "results file: $nameValueParams{'searchOutput'} : $!";
  }
  else {
    $NCBIFILE = $nameValueParams{'searchOutput'};
  }

  my $callbackFunc = undef;
  if ( defined $nameValueParams{'callback'}
       && ref( $nameValueParams{'callback'} ) == /CODE/ )
  {
    $callbackFunc = $nameValueParams{'callback'};
  }

  my $matrix;
  $matrix = $nameValueParams{'matrixName'}
      if ( defined $nameValueParams{'matrixName'} );

  my $resultColl = SearchResultCollection->new();

  # So simple it hurts, no parsing madness with tab separated value
  while ( <$NCBIFILE> ) {

    print "RMBLASTN: $_"
        if ( exists $nameValueParams{'debug'}
             && $nameValueParams{'debug'} > 8 );

    if ( /^\d+/ ) {
      s/[\n\r]+//g;
      my @flds = split(/\t/);
      #
      # 18 Field Format:
      #   -outfmt="6 score perc_sub perc_query_gap perc_db_gap qseqid qstart 
      #            qend qlen sstrand sseqid sstart send slen kdiv cpg_kdiv
      #            transi transv cpg_sites"
      # 20 Field Format:
      #   -outfmt="6 score perc_sub perc_query_gap perc_db_gap qseqid qstart 
      #            qend qlen sstrand sseqid sstart send slen kdiv cpg_kdiv 
      #            transi transv cpg_sites qseq sseq"
      #
      # Example:
      #  2430    17.63   8.33    6.51    qseq1  595     1530    1587    plus    
      #     dseq1    1722      2673    6004    13.21   12.91   73      9       
      #     91      AATTGTCACCAAACAAAT	AATTGTTACTAAATTTGTCA-ACAAAT
      #
      if ( @flds == 18 || @flds == 20 ){
         my $orient = "";
         my $sbjStart = $flds[10];
         my $sbjEnd = $flds[11];
         if ( $flds[8] eq "minus" )
         {
           $orient = "C";
           $sbjStart = $flds[11];
           $sbjEnd = $flds[10];
         }
         my $result = SearchResult->new(
                                     queryName      => $flds[4],
                                     queryStart     => $flds[5],
                                     queryEnd       => $flds[6],
                                     queryRemaining => ( $flds[7] - $flds[6] ),
                                     orientation    => $orient,
                                     subjName       => $flds[9],
                                     subjStart      => $sbjStart,
                                     subjEnd        => $sbjEnd,
                                     subjRemaining  => ( $flds[12] - $sbjEnd ),
                                     pctDiverge     => $flds[1],
                                     pctInsert      => $flds[3],
                                     pctDelete      => $flds[2],
                                     matrixName     => $matrix,
                                     score          => $flds[0],
                                     pctRawKimuraDiverge => $flds[13],
                                     pctKimuraDiverge    => $flds[14],  # CpG adjusted
                                     cpGSites            => $flds[17]    
        );
        if ( @flds == 20 ) {
          $result->setQueryString($flds[18]); 
          $result->setSubjString($flds[19]); 
        }

        if ( defined $callbackFunc ) {
          $callbackFunc->( $result );
        }
        else {
          $resultColl->add( $result );
        }


      }else {
        # Warn....bad number of fields
      }
      
    }else {
      # possibly a header row?
    }

# Reorient if this is a reverse strand
# hit.
#
#if ( $qryOrient eq "C" ) 
#
#          # Fix the sequence orientation so that it
#          # matches the SearchResult.pm convention of
#          # the query being in the forward direction.
#          $qrySeq = reverse $qrySeq;
#          $qrySeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
#          $sbjSeq = reverse $sbjSeq;
#          $sbjSeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
#          $orientation = "C";

  }
  close $NCBIFILE;

  return $resultColl;
}


##-------------------------------------------------------------------------##

=head2 parseReportOutput()

  Use: my $SearchResultCollection = NCBIBLASTSearchEngine::parseReportOutput(
                                     searchOutput => $filename|$FH,
                                     [matrixName => $matrixName],
                                     [excludeAlignments => 1]  );

  Parse the result of a search in NCBI Blast report alignment format (default)
  and return a SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub parseReportOutput {
  my $this;
  $this = shift if ( ref( $_[ 0 ] ) && UNIVERSAL::isa( $_[ 0 ], $CLASS ) );
  my %nameValueParams = @_;

  croak $CLASS. "::parseOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  my $NCBIFILE;
  if ( ref( $nameValueParams{'searchOutput'} ) !~ /GLOB|FileHandle|IO::File/ ) {
    print $CLASS
        . "::parseOutput() Opening file "
        . $nameValueParams{'searchOutput'} . "\n"
        if ( $nameValueParams{'debug'} );
    open $NCBIFILE, $nameValueParams{'searchOutput'}
        or croak $CLASS
        . "::parseOutput: Unable to open "
        . "results file: $nameValueParams{'searchOutput'} : $!";
  }
  else {
    $NCBIFILE = $nameValueParams{'searchOutput'};
  }

  my $callbackFunc = undef;
  if ( defined $nameValueParams{'callback'}
       && ref( $nameValueParams{'callback'} ) == /CODE/ )
  {
    $callbackFunc = $nameValueParams{'callback'};
  }

  my $inAlignState = 0;
  my $sbjID        = "";
  my $qryID        = "";
  my $absIndex     = 0;
  my $score        = 0;
  my $adjScore     = 0;
  my $sbjSeq       = "";
  my $qrySeq       = "";
  my $sbjOrient    = "";
  my $qryOrient    = "";    # TODO: Check to see if this ever happens
  my $sbjStart     = 0;
  my $sbjEnd       = 0;
  my $qryStart     = 0;
  my $qryEnd       = 0;
  my $qryLength    = 0;
  my $sbjLength    = 0;
  my $matrix;
  $matrix = $nameValueParams{'matrixName'}
      if ( defined $nameValueParams{'matrixName'} );

  my $resultColl = SearchResultCollection->new();

  my %gapAndX = (
    '-' => 1,
    'x' => 1,
    'X' => 1 );

  my %IUBMatchLookup = (
    "AA" => 1,
    "AC" => 0,
    "AG" => 0,
    "AT" => 0,
    "AB" => 0,
    "AD" => 1,
    "AH" => 1,
    "AV" => 1,
    "AR" => 1,
    "AY" => 0,
    "AK" => 0,
    "AM" => 1,
    "AS" => 0,
    "AW" => 1,
    "AN" => 1,

    "CA" => 0,
    "CC" => 1,
    "CG" => 0,
    "CT" => 0,
    "CB" => 1,
    "CD" => 0,
    "CH" => 1,
    "CV" => 1,
    "CR" => 0,
    "CY" => 1,
    "CK" => 0,
    "CM" => 1,
    "CS" => 1,
    "CW" => 0,
    "CN" => 1,

    "GA" => 0,
    "GC" => 0,
    "GG" => 1,
    "GT" => 0,
    "GB" => 1,
    "GD" => 1,
    "GH" => 0,
    "GV" => 1,
    "GR" => 1,
    "GY" => 0,
    "GK" => 1,
    "GM" => 0,
    "GS" => 1,
    "GW" => 0,
    "GN" => 1,

    "TA" => 0,
    "TC" => 0,
    "TG" => 0,
    "TT" => 1,
    "TB" => 1,
    "TD" => 1,
    "TH" => 1,
    "TV" => 0,
    "TR" => 0,
    "TY" => 1,
    "TK" => 1,
    "TM" => 0,
    "TS" => 0,
    "TW" => 1,
    "TN" => 1
  );

  while ( <$NCBIFILE> ) {

    print "RMBLASTN: $_"
        if ( exists $nameValueParams{'debug'}
             && $nameValueParams{'debug'} > 8 );

    #
    # Conditions for the end of a hit record:
    #   o Must have seen a score
    #   o Must be in the alignment state (ie. have seen Query: and
    #     Subj: recently)
    #   o Must see something which isn't either "Query:" "Subj:" or " "
    #     or must see the end of the file
    #
    if ( $inAlignState ) {
      if ( !/^(Query |Sbjct |\s{8}|\n|\r)/ || eof ) {

        #
        # Reorient if this is a reverse strand
        # hit.
        #
        my $orientation = "";
        if ( $qryOrient eq "C" ) {

          # Fix the sequence orientation so that it
          # matches the SearchResult.pm convention of
          # the query being in the forward direction.
          $qrySeq = reverse $qrySeq;
          $qrySeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
          $sbjSeq = reverse $sbjSeq;
          $sbjSeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
          $orientation = "C";
        }
        elsif ( $sbjOrient eq "C" ) {
          $orientation = "C";
        }

        #
        # Calculate percent divergence
        #           percent insertions
        #           percent deletions
        #
        my %baseFreq = ();
        my $mismatch = 0;
        for ( my $i = 0 ; $i < length( $qrySeq ) ; $i++ ) {
          my $qryBase = substr( $qrySeq, $i, 1 );
          my $sbjBase = substr( $sbjSeq, $i, 1 );
          # NYTProfiler found this to be extremely expensive!
          #next if ( $qryBase =~ /-|x/i || $sbjBase =~ /-|x/i );
          # This is faster...but
          #next if ( $qryBase eq '-' || $sbjBase eq '-' ||
          #          $qryBase eq'x' || $sbjBase eq 'x' ||
          #          $qryBase eq 'X' || $sbjBase eq 'X' );
          # Faster still:
          next if ( exists $gapAndX{$qryBase} || exists $gapAndX{$sbjBase} );
          $baseFreq{$qryBase}++;
          $mismatch++
              if (    !$IUBMatchLookup{ uc( $qryBase . $sbjBase ) }
                   && !$IUBMatchLookup{ uc( $sbjBase . $qryBase ) } );

        }
        my $percDiv =
            sprintf( "%4.2f", $mismatch * 100 / ( $qryEnd + 1 - $qryStart ) );
        my $qgap = $qrySeq =~ tr/-/-/;
        my $sgap = $sbjSeq =~ tr/-/-/;
        my $percIns =
            sprintf( "%4.2f", $sgap * 100 / ( ( $sbjEnd + 1 ) - $sbjStart ) );
        my $percDel =
            sprintf( "%4.2f", $qgap * 100 / ( ( $qryEnd + 1 ) - $qryStart ) );

        if ( defined( $nameValueParams{'excludeAlignments'} ) ) {
          $qrySeq = "";
          $sbjSeq = "";
        }

        my $result = SearchResult->new(
                                     queryName      => $qryID,
                                     queryStart     => $qryStart,
                                     queryEnd       => $qryEnd,
                                     queryRemaining => ( $qryLength - $qryEnd ),
                                     queryString    => $qrySeq,
                                     subjString     => $sbjSeq,
                                     orientation    => $orientation,
                                     subjName       => $sbjID,
                                     subjStart      => $sbjStart,
                                     subjEnd        => $sbjEnd,
                                     subjRemaining  => ( $sbjLength - $sbjEnd ),
                                     pctDiverge     => $percDiv,
                                     pctInsert      => $percIns,
                                     pctDelete      => $percDel,
                                     matrixName     => $matrix,
                                     score          => $score
        );

        if ( defined $callbackFunc ) {
          $callbackFunc->( $result );
        }
        else {
          $resultColl->add( $result );
        }

        $score        = "";
        $sbjSeq       = "";
        $qrySeq       = "";
        $qryOrient    = "";
        $sbjOrient    = "";
        $sbjStart     = 0;
        $sbjEnd       = 0;
        $qryStart     = 0;
        $qryEnd       = 0;
        $inAlignState = 0;
      }
    }

    #
    # Query alignment
    #
    if ( /^Query\s+(\d+)\s+(\S+)\s+(\d+)/ ) {
      $qrySeq .= uc( $2 );
      if ( $qryStart < 1 ) {
        $qryStart = _min( $1, $3 );
        $qryEnd   = _max( $1, $3 );
      }
      else {
        $qryStart = _min( _min( $1, $3 ), $qryStart );
        $qryEnd   = _max( _max( $1, $3 ), $qryEnd );
      }
      $qryOrient = "C" if ( $1 > $3 );
      $inAlignState = 1;
    }
    elsif ( /^Query\s+(\-+)\s*/ ) {

      # I have seen cases in RMBLASTN output where
      # coordinates are not given if the entire line
      # is only gap characters.
      $qrySeq .= uc( $1 );
    }

    #
    # Subject alignment
    #
    if ( /^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)/ ) {
      $sbjSeq .= uc( $2 );
      my $leftNum  = $1;
      my $rightNum = $3;
      if ( $sbjStart < 1 ) {
        $sbjStart = _min( $leftNum, $rightNum );
        $sbjEnd   = _max( $leftNum, $rightNum );
      }
      else {
        $sbjStart = _min( _min( $leftNum, $rightNum ), $sbjStart );
        $sbjEnd   = _max( _max( $leftNum, $rightNum ), $sbjEnd );
      }
      $inAlignState = 1;
    }
    elsif ( /^Sbjct\s+(\-+)\s*/ ) {

      # I have seen cases in RMBLASTN output where
      # coordinates are not given if the entire line
      # is only gap characters.
      $sbjSeq .= uc( $1 );
    }

    #
    # Query source name
    #
    if ( /^Query\s*=\s*(\S+)/ ) {
      $qryID     = $1;
      $qryLength = -1;
    }

    #
    # Length of qyery/database
    #
    if ( /^Length\s*=\s*(\d+)\s*$/ ) {
      if ( $qryLength > 0 ) {
        $sbjLength = $1;
      }
      else {
        $qryLength = $1;
      }
    }

    #
    # Score
    #
    if ( /Score\s*=\s*(\d+)/ ) {
      $score = $1;
    }

    #
    # Strand
    #
    if ( /Strand\s*=\s*(Plus|Minus)\/(Plus|Minus)/ ) {
      if ( $2 eq "Minus" ) {
        $sbjOrient = "C";
      }
    }

    #
    # Hit description line
    #
    if ( /^>\s*(\S+).*/ ) {
      $sbjID = $1;
    }

  }
  close $NCBIFILE;

  return $resultColl;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _min( $num1, $num2 );
##
##              $num1   :       A number to be compared
##              $num2   :       A number to be comprared
##
##      Returns:                The minimum of the two numbers
##
##-------------------------------------------------------------------------##
sub _min {
  my ( $num1, $num2 ) = @_;
  if ( $num1 < $num2 ) {
    return ( $num1 );
  }
  else {
    return ( $num2 );
  }
}

##-------------------------------------------------------------------------##
## Use: my _max( $num1, $num2 );
##
##              $num1   :       A number to be compared
##              $num2   :       A number to be comprared
##
##      Returns:                The maximum of the two numbers
##
##-------------------------------------------------------------------------##
sub _max {
  my ( $num1, $num2 ) = @_;
  if ( $num1 < $num2 ) {
    return ( $num2 );
  }
  else {
    return ( $num1 );
  }
}

##-------------------------------------------------------------------------##
## Use: my _ucFirst( $string );
##
##   Uppercases the first character in a string and returns it.
##
##-------------------------------------------------------------------------##
sub _ucFirst {
  my $string = shift;

  if ( defined $string && $string ne "" ) {
    substr( $string, 0, 1 ) = uc( substr( $string, 0, 1 ) );
  }
  return $string;
}

##-------------------------------------------------------------------------##
## Serialization & Debug Routines
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $string = toString([$this]);
##
##      $this         : Normally passed implicitly
##
##  Returns
##
##      Uses the Data::Dumper to create a printable reprentation
##      of a data structure.  In this case the object data itself.
##
##-------------------------------------------------------------------------##
sub toString {
  my $this = shift;
  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  return $data_dumper->Dump();
}

##-------------------------------------------------------------------------##
## Use: my serializeOUT( $filename );
##
##	  $filename	: A filename to be created
##
##  Returns
##
##	Uses the Data::Dumper module to save out the data
##	structure as a text file.  This text file can be
##	read back into an object of this type.
##
##-------------------------------------------------------------------------##
sub serializeOUT {
  my $this     = shift;
  my $fileName = shift;

  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  open OUT, ">$fileName";
  print OUT $data_dumper->Dump();
  close OUT;
}

##-------------------------------------------------------------------------##
## Use: my serializeIN( $filename );
##
##	$filename	: A filename containing a serialized object
##
##  Returns
##
##	Uses the Data::Dumper module to read in data
##	from a serialized PERL object or data structure.
##
##-------------------------------------------------------------------------##
sub serializeIN {
  my $this         = shift;
  my $fileName     = shift;
  my $fileContents = "";
  my $oldSep       = $/;
  undef $/;
  my $in;
  open $in, "$fileName";
  $fileContents = <$in>;
  $/            = $oldSep;
  close $in;
  return eval( $fileContents );
}

1;
