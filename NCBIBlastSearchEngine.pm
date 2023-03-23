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

  my $result = `$value -version 2>&1`;
  if ( $result =~ /rmblast[n]\s*:?\s*(\S.*)/ ) {
    $this->{'engineName'} = "rmblastn";
    $this->{'version'}    = $1;
    if ( $this->{'version'} =~ /(\d+)\.(\d+)\.(\d+)\+/ ) {
      my $majorVer = $1;
      my $minorVer = $2;
      my $revision = $3;
      if ( $majorVer > 2 || ($majorVer == 2 && $minorVer >= 13 )) {
         # Since the 2.13.0+ release of RMBlast we now have:
         #    - pre-computed Kimura divergence, Kimura CpG adjusted,
         #      transitions, transversions, cpg_sites, and the cross_match
         #      stats (perc_sub, perc_query_gap, perc_subj_gap).
         #    - Ability to thread on the query sequences rather than just
         #      the subject sequences.
         #    - The ability to output tab delimited format with all the above
         #      fields.
         $this->{'hasQueryThreading'} = 1;
         $this->{'hasTabFormat'} = 1;
      }
    }
  }
  else {
    croak $CLASS
        . "::setPathToEngine( $value ): Cannot determine "
        . "engine variant and version!\n";
  }

  my $oldValue = $this->{'pathToEngine'};
  $this->{'pathToEngine'} = $value;

  return $oldValue;
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

  # Generate parameter line
  my $parameters    = " -num_alignments 9999999";
  my $spanParameter = "";
  my $value;
  if ( ( $value = $this->getSubject() ) ) {

    # Make sure we have the compressed form of the database handy
    if (    -f "$value.nin"
         || -f "$value.nhr"
         || -f "$value.nsq" )
    {
      $parameters .= " -db $value";
    }
    else {
      croak $CLASS
          . "::search: Error...compressed subject "
          . "database ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error subject undefined!\n";
  }

  if ( ( $value = $this->getQuery() ) ) {
    if ( -f $value ) {
      $parameters .= " -query $value";
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
    $parameters .= " -gapopen " . abs( $value - $this->getInsGapExt() );
  }
  else {
    $parameters .= " -gapopen 12";
  }

  if ( defined( $value = $this->getInsGapExt() )
       && $value =~ /\d+/ )
  {
    $parameters .= " -gapextend " . abs( $value );
  }
  else {
    $parameters .= " -gapextend 2";
  }
  if ( ( $value = $this->getMaskLevel() ) ) {
    $parameters .= " -mask_level $value" if ( $value > 0 );
  }
  if (    ( $value = $this->getScoreMode() )
       && ( $value == SearchEngineI::basicScoreMode ) )
  {

    # Do nothing
  }
  else {
    $parameters .= " -complexity_adjust ";
  }

  if ( defined( $value = $this->getMinMatch() )
       && $value =~ /\d+/ )
  {
    $parameters .= " -word_size $value";
  }
  else {
    $parameters .= " -word_size 14";
  }

  if ( defined( $value = $this->getMinScore() )
       && $value =~ /\d+/ )
  {
    my $minScore = $value;
    ## In some cases we want to generate a pseudo
    ## global alignment.  I.e when refining a previously
    ## aligned region with new consensi.
    ## Currently this special case is only used by
    ## the "align.pl" utility.  NOTE: To be safe, the
    ## only way this is invoked is with a negative value
    ## e.g -bandwidth -29 ( cm bandwidth default is 14 [14*2+1])
    if ( defined( $value = $this->getBandwidth() )
         && $value < 0 )
    {
      # Document diff between xdrops for NCBI Blast
      # xdrop_ungap: 
      # xdrop_gap:
      # xdrop_gap_final:
      $parameters .= " -xdrop_ungap " . ( $minScore * 2 ) . " -xdrop_gap_final "
          # Ins/Del extension penalties are the same for RMBlast ( only cm differentiates )
          # The tolerated gapped xdrop should tolerate a gap of bandwidth #.  So 
          # Gap init penalty + (extension penalty * bandwidth )
          . ( ( abs( $value ) * abs( $this->getInsGapExt() ) ) +
              abs( $this->getGapInit() ) )
          . " -xdrop_gap "
          . int( $minScore / 2 ) . " ";
    }
    else {

      # These are inherited from MaskerAid.  It's a strange choice as
      # it creates a side effect on indel size.  Lower minscores
      # reduce the allowable indel length whereas higher minscores
      # allow really large indel sizes ( assuming they also reach
      # the score threshold ).
      $parameters .=
            " -xdrop_ungap "
          . ( $minScore * 2 )
          . " -xdrop_gap_final "
          . ( $minScore )
          . " -xdrop_gap "
          . int( $minScore / 2 ) . " ";
    }
    $parameters .= " -min_raw_gapped_score $minScore -dust no ";
  }

  if ( exists $this->{'hasTabFormat'} && !$this->{'forceLegacyParser'} ) 
  {
    if ( $this->getGenerateAlignments() )
    {
      $parameters .= " -outfmt=\"6 score perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen sstrand sseqid sstart send slen kdiv cpg_kdiv transi transv cpg_sites qseq sseq\" ";
    }else {
      $parameters .= " -outfmt=\"6 score perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen sstrand sseqid sstart send slen kdiv cpg_kdiv transi transv cpg_sites\" ";
    }
  }

  #
  # TODO: Is there some way we can check to see if this
  #       is a MT version of rmblastn?  Also a good way
  #       to know if we should call with threads turned
  #       on?
  if ( defined( $value = $this->getCores() ) ) {
    $parameters .= " -num_threads $value ";
  }
  else {
    $parameters .= " -num_threads 4 ";
  }

  if ( defined( $value = $this->getThreadByQuery() ) && $value > 0 ) {
    $parameters .= " -mt_mode 1 ";
  }

  if ( defined( $value = $this->getMatrix() ) ) {

    # Test if matrix exists
    if ( -f $value ) {

      # NCBIBLAST requires that the matrix filename parameter
      # be relative to a directory path specified in
      # environment variables.
      my @path = split( /[\\\/]/, $value );
      my $matrix = pop @path;
      $parameters .= " -matrix $matrix";

    }
    else {
      croak $CLASS. "::search: Error...matrix ($value) does not exist!\n";
    }
  }

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
    $searchResultsCollection = parseOutput( %parseParams );
  }else {
    $parseParams{'searchOutput'} = $POUTPUT;
    $searchResultsCollection = parseOutput( %parseParams );
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
  my %nameValueParams = @_;

  croak $CLASS. "::parseOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  if ( exists $nameValueParams{'format'} &&  $nameValueParams{'format'} eq "tab" )
  {
    return ( &parseTabOutput(%nameValueParams) );
  }else {
    return ( &parseReportOutput(%nameValueParams) );
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
