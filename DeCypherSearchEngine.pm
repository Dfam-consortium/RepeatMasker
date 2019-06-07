#!/usr/bin/env perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) DeCypherSearchEngine.pm
##  Author:
##      Michael Sievers sievers@timelogic.com
##  Description:
##      An implementation of SearchEngineI for the
##      the TimeLogic DeCypher search engine.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2004 Developed by
#* Arian Smit and Robert Hubley.
#* Modified by Michael Sievers for TimeLogic 2005
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
# bless(
#      'DeCypherSearchEngine' );
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

DeCypherSearchEngine

=head1 SYNOPSIS

use DeCypherSearchEngine

Usage: 

  use SearchEngineI;
  use DeCypherSearchEngine;
  use SearchResultCollection;

  my $DeCypherEngine = DeCypherSearchEngine->new( 
                    pathToEngine=>"$DECYPHER" );

  $DeCypherEngine->setMatrix( "/users/bob/simple.matrix" );
  $DeCypherEngine->setQuery( "/users/bob/query.fasta" );
  $DeCypherEngine->setSubject( "/users/bob/subject.fasta" );
  my $searchResults = $DeCypherEngine->search();

=head1 DESCRIPTION

  A concrete implementation of the abstract class / interface SearchEngineI
  which use the TimeLogic DeCypher sequence search engine.

=head1 INSTANCE METHODS

=cut 

package DeCypherSearchEngine;
use strict;
use SearchEngineI;
use SearchResultCollection;
use Data::Dumper;
use FileHandle;
use File::Basename;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter SearchEngineI);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Version
#
my $VERSION  = 0.1;
my $CLASS    = "DeCypherSearchEngine";
my $DEBUG    = 0;
my $DECYPHER = "/cygdrive/c/dc_local/bin/";

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  croak $CLASS
      . "::new: Missing path to search engine!\n\n"
      . "use \$searchEngine = $CLASS->new( pathToEngine=>\"$DECYPHER/"
      . "dc_template\")\n"
      if ( not defined $nameValuePairs{'pathToEngine'} );

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  $this->setPathToEngine( $nameValuePairs{'pathToEngine'} );

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

  croak $CLASS. "::setPathToEngine( $value ): Program does not exist!"
      if ( not -x $value || `which $value` );

  my $result = `$value 2>&1`;
  if ( $result =~ /^BLASTP (\S+) \[.*/ ) {
    $this->{'version'} = $1;
  }

  my $oldValue = $this->{'pathToEngine'};
  $this->{'pathToEngine'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## Instance Methods
##-------------------------------------------------------------------------##

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

  # Test if engine is available
  my $engine = $this->getPathToEngine();
  if ( !defined $engine || !-f "$engine" ) {
    croak $CLASS
        . "::search: The path to the search engine is undefined or\n"
        . "is set incorrectly: $engine\n";
  }

  # Generate parameter line
  my $parameters = "";
  my $value;
  my $targStr = "";
  my $xtnd;

  if ( ( $value = $this->getSubject() ) ) {
    $value =~ /(.*\/)*(.*)/;    #Get the path-free target name
    $value = $2;
    if ( $value =~ /simple.lib/ || $value =~ /at.lib/ ) {
      $parameters = "-span1 t ";
    }

    open TARG, "dc_target -aa |";
    while ( $targStr = <TARG> ) {

      # Look for database name in target listing
      if ( $targStr =~ /$value/ ) {
        $parameters .= "-target $value ";
        last;
      }
    }
    if ( $parameters !~ /-target/ ) {
      croak $CLASS. "::search: Error...subject ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error subject undefined!\n";
  }

  my $outputDirName;
  if ( ( $value = $this->getQuery() ) ) {
    $value =~ s/\/cygdrive\/(\w)/$1:/;    # cygwin fix
    if ( -f $value ) {
      $parameters .= "-query $value -template RM-blastp -quiet ";
    }
    else {
      croak $CLASS. "::search: Error...query ($value) does not exist!\n";
    }
    $outputDirName = dirname( $value );
  }
  else {
    croak $CLASS. "::search: Error query undefined!\n";
  }

  if ( defined( $value = $this->getInsGapExt() )
       && $value =~ /\d+/ )
  {
    $xtnd = abs( $value );
    $parameters .= "-xtnd $xtnd ";
  }

  if ( defined( $value = $this->getGapInit() )
       && $value =~ /\d+/ )
  {
    $parameters .= "-open "
        . ( abs( $value ) - $xtnd )
        . " ";    # need this to correct open penalty to match WU behavior
  }

  if ( defined( $value = $this->getMinMatch() )
       && $value =~ /\d+/ )
  {
    $parameters .= "-word $value ";
  }
  else {
    $parameters .= "-word 14 ";
  }

  if ( defined( $value = $this->getMinScore() )
       && $value =~ /\d+/ )
  {
    $parameters .= "-thresh \"rawscore=$value\" ";
    $parameters .=
        "-xthresh " . int( $value * .45 ) . " ";    # extension threshold
    $parameters .= "-X $value ";                    # X dropoff
  }

  my ( $matrixRef, $matrixAlphabet, $matrixFreqRef ) = ();

  if ( defined( $value = $this->getMatrix() ) ) {

    # Test if matrix exists
    if ( -f $value ) {

      # Read in the matrix if we need to
      # adjust scores
      if ( $this->getScoreMode() == SearchEngineI::complexityAdjustedScoreMode )
      {
        ( $matrixRef, $matrixAlphabet, $matrixFreqRef ) = _readMatrix( $value );
      }
      $value =~ s/.*\/(.*)\.matrix$/$1/;    # need only the root name
      $parameters .= "-matrix %MATRIXPATH%/RM/$value ";
    }
    else {
      croak $CLASS. "::search: Error...matrix ($value) does not exist!\n";
    }
  }

  # Invoke engine and handle errors
  my $POUTPUT = new FileHandle;
  my $errFile;
  do {
    $errFile = $outputDirName . "/deResults-" . time() . ".err";
  } while ( -f $errFile );
  my $pid;
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

  print $CLASS
      . "::search() Invoking search engine as: $engine "
      . "$runParameters 2>$errFile |\n"
      if ( $DEBUG );

  $pid = open( $POUTPUT, "$engine $runParameters 2>$errFile |" );

  # Create SearchResultCollection object from
  # the engine results.
  my $searchResultsCollection;
  if ( $DEBUG ) {
    ## Create a debug file
    my $outFile;
    do {
      $outFile = $outputDirName . "/deResults-" . time() . ".out";
    } while ( -f $outFile );
    open OUT, ">$outFile";
    while ( <$POUTPUT> ) {
      print OUT $_;
    }
    close OUT;
    close $POUTPUT;
    $searchResultsCollection = parseOutput( searchOutput => $outFile );
  }
  else {

    # Just pipe to the parser
    $searchResultsCollection = parseOutput( searchOutput => $POUTPUT );
  }
  close $POUTPUT;

  my $resultCode = ( $? >> 8 );
  print "DeCypher returned a the following result code >$resultCode<\n"
      if ( $DEBUG );

  # Check return error code
  if ( $resultCode != 0 ) {
    open ERR, "<$errFile";
    my $errOk = 1;
    while ( <ERR> ) {
      if ( /^FATAL:/ ) {
        $errOk = 0 unless ( /context|cpus|P option|uses P|shorter/i );
      }
    }
    close ERR;
    $resultCode = 0 if ( $errOk == 1 );
  }
  unlink $errFile unless ( $DEBUG );

  #
  # Postprocess the results
  #
  if ( defined $searchResultsCollection
       && $searchResultsCollection->size() > 0 )
  {

    #
    # Calculate the complexity adjusted score if necessary
    #
    my $minScore = $this->getMinScore();
    if ( $this->getScoreMode() == SearchEngineI::complexityAdjustedScoreMode ) {

     # We need a defined matrix in order to calculate the complexity adjusted
     # score.  If the user overrides the parameters then we cannot be sure which
     # matrix was used.
      if ( !defined $matrixRef ) {
        croak "$CLASS::search: Evidently you overrode the parameters to "
            . "this SearchEngine but tried to use complexityAdjustedScoreMode "
            . "bad bad bad!\n";
      }

      my $lambda = _calculateLambda( $matrixRef, $matrixFreqRef );

      for ( my $i = $searchResultsCollection->size() - 1 ; $i >= 0 ; $i-- ) {
        my $adjScore = _complexityAdjust(
                          $searchResultsCollection->get( $i )->getScore(),
                          $searchResultsCollection->get( $i )->getQueryString(),
                          $searchResultsCollection->get( $i )->getSubjString(),
                          $lambda,
                          $matrixAlphabet,
                          $matrixRef,
                          $matrixFreqRef
        );

        if ( defined $minScore && $adjScore < $minScore ) {
          $searchResultsCollection->remove( $i );
        }
        else {
          $searchResultsCollection->get( $i )->setScore( $adjScore );
        }
      }
    }

    #
    # Filter set using "masklevel" concept
    #

    # Sort scores high to low
    $searchResultsCollection->sort(
      sub ($$) {
        $_[ 1 ]->getScore() <=> $_[ 0 ]->getScore();
      }
    );
    my $maskLevel;
    if ( defined( $maskLevel = $this->getMaskLevel() )
         && $maskLevel < 101 )
    {
      my %deleteHash = ();    # Hash to hold indexes of filtered entries
      for ( my $i = 0 ; $i < $searchResultsCollection->size() ; $i++ ) {
        for ( my $j = $i + 1 ; $j < $searchResultsCollection->size() ; $j++ ) {
          my $overlap = 0;
          my $perc    = 0;
          my $result1 = $searchResultsCollection->get( $i );
          my $result2 = $searchResultsCollection->get( $j );
          next if ( $result1->getQueryName() ne $result2->getQueryName() );

          # Get members
          my $begin1 = $result1->getQueryStart();
          my $begin2 = $result2->getQueryStart();
          my $end1   = $result1->getQueryEnd();
          my $end2   = $result2->getQueryEnd();

          # Check if they overlap
          next if ( $begin2 > $end1 || $begin1 > $end2 );

          # Calc overlap
          $overlap = $begin1 - $begin2 if ( $begin2 < $begin1 );
          $overlap += $end2 - $end1 if ( $end2 > $end1 );
          $perc = ( $overlap / ( $end2 - $begin2 + 1 ) ) * 100 if ( $overlap );
          if ( $perc < ( 100 - $maskLevel ) ) {
            $deleteHash{$j} = 1;
          }
        }
      }

      # Remove all hits which were filtered above
      if ( keys( %deleteHash ) ) {
        foreach my $index ( sort { $b <=> $a } keys( %deleteHash ) ) {
          $searchResultsCollection->remove( $index );
        }
      }
    }

    # Secondary sort by query position
    $searchResultsCollection->sort(
      sub ($$) {
        $_[ 0 ]->getQueryName cmp $_[ 1 ]->getQueryName()
            || $_[ 0 ]->getQueryStart() <=> $_[ 1 ]->getQueryStart();
      }
    );

  }

  return ( $resultCode, $searchResultsCollection );
}

##-------------------------------------------------------------------------##

=head1 Class Methods

=cut

##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 parseOutput()

  Use: my $SearchResultCollection = DeCypherSearchEngine::parseOutput(
                                     searchOutput => $filename|$FH,
                                     [excludeAlignments => 1],
                                     [scoreHighThresh => #],
                                     [scoreLowThresh => #],
                                     [subjPattern => ""],
                                     [queryPattern => ""] );

  Parse the result of a search and return a SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub parseOutput {
  my %nameValueParams = @_;

  croak $CLASS. "::parseOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  my $DEFILE;
  if ( ref( $nameValueParams{'searchOutput'} ) !~ /GLOB|FileHandle/ ) {
    print $CLASS
        . "::parseOutput() Opening file "
        . $nameValueParams{'searchOutput'} . "\n"
        if ( $DEBUG );
    open $DEFILE, $nameValueParams{'searchOutput'}
        or die $CLASS
        . "::parseOutput: Unable to open "
        . "results file: $nameValueParams{'searchOutput'} : $!";
  }
  else {
    $DEFILE = $nameValueParams{'searchOutput'};
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

  my $resultColl = SearchResultCollection->new();

  while ( <$DEFILE> ) {

    #
    # Conditions for the end of a hit record:
    #   o Must have seen a score
    #   o Must be in the alignment state (ie. have seen Query: and
    #     Subj: recently)
    #   o Must see something which isn't either "Query:" "Subj:" or " "
    #     or must see the end of the file
    #
    if ( $inAlignState ) {
      if ( !/^(Query:|Sbjct:|\s{8}|\n|\r)/ || eof ) {

        $qryLength =~ s/,//g;

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

          # Fix the subject coordinates since they
          # are based on a pre-reversed sequence.
          my $tmpEnd = $sbjEnd;
          $sbjEnd      = $sbjLength - $sbjStart + 1;
          $sbjStart    = $sbjLength - $tmpEnd + 1;
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
          next if ( $qryBase =~ /-|x/i || $sbjBase =~ /-|x/i );
          $baseFreq{$qryBase}++;
          $mismatch++ if ( $qryBase ne $sbjBase );
        }
        my $percDiv =
            sprintf( "%4.2f", $mismatch * 100 / ( $qryEnd + 1 - $qryStart ) );
        my $qgap = $qrySeq =~ tr/-/-/;
        my $sgap = $sbjSeq =~ tr/-/-/;
        my $percIns =
            sprintf( "%4.2f", $sgap * 100 / ( ( $sbjEnd + 1 ) - $sbjStart ) );
        my $percDel =
            sprintf( "%4.2f", $qgap * 100 / ( ( $qryEnd + 1 ) - $qryStart ) );

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
                                     queryString    => $qrySeq,
                                     pctDiverge     => $percDiv,
                                     pctInsert      => $percIns,
                                     pctDelete      => $percDel,
                                     score          => $score
        );

        $resultColl->add( $result );

        $score        = "";
        $sbjSeq       = "";
        $qrySeq       = "";
        $qryOrient    = "";
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
    if ( /Query:\s+(\d+)\s+(\S+)\s+(\d+)/ ) {
      $qrySeq .= $2;
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

    #
    # Subject alignment
    #
    if ( /Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/ ) {
      $sbjSeq .= $2;
      if ( $1 > $3 ) {
        croak "$CLASS::parseOutput: Strange blast alignment format.  "
            . "The subject sequence should not contain descending "
            . "indexes: ( $_ )!\n";
      }
      if ( $sbjStart < 1 ) {
        $sbjStart = _min( $1, $3 );
        $sbjEnd   = _max( $1, $3 );
      }
      else {
        $sbjStart = _min( _min( $1, $3 ), $sbjStart );
        $sbjEnd   = _max( _max( $1, $3 ), $sbjEnd );
      }
      $inAlignState = 1;
    }

    #
    # Query source name
    #
    if ( /^Query=\s+(\S+)/ ) {
      $qryID = $1;
    }

    #
    # Query length
    #
    $qryLength = $1 if ( /^\s+\(([,\d]+) letters/ );

    #
    # Database source name
    #
    #$this->{databaseName} = $1 if ( /Database:\s+(\S+)/ );

    if ( /^\s+Length = (\d+)\s*$/ ) {
      $sbjLength = $1;
    }

    #
    # Score
    #
    if ( /Score.*\((\d+)\)/ ) {    # Get raw score from score line
      $score = $1;
    }

    #
    # Hit description line
    #
    if ( /^>(.*)/ ) {
      $sbjID = $1;
      ## Take note!  This is a very special case here.
      ## In order to search DNA with an expanded alphabet
      ## we have used blastp to search the DNA sequence
      ## as if it were amino acids.  Since blastp does
      ## not search both strands we seed the database with
      ## both the forward and reverse strands labeling the
      ## reverse copies with the string "(anti)".  If this
      ## is located in the results treat this hit as
      ## occuring on the revsere complement of the database
      ## entry.
      if ( $sbjID =~ /anti/ ) {
        $sbjOrient = "C";
      }
      else {
        $sbjOrient = "";
      }
      ( $sbjID ) = $sbjID =~ /(\S+)/;
    }

    #
    # Grab the matrix from the parameter list
    #
    if ( /matrix=\s*(\S+)/ ) {
      my @path = split( /[\\\/]/ );
      $matrix = $path[ $#path ];

      # TODO: Go back and set matrix for all previous hits
    }
  }
  close $DEFILE;

  return $resultColl;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _complexityAdjust( $swnScore, $qrySeq, $sbjSeq, $matLambda,
##                            $matAlphabet, $matScoresRef, $matFreqsRef );
##
##      $swnScore    : Smith-Waterman Raw Alignment Score
##      $qrySeq      : Query string from the alignment
##      $sbjSeq      : Subject string from the alignment
##      $matLambda   : Matrix Lambda parameter
##      $matAlphabet : Matrix alphabet string (in column order)
##      $matScoresRef: Score Matrix
##      $matFreqsRef : Matrix alphabet freq vector
##
## Returns
##      $adj_score : The complexity adjusted score according to
##                   to Phil Green's swat/cross_match program.
##-------------------------------------------------------------------------##
sub _complexityAdjust {
  my ( $swnScore, $qrySeq, $sbjSeq, $matLambda, $matAlphabet, $matScoresRef,
       $matFreqsRef )
      = @_;
  my $mCol      = 0;
  my $t_factor  = 0;
  my $t_sum     = 0;
  my $t_counts  = 0;
  my $n_letters = 0;
  my $baseIndex = 0;
  my $qBase     = "";
  my $sBase     = "";
  my $pos_score = 0;
  my @matCounts = ();
  my $adj_score = 0;

  # Creates a vector of base counts from the query sequence.  It
  # ignores base instances when they are part of a deletion or
  # are insertion characters "-".
  for ( $baseIndex = 0 ; $baseIndex < length( $qrySeq ) ; $baseIndex++ ) {
    $qBase = substr( $qrySeq, $baseIndex, 1 );
    $sBase = substr( $sbjSeq, $baseIndex, 1 );
    if ( $qBase ne "-" && $sBase ne "-" ) {
      $pos_score +=
          $$matScoresRef[ index( $matAlphabet, $qBase ) ]
          [ index( $matAlphabet, $sBase ) ];
      $matCounts[ index( $matAlphabet, $qBase ) ]++
          if ( ( index $matAlphabet, $qBase ) >= 0 );
    }
  }

  #
  #
  #
  for ( $mCol = 0 ; $mCol < length( $matAlphabet ) ; $mCol++ ) {
    if ( defined $matCounts[ $mCol ] && $matCounts[ $mCol ] > 0 ) {
      if ( $$matFreqsRef[ $mCol ] > 0 && log( $$matFreqsRef[ $mCol ] ) != 0 ) {
        $t_factor += $matCounts[ $mCol ] * log( $matCounts[ $mCol ] );
        $t_sum    += $matCounts[ $mCol ] * log( $$matFreqsRef[ $mCol ] );
        $t_counts += $matCounts[ $mCol ];
        $n_letters++;
      }
    }
  }

  #
  #
  #
  $t_factor -= $t_counts * log( $t_counts );
  $t_sum    -= $t_factor;

  #
  # Looks like Phil changed his mind here
  #
  #my $complexity_factor = 0.25;
  #if ( $n_letters > ( 1 / $complexity_factor )) {
  #  $t_factor /= $t_counts * log(1 / $n_letters);
  #}else {
  #  $t_factor /= $t_counts * log($complexity_factor);
  #}
  #
  #$old_adj_score = $swnScore + ($t_factor * $pos_score) - $pos_score + .5;

  #
  #
  #
  $adj_score = sprintf( "%0.0d", $swnScore + $t_sum / $matLambda + .999 );

  $adj_score = 0 if ( !( $adj_score =~ /\d+/ ) || $adj_score < 0 );
  return ( $adj_score );
}

##-------------------------------------------------------------------------##
## Use: my _calculateLambda( $matScoresRef, $matFreqsRef );
##
##      $matScoresRef: Score Matrix
##      $matFreqsRef : Matrix alphabet freq vector
##
## Returns
##      $lambda : The lambda parameter derived from the matrix
##                and the matrix alphabet frequencies.  This
##                is derived from Phil Green's swat/cross_match
##                programs.
##-------------------------------------------------------------------------##
sub _calculateLambda {
  my ( $matScoresRef, $matFreqsRef ) = @_;
  my $lambda_upper = 0;
  my $lambda_lower = 0;
  my $lambda       = 0.5;
  my $sum          = 0;

  do {
    $sum = _getSum( $lambda, $matScoresRef, $matFreqsRef );
    if ( $sum < 1.0 ) {
      $lambda_lower = $lambda;
      $lambda *= 2.0;
    }
  } while ( $sum < 1.0 );

  $lambda_upper = $lambda;

  while ( $lambda_upper - $lambda_lower > .00001 ) {
    $lambda = ( $lambda_lower + $lambda_upper ) / 2.0;
    $sum = _getSum( $lambda, $matScoresRef, $matFreqsRef );
    if ( $sum >= 1.0 ) {
      $lambda_upper = $lambda;
    }
    else {
      $lambda_lower = $lambda;
    }
  }
  return ( $lambda );
}

##-------------------------------------------------------------------------##
## Use: my _getSum( $lambda, $matScoresRef, $matFreqsRef );
##
##      $matScoresRef: Score Matrix
##      $matFreqsRef : Matrix alphabet frequency vector
##
## Returns
##      $sum : Good question....???  This is used by the
##             lambda estimation procedure.  _getSum is
##             derived from Phil Green's swat/cross_match
##             programs.
##-------------------------------------------------------------------------##
sub _getSum {
  my ( $lambda, $matScoresRef, $matFreqsRef ) = @_;
  my $check = 0;
  my $total = 0;
  my $i     = 0;
  my $j     = 0;

  for ( $i = 0 ; $i <= $#$matFreqsRef ; $i++ ) {
    for ( $j = 0 ; $j <= $#$matFreqsRef ; $j++ ) {
      if ( $$matFreqsRef[ $i ] && $$matFreqsRef[ $j ] ) {
        $total += $$matFreqsRef[ $i ] * $$matFreqsRef[ $j ] *
            exp( $lambda * $$matScoresRef[ $i ][ $j ] );
        $check += $$matFreqsRef[ $i ] * $$matFreqsRef[ $j ];
      }
    }
  }
  die "error in _getSum!! check=$check\n"
      if (    $check > 1.001
           || $check < .999 );

  return ( $total );
}

## TODO: Use SequenceSimilarityMatrix.pm here
##-------------------------------------------------------------------------##
## Use: my ( \@matrix, $alphabet, \@freqArray ) =
##                                        _readMatrix( $matrixFileName );
##
##      $matFileName : Name of the file containing the Blast Matrix
##
## Returns
##      This is a very specific matrix reader for Blast matrices
##      which have been created for RepeatMasker.  This routine
##      assumes that the standard crossmatch matrix frequency line
##      is included in the comment section. Ie.:
##
##        # FREQS A 0.325 C 0.175 G 0.175 T 0.325
##
##      The routine returns the matrix stored as a 2-d array.  It also
##      returns the alphabet as a compressed string and the frequencies
##      as an array ( sizeof alphabet ).
##
##-------------------------------------------------------------------------##
sub _readMatrix {
  my ( $matrixFileName ) = @_;
  my %freqHash           = ();
  my @matrix             = ();
  my @rowValues          = ();
  my $row                = 0;
  my $alphabet           = "";
  my $i                  = 0;
  my @freqArray          = ();

  open MATRIX, "<$matrixFileName" || die "Can't open $matrixFileName!\n";

  while ( <MATRIX> ) {
    $_ = uc( $_ );
    chomp;
    if ( /FREQS\s+(.*)/ ) {
      %freqHash = split " ", $1;
    }
    if ( /^\s*[A-Z]\s+[A-Z]\s+[A-Z]\s+[A-Z]\s+/ ) {
      s/ //g;
      $alphabet = $_;
    }
    elsif ( $alphabet
            && /^\s*\S\s+[\d-]+\s+[\d-]+\s+[\d-]+\s+[\d-]+\s+/ )
    {
      @rowValues = split;
      shift @rowValues;
      $matrix[ $row++ ] = [ @rowValues ];
    }
  }
  close MATRIX;

  if ( scalar( keys %freqHash ) == 0 ) {
    $freqHash{"A"} = 0.25;
    $freqHash{"C"} = 0.25;
    $freqHash{"G"} = 0.25;
    $freqHash{"T"} = 0.25;
  }

  for ( $i = 0 ; $i < length( $alphabet ) ; $i++ ) {
    $freqArray[ $i ] = $freqHash{ substr( $alphabet, $i, 1 ) } || 0;
  }

  return ( \@matrix, $alphabet, \@freqArray );
}

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
