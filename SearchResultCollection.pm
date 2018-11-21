#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) SearchResultCollection.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An datastructure for holding search results
##      from various search engines.
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
# bless( { 'results' => [ $SearchResultRef1,
#                         $SearchResultRef2,
#                         ..
#                       ]
#        }
#     }, 'SearchResultCollection' );
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

SearchResultCollection

=head1 SYNOPSIS

use SearchResultCollection

Usage: 

    $SearchResultsCollection = SearchResultCollection->new();

=head1 DESCRIPTION

A class for storing the results from a sequence search engine.
NOTE: This is basically an ArrayList with a specialized write
method.  See ArrayList.pm for accessors methods.

=head1 INSTANCE METHODS

=cut 

package SearchResultCollection;
use strict;
use SearchResult;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require ArrayList;
require Exporter;

@ISA = qw(ArrayList Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "SearchResultCollection";
my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class = shift;

  # TODO: Make this handle filenames!!!!
  #       ie. move parser back!
  my $this = $class->SUPER::new( @_ );

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 write()

  Use: $obj->write( $file, $alignmentMode );

    $file : String
    $alignmentMode : SearchResult::NoAlign
                     SearchResult::AlignWithQuerySeq
                     SearchResult::AlignWithSubjSeq

  Write object to a file.  If $printAlignments
  then include the alignment data.

=cut

##-------------------------------------------------------------------------##
sub write {
  my $this          = shift;
  my $filename      = shift;
  my $alignmentMode = shift;

  open OUT, ">$filename"
      or croak "SearchResultCollection::write - "
      . "Error...could not open file $filename "
      . "for writing!\n";

  for ( my $i = 0 ; $i < $this->size() ; $i++ ) {
    print OUT $this->get( $i )->toStringFormatted( $alignmentMode );
  }
  close OUT;
}

##-------------------------------------------------------------------------##

=head2 toString()

  Use: my $str = $obj->toString( $alignmentMode );

    $alignmentMode : SearchResult::NoAlign
                     SearchResult::AlignWithQuerySeq
                     SearchResult::AlignWithSubjSeq

  Create a monolithic string with all the object data
  formatted for output.

=cut

##-------------------------------------------------------------------------##
sub toString() {
  my $this          = shift;
  my $alignmentMode = shift;

  my $str = "";
  for ( my $i = 0 ; $i < $this->size() ; $i++ ) {
    $str .= $this->get( $i )->toStringFormatted( $alignmentMode );
  }
  return $str;
}

##-------------------------------------------------------------------------##

=head2 maskLevelFilter()

  Use: $numRemoved = $obj->maskLevelFilter(
                                -value => 80,
                                [-inputSeq => SearchResult::Query
                                            SearchResult::Subject] );

    value : Integer - masklevel
    inputSeq: Integer - SearchResult::Query (default) or 
                      SearchResult::Subject.  

  The masklevel controls the reporting of matches based on the 
  overlap of aligned bases.  This routine implements the following
  method: 

      A match is reported if and only if at least (100 - masklevel)% 
      of the aligned bases are not contained within the domain of any
      higher-scoring alignments.

  We achieve this using a greedy algorithm that processes the alignments
  in score ( high to low ) sort order.  Notably this differs from the 
  current way crossmatch implements it's "masklevel" functionality.  
  Crossmatch will remove a lower scoring hit if either of the alignments
  does not have significant non-overlapping extension ( ie. a longer
  lower scoring hit will be removed to favor a contained short higher 
  scoring hit ).  

  Returns the count of entries filtered out. Does *not* alter
  the sort order of the remaining SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub maskLevelFilter {
  my $this       = shift;
  my %parameters = @_;

  croak $CLASS
      . "::maskLevelFilter: Missing or invalid mask level value "
      . "parameter ( value = "
      . $parameters{'value'} . " )!\n"
      if ( !( $parameters{'value'} =~ /\d+/ ) );
  my $maskLevel = $parameters{'value'};

  my $DEBUG = 0;

  croak $CLASS
      . "::maskLevelFilter: Invalid inputSeq parameter "
      . "parameter ( inputSeq = "
      . $parameters{'inputSeq'}
      . ")!\n"
      if (
           defined $parameters{'inputSeq'}
           && (    $parameters{'inputSeq'} != SearchResult::Query
                || $parameters{'inputSeq'} != SearchResult::Subject )
      );

  #
  # Filter set using "masklevel" concept
  #
  return ( 0 ) if ( $maskLevel >= 101 );

  my $inputSeqName  = "getQueryName";
  my $inputSeqStart = "getQueryStart";
  my $inputSeqEnd   = "getQueryEnd";
  if ( defined $parameters{'inputSeq'}
       && $parameters{'inputSeq'} == SearchResult::Subject )
  {
    $inputSeqName  = "getSubjName";
    $inputSeqStart = "getSubjStart";
    $inputSeqEnd   = "getSubjEnd";
  }

  # Sort by inputSeq and then scores high to low
  my @indices = sort {
    $this->get( $a )->$inputSeqName() cmp $this->get( $b )->$inputSeqName()
        || $this->get( $b )->getScore <=> $this->get( $a )->getScore
        || (
       ( $this->get( $a )->getQueryEnd() - $this->get( $a )->getQueryStart() )
       <=> (
         $this->get( $b )->getQueryEnd() - $this->get( $b )->getQueryStart() ) )
  } ( 0 .. $this->size() - 1 );

  my %deleteHash = ();    # Hash to hold indexes of filtered entries
  my %domains    = ();
  for ( my $i = 0 ; $i <= $#indices ; $i++ ) {
    next if ( defined $deleteHash{ $indices[ $i ] } );
    my $result = $this->get( $indices[ $i ] );
    my $name1  = $result->$inputSeqName();
    my $begin1 = $result->$inputSeqStart();
    my $end1   = $result->$inputSeqEnd();

    if ( $DEBUG ) {
      print "Processing: " . $result->getScore() . " $name1:$begin1-$end1\n";
    }

    my $intersected = 0;
    my $j           = 0;
    for ( $j = 0 ; $j <= $#{ $domains{$name1} } ; $j++ ) {
      my $domain = $this->get( $domains{$name1}->[ $j ] );

      # Get members
      my $domainBeg = $domain->$inputSeqStart();
      my $domainEnd = $domain->$inputSeqEnd();
      print "   -- vs $domainBeg-$domainEnd : "
          . ( $#{ $domains{$name1} } - $j )
          . " domains remaining\n"
          if ( $DEBUG );

      # keep domain list sorted by query start
      last if ( $domainBeg > $end1 );

      # Check if they overlap
      next if ( $domainBeg > $end1 || $domainEnd < $begin1 );

      # Calc overlap extension
      my $highestStart = $begin1;
      $highestStart = $domainBeg if ( $domainBeg > $begin1 );
      my $lowestEnd = $end1;
      $lowestEnd = $domainEnd if ( $domainEnd < $end1 );
      my $intersection      = $lowestEnd - $highestStart + 1;
      my $lowScoreExtension = ( $end1 - $begin1 + 1 ) - $intersection;

      # % of the low scoring HSP outside the domain of the
      # high scoring HSP
      my $lowScoreExtPerc = sprintf(
                                     "%0.0f",
                                     (
                                       $lowScoreExtension /
                                           ( $end1 - $begin1 + 1 )
                                         ) * 100
      );

      if ( $DEBUG ) {
        print "    ---> Intersection: "
            . " $intersection  LowScoreExt: $lowScoreExtension LowScoreExtPerc: $lowScoreExtPerc\n";
      }

      if ( $lowScoreExtPerc <= ( 100 - $maskLevel ) ) {
        print "  Deleting hit $begin1-$end1 ("
            . $result->getScore()
            . ") because of domain $domainBeg-$domainEnd ("
            . $domain->getScore() . ")\n"
            if ( $DEBUG );
        $deleteHash{ $indices[ $i ] } = 1;
        $intersected = 1;
      }
    }    # for ( my $j = $i + 1...
    if ( !$intersected ) {
      my $lastDomainIdx = $#{ $domains{$name1} };
      if (    $lastDomainIdx > 0
           && $begin1 >
           $this->get( $domains{$name1}->[ $lastDomainIdx ] )->getQueryStart() )
      {
        print "appending domain\n" if ( $DEBUG );
        push @{ $domains{$name1} }, $indices[ $i ];
      }
      else {
        while ( $j > 0 ) {
          last
              if ( $begin1 >
                  $this->get( $domains{$name1}->[ $j - 1 ] )->getQueryStart() );
          $j--;
        }
        print "Splicing in at $j\n" if ( $DEBUG );
        splice( @{ $domains{$name1} }, $j, 0, $indices[ $i ] );
      }
    }
  }    # for ( my $i = 0...
  undef @indices;

  if ( $DEBUG ) {
    foreach my $name ( %domains ) {
      for ( my $j = 0 ; $j < $#{ $domains{$name} } ; $j++ ) {
        my $domain      = $this->get( $domains{$name}->[ $j ] );
        my $next_domain = $this->get( $domains{$name}->[ $j ] );
        if ( $domain->getQueryStart() > $next_domain->getQueryStart() ) {
          warn "Ooops the domain structure is not in sorted order!\n";
        }
      }
    }
  }

  my $numRemoved = scalar( keys( %deleteHash ) );
  if ( keys( %deleteHash ) ) {
    foreach my $index ( sort { $b <=> $a } keys( %deleteHash ) ) {
      $this->remove( $index );
    }
  }
  undef %deleteHash;

  return ( $numRemoved );

}    # sub maskLevelFilter

##-------------------------------------------------------------------------##

=head2 filterContainedResults()

  Use: $numRemoved = $obj->filterContainedResults(
                                -value => 80,
                                [-inputSeq => SearchResult::Query
                                            SearchResult::Subject] );

    value : Integer - masklevel
    inputSeq: Integer - SearchResult::Query (default) or 
                      SearchResult::Subject.  

  Returns the count of entries filtered out. Does *not* alter
  the sort order of the remaining SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub filterContainedResults {
  my $this       = shift;
  my %parameters = @_;

  my $DEBUG      = 0;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak $CLASS
      . "::$subroutine: Missing or invalid mask level value "
      . "parameter ( value = "
      . $parameters{'value'} . " )!\n"
      if ( !( $parameters{'value'} =~ /\d+/ ) );
  my $maskLevel = $parameters{'value'};

  croak $CLASS
      . "::$subroutine: Invalid inputSeq parameter "
      . "parameter ( inputSeq = "
      . $parameters{'inputSeq'}
      . ")!\n"
      if (
           defined $parameters{'inputSeq'}
           && (    $parameters{'inputSeq'} != SearchResult::Query
                || $parameters{'inputSeq'} != SearchResult::Subject )
      );

  #
  # Filter set using "masklevel" concept
  #
  return ( 0 ) if ( $maskLevel >= 101 );

  my $inputSeqName  = "getQueryName";
  my $inputSeqStart = "getQueryStart";
  my $inputSeqEnd   = "getQueryEnd";
  if ( defined $parameters{'inputSeq'}
       && $parameters{'inputSeq'} == SearchResult::Subject )
  {
    $inputSeqName  = "getSubjName";
    $inputSeqStart = "getSubjStart";
    $inputSeqEnd   = "getSubjEnd";
  }

  # Sort by inputSeq, start position and longest
  my @indices = sort {
    $this->get( $a )->$inputSeqName() cmp $this->get( $b )->$inputSeqName()
        || $this->get( $a )->$inputSeqStart()
        <=> $this->get( $b )->$inputSeqStart()
        || $this->get( $b )->$inputSeqEnd() <=> $this->get( $a )->$inputSeqEnd()
  } ( 0 .. ( $this->size() - 1 ) );

  my %deleteHash = ();    # Hash to hold indexes of filtered entries
  for ( my $i = 0 ; $i <= $#indices ; $i++ ) {
    my $result1 = $this->get( $indices[ $i ] );
    my $name1   = $result1->$inputSeqName();
    my $begin1  = $result1->$inputSeqStart();
    my $end1    = $result1->$inputSeqEnd();
    my $score1  = $result1->getScore();

    my $j = 0;
    for ( $j = $i + 1 ; $j <= $#indices ; $j++ ) {
      next if ( $deleteHash{ $indices[ $j ] } );
      my $result2 = $this->get( $indices[ $j ] );

      # Get members
      my $begin2 = $result2->$inputSeqStart();
      my $end2   = $result2->$inputSeqEnd();
      my $score2 = $result2->getScore();

      # keep domain list sorted by query start
      last if ( $begin2 > $end1 );

      if ( $end2 < $end1 && $score1 > $score2 ) {

       # contained and lower scoring
       # print "Delete: $begin2-$end2 ( " . $result2->getScore() .
       #       " ) because of $begin1-$end1 ( " . $result1->getScore() . " )\n";
        $deleteHash{ $indices[ $j ] } = 1;
      }
    }
  }
  undef @indices;

  my $numRemoved = scalar( keys( %deleteHash ) );
  if ( keys( %deleteHash ) ) {
    foreach my $index ( sort { $b <=> $a } keys( %deleteHash ) ) {
      $this->remove( $index );
    }
  }
  undef %deleteHash;

  return ( $numRemoved );

}    # sub filterContainedResults

##-------------------------------------------------------------------------##

=head2 maskLevelFilterMaskerAid()

  Use: $numRemoved = $obj->maskLevelFilter(
                                -value => 80,
                                [-inputSeq => SearchResult::Query
                                            SearchResult::Subject] );

    value : Integer - masklevel
    inputSeq: Integer - SearchResult::Query (default) or 
                      SearchResult::Subject.  

  The masklevel controls the reporting of matches based on the 
  overlap of aligned bases.  Typically a match is reported only
  if at least (100 - masklevel)% of the bases in its "domain"
  (the part of the query that is aligned) are not contained within
  the domain of any higher-scoring match.

  Returns the count of entries filtered out. Does *not* alter
  the sort order of the remaining SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub maskLevelFilterMaskerAid {
  my $this       = shift;
  my %parameters = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak $CLASS
      . "::$subroutine: Missing or invalid mask level value "
      . "parameter ( value = "
      . $parameters{'value'} . " )!\n"
      if ( !( $parameters{'value'} =~ /\d+/ ) );
  my $maskLevel = $parameters{'value'};

  croak $CLASS
      . "::$subroutine: Invalid inputSeq parameter "
      . "parameter ( inputSeq = "
      . $parameters{'inputSeq'}
      . ")!\n"
      if (
           defined $parameters{'inputSeq'}
           && (    $parameters{'inputSeq'} != SearchResult::Query
                || $parameters{'inputSeq'} != SearchResult::Subject )
      );

  #
  # Filter set using "masklevel" concept
  #
  return ( 0 ) if ( $maskLevel >= 101 );

  my $inputSeqName  = "getQueryName";
  my $inputSeqStart = "getQueryStart";
  my $inputSeqEnd   = "getQueryEnd";
  if ( defined $parameters{'inputSeq'}
       && $parameters{'inputSeq'} == SearchResult::Subject )
  {
    $inputSeqName  = "getSubjName";
    $inputSeqStart = "getSubjStart";
    $inputSeqEnd   = "getSubjEnd";
  }

  # Sort by inputSeq and then scores high to low
  my @indices = sort {
    $this->get( $a )->$inputSeqName() cmp $this->get( $b )->$inputSeqName()
        || $this->get( $b )->getScore <=> $this->get( $a )->getScore
  } ( 0 .. $this->size() - 1 );

  my %deleteHash = ();    # Hash to hold indexes of filtered entries
  for ( my $i = 0 ; $i <= $#indices ; $i++ ) {
    next if ( defined $deleteHash{ $indices[ $i ] } );

    for ( my $j = $i + 1 ; $j <= $#indices ; $j++ ) {
      next if ( defined $deleteHash{ $indices[ $j ] } );

      my $result1 = $this->get( $indices[ $i ] );
      my $result2 = $this->get( $indices[ $j ] );

      last if ( $result1->$inputSeqName() ne $result2->$inputSeqName() );

      # Get members
      my $begin1 = $result1->$inputSeqStart();
      my $begin2 = $result2->$inputSeqStart();
      my $end1   = $result1->$inputSeqEnd();
      my $end2   = $result2->$inputSeqEnd();

      # Check if they overlap
      next if ( $begin2 > $end1 || $end2 < $begin1 );
      next if ( $begin2 > $end1 || $end2 < $begin1 );

      # Calc overlap extension
      my $extension = 0;
      $extension = $begin1 - $begin2 if ( $begin2 < $begin1 );
      $extension += $end2 - $end1 if ( $end2 > $end1 );

      # % of the low scoring HSP outside the domain of the
      # high scoring HSP
      my $perc = ( $extension / ( $end2 - $begin2 + 1 ) ) * 100;
      if ( $perc < ( 100 - $maskLevel ) ) {

        #if ( $perc == 0 )
        #{
        #  $deleteHash{ $indices[ $j ] } = 2;
        #} else
        #{
        $deleteHash{ $indices[ $j ] } = 1;

        #}
      }
    }    # for ( my $j = $i + 1...
  }    # for ( my $i = 0...
  undef @indices;

  # Remove all hits which were filtered above
  my $numRemoved = scalar( keys( %deleteHash ) );
  if ( keys( %deleteHash ) ) {
    foreach my $index ( sort { $b <=> $a } keys( %deleteHash ) ) {
      $this->remove( $index );
    }
  }
  undef %deleteHash;

  return ( $numRemoved );

}    # sub maskLevelFilterMaskerAid

##-------------------------------------------------------------------------##

=head2 crossmatchMaskLevelFilter()

  Use: $numRemoved = $obj->crossmatchMaskLevelFilter(
                                -value => 80,
                                [-inputSeq => SearchResult::Query
                                            SearchResult::Subject] );

    value : Integer - masklevel
    inputSeq: Integer - SearchResult::Query (default) or 
                      SearchResult::Subject.  

  Crossmatch has the ability to filter repeats based on a
  overlap filter called the masklevel.  From the phrap docs:

     "Masklevel controls the grouping of matches according to 
      their domains (the segment of the query that is aligned), which 
      is useful when different portions of the query may match
      different subject sequence regions (e.g. cDNA queries vs 
      genomic subject, or chimeric sequencing read queries).  
      Two matches for the same query are considered to be in the 
      same "query domain group" if at least masklevel% of the bases 
      in the domain of either one of them is contained within the
      domain of the other. Thus for the default value of 80, matches
      are assigned to the same group if the domain of either one
      is at least 80% contained within the domain of the other."

  Thus theoretically this would remove the longer lower scoring 
  hit in this set at masklevel 99:

  460  6.67 0.00 0.00  Seq1 1 60 (40) Cons2 1  60 (0)  <==DELETE
  484  0.00 0.00 0.00  Seq1 1 50 (50) Cons1 11 60 (0)  

  i.e:
  -------------------DELETE-----------------------------------> 
  -------------------------------------------------->

  despite the fact that the low scoring hit has 16% of its bases 
  outside the query domain.  This is because the other hit is fully
  contained ie it has 0% of it's bases outside the query domain.
  This is not always always the desired behaviour for an result
  filter, in RepeatMasker ( and the original maskeraid ) we
  prefer longer matches or at least to not to remove them from the result
  set.  In RM's implementation of masklevel both the results above will
  be printed.  In most cases the results will be the same as the 
  longer alignments tend to have higher scores. 

  Until recently we didn't appreciate quite fully how this was 
  implemented.  This routine attempts to match closely the
  actual crossmatch implementation.  

  Crossmatch uses a greedy algorithm, as most tiling path methods
  do, to implement the masklevel filter.  As with most greedy
  algorithms the order of operations is important and can 
  significantly change the results.  Instead of using a high to
  low sort order crossmatch filters in order of received hits.  
  This works out to be a sort order like:

        1. Query Sequence
        2. ID order in subject Library
        3. Orientation 
                - Reverse inputSeq first, query start low to high
                - Forward inputSeq second, query high to low

  This routine requires that the search result collection be
  pre-ordered as stated above.  Here is an example of how to 
  do this given a hash of subject IDs containing the numbered
  position in the subject file:

   $searchResults->sort(
      sub ($$) {
        $_[ 0 ]->getQueryName cmp $_[ 1 ]->getQueryName()
      || $idOrder{$_[ 0 ]->getSubjName()} <=> $idOrder{$_[ 1 ]->getSubjName()}
      || $_[ 1 ]->getOrientation() cmp $_[ 0 ]->getOrientation()
      || ( $_[ 1 ]->getOrientation() eq "C" && 
           ( $_[ 0 ]->getQueryStart() <=> $_[ 1 ]->getQueryStart() ) )
      || ( $_[ 1 ]->getOrientation() ne "C" && 
           (  $_[ 1 ]->getQueryStart() <=> $_[ 0 ]->getQueryStart() ) );
      }
    );

  NOTE: This is an experimental routine and should only be used for
  testing at this time. 

=cut

##-------------------------------------------------------------------------##
sub crossmatchMaskLevelFilter {
  my $this       = shift;
  my %parameters = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak $CLASS
      . "::$subroutine: Missing or invalid mask level value "
      . "parameter ( value = "
      . $parameters{'value'} . " )!\n"
      if ( !( $parameters{'value'} =~ /\d+/ ) );
  my $maskLevel = $parameters{'value'};

  croak $CLASS
      . "::$subroutine: Invalid inputSeq parameter "
      . "parameter ( inputSeq = "
      . $parameters{'inputSeq'}
      . ")!\n"
      if (
           defined $parameters{'inputSeq'}
           && (    $parameters{'inputSeq'} != SearchResult::Query
                || $parameters{'inputSeq'} != SearchResult::Subject )
      );

  #
  # Filter set using "masklevel" concept
  #
  return ( 0 ) if ( $maskLevel >= 101 );

  my $inputSeqName  = "getQueryName";
  my $inputSeqStart = "getQueryStart";
  my $inputSeqEnd   = "getQueryEnd";
  if ( defined $parameters{'inputSeq'}
       && $parameters{'inputSeq'} == SearchResult::Subject )
  {
    $inputSeqName  = "getSubjName";
    $inputSeqStart = "getSubjStart";
    $inputSeqEnd   = "getSubjEnd";

  }

  my @indices = ( 0 .. $this->size() - 1 );

  my @queryDomains = ( $this->get( $indices[ 0 ] ) );

  for ( my $i = 1 ; $i <= $#indices ; $i++ ) {
    my $intersected = 0;
    my $result1     = $this->get( $indices[ $i ] );
    my $begin1      = $result1->$inputSeqStart();
    my $end1        = $result1->$inputSeqEnd();
    for ( my $j = 0 ; $j <= $#queryDomains ; $j++ ) {
      my $result2 = $queryDomains[ $j ];

      # Get members
      my $begin2 = $result2->$inputSeqStart();
      my $end2   = $result2->$inputSeqEnd();

      # Check if they overlap
      next if ( $begin2 > $end1 || $end2 < $begin1 );
      next if ( $begin2 > $end1 || $end2 < $begin1 );

      # st set to larger begin
      my $st = $begin2;
      $st = $begin1 if ( $begin1 > $begin2 );
      my $ed = $end1;
      $ed = $end2 if ( $end1 > $end2 );
      my $intersect = ( $ed - $st + 1 );
      $intersect = $intersect * ( 100.0 / $maskLevel );

      if (    $intersect >= ( $end1 - $begin1 + 1 )
           || $intersect >= ( $end2 - $begin2 + 1 ) )
      {
        if ( $result1->getScore() > $result2->getScore() ) {

          # Hit is better
          $queryDomains[ $j ] = $result1;

       #print "Mrg: $begin1-$end1 ( " . $result1->getScore() .
       #      " ) subsumed by $begin2-$end2 ( " . $result2->getScore() . " )\n";
          $intersected = 1;
        }
        else {

          # Query domain is better
          $intersected = 1;

    #print "Mrg: $begin2-$end2 ( " . $result2->getScore() .
    #      " ) eats next best $begin1-$end1 ( " . $result1->getScore() . " )\n";
        }
        last;
      }
    }
    if ( !$intersected ) {

      #print "Mrg: $begin1-$end1 ( " . $result1->getScore() . " )\n";
      push @queryDomains, $result1;
    }
  }

  @queryDomains =
      sort { $a->getQueryStart() <=> $b->getQueryStart() } @queryDomains;

  my %deleteHash = ();    # Hash to hold indexes of filtered entries
  for ( my $i = 0 ; $i <= $#queryDomains ; $i++ ) {
    next if ( defined $deleteHash{$i} );
    for ( my $j = $i + 1 ; $j <= $#queryDomains ; $j++ ) {
      next if ( defined $deleteHash{ $indices[ $j ] } );

      my $result1 = $queryDomains[ $i ];
      my $result2 = $queryDomains[ $j ];

      last if ( $result1->$inputSeqName() ne $result2->$inputSeqName() );

      # Get members
      my $begin1 = $result1->$inputSeqStart();
      my $begin2 = $result2->$inputSeqStart();
      my $end1   = $result1->$inputSeqEnd();
      my $end2   = $result2->$inputSeqEnd();

      # Check if they overlap
      next if ( $begin2 > $end1 || $end2 < $begin1 );
      next if ( $begin2 > $end1 || $end2 < $begin1 );

      # st set to larger begin
      my $st = $begin2;
      $st = $begin1 if ( $begin1 > $begin2 );
      my $ed = $end1;
      $ed = $end2 if ( $end1 > $end2 );
      my $intersect = ( $ed - $st + 1 );
      $intersect = $intersect * ( 100.0 / $maskLevel );

      if (    $intersect >= ( $end1 - $begin1 + 1 )
           || $intersect >= ( $end2 - $begin2 + 1 ) )
      {
        if ( $result1->getScore() > $result2->getScore() ) {
          $deleteHash{$j} = 1;
        }
        else {
          $deleteHash{$i} = 1;
        }
        last;
      }
    }    # for ( my $j = $i + 1...
  }    # for ( my $i = 0...
  undef @indices;

  my $origSize = $this->size();
  $this->clear();
  for ( my $i = 0 ; $i <= $#queryDomains ; $i++ ) {
    next if ( $deleteHash{$i} );
    $this->add( $queryDomains[ $i ] );
  }

  return ( $origSize - $this->size() );

}

1;
