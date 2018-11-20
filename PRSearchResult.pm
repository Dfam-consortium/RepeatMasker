#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) PRSearchResult.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An object for holding ProcessRepeats meta data related to a
##      search result.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2011 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
# ChangeLog
#
#     $Log$
#
###############################################################################
# To Do:
#

=head1 NAME

PRSearchResult

=head1 SYNOPSIS

use PRSearchResult

Usage: 

    my $PRSearchResult = PRSearchResult->new();

  or 

    my $PRSearchResultsCollection = PRSearchResult->new( 
                                 queryName=>value, subjName=>value,
                                 pctInsert=>value, pctDelete=>value, 
                                 queryStart=>value, queryEnd=>value, 
                                 score=>value, pctDiverge=>value, 
                                 subjRemaining=>value, subjType=>value,
                                 queryRemaining=>value, id=>value,
                                 orientation=>value, queryString=>value,
                                 subjString=>value, matrix=>value,
                                 id=>value, lineageId=>value );

=head1 DESCRIPTION

A subclass of SearchResult which holds a ProcessRepeats specific search
result.

=head1 SEE ALSO

=back

=head1 COPYRIGHT

Copyright 2011 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package PRSearchResult;
use strict;
use SearchResult;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(SearchResult Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "PRSearchResult";
my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Constructor:
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  my $this = $class->SUPER::new( @_ );

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

# Allow import of values
#if ( %nameValuePairs )
#{
#  while ( my ( $name, $value ) = each( %nameValuePairs ) )
#  {
#    my $method = "set" . _ucFirst( $name );
#    unless ( $this->can( $method ) )
#    {
#      croak(
#           "PRSearchResult::add: Instance variable $name doesn't exist." . "" );
#    }
#    $this->$method( $value );
#  }
#}

  return $this;
}

##-------------------------------------------------------------------------##

=head2 clone()

  Use: my $newObj = $obj->clone();

  Clone a PRSearchResult *duplicating* all the values of the old
  object in the new one.

=cut

##-------------------------------------------------------------------------##
sub clone {
  my $this = shift;

  my %newHash = %{$this};
  my $newObj  = \%newHash;

  bless $newObj, ref( $this );

  return $newObj;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

# Just these fields for now
sub setFrom {
  my $this   = shift;
  my $source = shift;

  $this->setScore( $source->getScore() );
  $this->setPctDiverge( $source->getPctDiverge() );
  $this->setPctKimuraDiverge( $source->getPctKimuraDiverge() );
  $this->setPctDelete( $source->getPctDelete() );
  $this->setPctInsert( $source->getPctInsert() );
  $this->setQueryName( $source->getQueryName() );
  $this->setQueryStart( $source->getQueryStart() );
  $this->setQueryEnd( $source->getQueryEnd() );
  $this->setQueryRemaining( $source->getQueryRemaining() );
  $this->setOrientation( $source->getOrientation() );
  $this->setHitName( $source->getHitName() );
  $this->setClassName( $source->getClassName() );
  $this->setSubjStart( $source->getSubjStart() );
  $this->setSubjEnd( $source->getSubjEnd() );
  $this->setSubjRemaining( $source->getSubjRemaining() );
  $this->setLineageId( $source->getLineageId() );
  $this->setId( $source->getId() );
  $this->setPRID( $source->getPRID() );
  $this->setEquivHash( $source->getEquivHash() );
  $this->setSubjString( $source->getSubjString() );
  $this->setQueryString( $source->getQueryString() );
  $this->setMatrixName( $source->getMatrixName() );

  # WARNING: Do not inherit "derivedFromAnnots".
  # These are the fields we are ignoring currently
  #    setStatus
  #    setLeftLinkedHit
  #    setRightLinkedHit
  #    ..

  return $this;

}

## DEPRECATED
# The only use of this by PR was to clear the data to save memory
sub setAlignData {
  my $this = shift;
  $this->setSubjString( "" );
  $this->setQueryString( "" );
}

## DEPRECATED
sub getAlignData {
  my $this = shift;
  my $str  = $this->toStringFormatted( SearchResult::AlignWithQuerySeq );
  $str =~ s/^\s*\d+[^\n\r]+$//;
  return $str;
}

sub getSeq2Len {
  my $this = shift;
  return ( $this->getSubjEnd() + $this->getSubjRemaining() );
}

## beware the orientation is "+" "C"

sub getEquivHash {
  my $this = shift;
  return $this->{'equivHash'};
}

sub setEquivHash {
  my $this  = shift;
  my $value = shift;
  $this->{'equivHash'} = $value;
}

# *.out: id lineageId overlap
# *.cat: stage refid  ( getLastField, getCatID )

# get/set prID
# get/set lineageId

# The final ID tag in Process Repeats
sub setPRID {
  my $this  = shift;
  my $value = shift;
  $this->{'prID'} = $value;
}

sub getPRID {
  my $this = shift;
  return $this->{'prID'};
}

sub isCut {
  my $this = shift;
  if ( $this->getLineageId() =~ /c_b\d+s\d+i\d+/ ) {
    return 1;
  }
  return 0;
}

sub isMasked {
  my $this = shift;
  if ( $this->getLineageId() =~ /m_b\d+s\d+i\d+/ ) {
    return 1;
  }
  return 0;
}

sub getStage {
  my $this = shift;
  if ( $this->getLineageId() =~ /[mc]_b\d+s(\d+)i\d+/ ) {
    return ( $1 );
  }
  return -1;
}

## This is used in PR to hold the state of an annotation.
## This is most often used to flag records to be deleted.
sub getStatus {
  my $this = shift;
  return $this->{'status'};
}

sub setStatus {
  my $this  = shift;
  my $value = shift;
  $this->{'status'} = $value;
}

sub getBlockLeftJoin {
  my $this = shift;
  return $this->{'blockLeftJoin'};
}

sub setBlockLeftJoin {
  my $this  = shift;
  my $value = shift;
  $this->{'blockLeftJoin'} = $value;
}

sub getBlockRightJoin {
  my $this = shift;
  return $this->{'blockRightJoin'};
}

sub setBlockRightJoin {
  my $this  = shift;
  my $value = shift;
  $this->{'blockRightJoin'} = $value;
}

sub getLeftLinkedHit {
  my $this = shift;
  return $this->{'leftHit'};
}

sub setLeftLinkedHit {
  my $this  = shift;
  my $value = shift;
  $this->{'leftHit'} = $value;
}

#returns element ( if set ) who is linked to our right edge ( Seq1End )
#
sub getRightLinkedHit {
  my $this = shift;
  return $this->{'rightHit'};
}

sub setRightLinkedHit {
  my $this  = shift;
  my $value = shift;
  $this->{'rightHit'} = $value;
}

sub getClassName {
  my $this = shift;
  my ( $hitName, $className ) =
      split( /\#/, $this->getSubjName() );
  return $className;
}

# TODO Add some error checking
sub setClassName {
  my $this         = shift;
  my $newClassName = shift;
  my ( $hitName, $className ) =
      split( /\#/, $this->getSubjName() );
  $this->setSubjName( $hitName . "#" . $newClassName );
}

sub getHitName {
  my $this = shift;
  my ( $hitName, $className ) =
      split( /\#/, $this->getSubjName() );
  return $hitName;
}

sub setHitName {
  my $this       = shift;
  my $newHitName = shift;
  my ( $hitName, $className ) =
      split( /\#/, $this->getSubjName() );
  $this->setSubjName( $newHitName . "#" . $className );
}

sub getUniqID {
  my $this = shift;
  return (   $this->getScore()
           . $this->getPctDiverge()
           . $this->getPctDelete()
           . $this->getPctInsert()
           . $this->getQueryName()
           . $this->getQueryStart()
           . $this->getQueryEnd()
           . $this->getQueryRemaining()
           . $this->getOrientation()
           . $this->getSubjName()
           . $this->getSubjStart()
           . $this->getSubjEnd()
           . $this->getSubjRemaining() );
}

sub getDerivedFromAnnot {
  my $this = shift;
  return ( $this->{'derivedFrom'} );
}

sub setDerivedFromAnnot {
  my $this   = shift;
  my $member = shift;
  my $annot  = new PRSearchResult;
  $annot->setFrom( $member );
  $this->{'derivedFrom'} = [ $member ];
}

sub addDerivedFromAnnot {
  my $this   = shift;
  my $member = shift;

  my $annot = new PRSearchResult;
  $annot->setFrom( $member );
  ## This needs to be a deep copy
  if ( $member->getDerivedFromAnnot() ) {
    foreach my $subAnnot ( @{ $member->getDerivedFromAnnot() } ) {
      $annot->addDerivedFromAnnot( $subAnnot );
    }
  }
  push @{ $this->{'derivedFrom'} }, $annot;
}

sub print {
  my $this       = shift;
  my $noLineTerm = shift;
  my $lineTerm   = "\n" if ( !defined $noLineTerm );
  my $outStr = sprintf(
                        "%4d %3.1f %3.1f %3.1f %-18.18s %8d %8d "
                            . "%8d %1s %-18.18s %8d %8d %8d %8s %6s$lineTerm",
                        $this->getScore(),         $this->getPctDiverge(),
                        $this->getPctDelete(),     $this->getPctInsert(),
                        $this->getQueryName(),     $this->getQueryStart(),
                        $this->getQueryEnd(),      $this->getQueryRemaining(),
                        $this->getOrientation(),   $this->getSubjName(),
                        $this->getSubjStart(),     $this->getSubjEnd(),
                        $this->getSubjRemaining(), $this->getId(),
                        $this->getStatus()
  );
  print $outStr;
}

sub printBrief {
  my $this = shift;
  my $outStr = sprintf(
    "%4d %-18.18s %8d %8d " . "%1s %-18.18s %8d %8d %20s\n",
    $this->getScore(),
    $this->getQueryName(),   $this->getQueryStart(),
    $this->getQueryEnd(),
    $this->getOrientation(), $this->getSubjName(),
    $this->getSubjStart(), $this->getSubjEnd(), $this->getLineageId()

  );
  print $outStr;
}

sub printLeftRightLinks {
  my $this = shift;
  if ( $this->getLeftLinkedHit() != undef ) {
    if ( $this->getLeftLinkedHit() == $this ) {
      print "    BLOCKED ( links to self )\n";
    }
    else {
      print "    ";
      $this->getLeftLinkedHit()->printBrief();
    }
  }
  else {
    print "    UNDEF\n";
  }
  print " -->";
  $this->printBrief();
  if ( $this->getRightLinkedHit() != undef ) {
    if ( $this->getRightLinkedHit() == $this ) {
      print "    BLOCKED ( links to self )\n";
    }
    else {
      print "    ";
      $this->getRightLinkedHit()->printBrief();
    }
  }
  else {
    print "    UNDEF\n";
  }
}

sub sanityCheckConsPos {
  my $this    = shift;
  my $hashRef = shift;

  # Find begining
  print "Sanity Checking:\n";
  my $firstInChain = $this;
  my $detectLoop   = 0;
  while (    $firstInChain->getLeftLinkedHit() != undef
          && $detectLoop < 50 )
  {

    # Ok to link to self.
    if ( $firstInChain == $firstInChain->getLeftLinkedHit() ) {
      last;
    }
    $firstInChain = $firstInChain->getLeftLinkedHit();
    $detectLoop++;
  }

  if ( $detectLoop >= 50 ) {
    print "WARNING! This chain contains a loop!!!!\n";
  }

  # Now print

  my $nextInChain = $firstInChain;
  do {
    if ( $nextInChain == $this ) {
      print "--> (";
    }
    else {
      print "    (";
    }
    print "" . $hashRef->{ $nextInChain->getId() } . "): ";
    $nextInChain->printBrief();

    # Link to itself is also allowed
    if ( $nextInChain == $nextInChain->getRightLinkedHit() ) {
      last;
    }
    $nextInChain = $nextInChain->getRightLinkedHit();
      } while ( $nextInChain )

}

#
# Print an element and ( if exists ) all of it chain members in order
#
sub printLinks {
  my $this = shift;

  # Find begining
  my $firstInChain = $this;
  my $detectLoop   = 0;
  while (    $firstInChain->getLeftLinkedHit() != undef
          && $firstInChain != $firstInChain->getLeftLinkedHit()
          && $detectLoop < 50 )
  {
    $firstInChain = $firstInChain->getLeftLinkedHit();
    $detectLoop++;
  }

  if ( $detectLoop >= 50 ) {
    print STDERR "WARNING! This chain contains a loop!!!!\n";
  }

  # Now print
  my $nextInChain = $firstInChain;
  do {
    if ( $nextInChain == $this ) {
      print "--> ";
    }
    else {
      print "    ";
    }
    $nextInChain->printBrief();

    # Link to itself is also allowed
    if ( $nextInChain == $nextInChain->getRightLinkedHit() ) {
      last;
    }
    $nextInChain = $nextInChain->getRightLinkedHit();
      } while ( $nextInChain )

}

# Intelligently relink an element's edges prior to removal
sub removeFromJoins {
  my $this = shift;

  my $DEBUG = 0;
  my $left  = $this->getLeftLinkedHit();
  my $right = $this->getRightLinkedHit();

  if ( $DEBUG ) {
    print "Remove from joins:\n";
    $this->printLinks();
  }

  # Just reset these since we are being removed
  $this->setLeftLinkedHit( undef )
      unless ( $this == $left );    # don't removedlinks to itself.
  $this->setRightLinkedHit( undef )
      unless ( $this == $right );    # don't remove links to itself.

  $left  = undef if ( $this == $left );
  $right = undef if ( $this == $right );

  if ( defined $left && !defined $right ) {

    # Remove left's edge reference
    $left->setRightLinkedHit( undef );
  }
  elsif ( defined $left && defined $right ) {

    # Join our neighbors together
    $left->setRightLinkedHit( $right );
    $right->setLeftLinkedHit( $left );
  }
  elsif ( !defined $left && defined $right ) {

    # Remove rights's edge reference
    $right->setLeftLinkedHit( undef );
  }
}

# Intelligently relink an elements edges for insertion
sub join {
  my $this    = shift;
  my $partner = shift;

  if ( $DEBUG ) {
    print "join this:\n";
    $this->printLinks();
    print "to partner\n";
    $partner->printLinks();
  }

  die "join(): Invalid join!  \$this == \$partner" if ( $this == $partner );

  my @cluster = ();
  foreach my $chain ( $this, $partner ) {
    push @cluster, $chain;
    my $nextInChain = $chain;
    my %seen        = ();
    while (    $nextInChain->getLeftLinkedHit()
            && $nextInChain->getLeftLinkedHit() != $nextInChain )
    {
      if ( $seen{ $nextInChain->getLeftLinkedHit() } ) {
        warn "WARN: breaking loop in fragments\n";
        $nextInChain->setLeftLinkedHit( undef );
        last;
      }
      $nextInChain = $nextInChain->getLeftLinkedHit();
      $seen{$nextInChain}++;
      push @cluster, $nextInChain;
    }
    $nextInChain = $chain;
    while (    $nextInChain->getRightLinkedHit()
            && $nextInChain->getRightLinkedHit() != $nextInChain )
    {
      if ( $seen{ $nextInChain->getRightLinkedHit() } ) {
        warn "WARN: breaking loop in fragments\n";
        $nextInChain->setRightLinkedHit( undef );
        last;
      }
      $nextInChain = $nextInChain->getRightLinkedHit();
      $seen{$nextInChain}++;
      push @cluster, $nextInChain;
    }
  }

  # Sort cluster
  @cluster = sort { $a->comparePositionOrder( $b ) } ( @cluster );

  my $lastAnnot = undef;
  foreach my $annot ( @cluster ) {
    next if ( $annot == $lastAnnot );

    $annot->setLeftLinkedHit( $lastAnnot );
    $annot->setRightLinkedHit( undef );
    if ( $lastAnnot ) {
      $lastAnnot->setRightLinkedHit( $annot );
    }
    $lastAnnot = $annot;
  }

  if ( $DEBUG ) {
    print "Now look what we did with this thing:\n";
    $this->printLinks();
  }

}

# Intelligently relink an elements edges for insertion
sub resortJoins {
  my $this = shift;

  my @cluster = ();
  push @cluster, $this;
  my $nextInChain = $this;
  while (    $nextInChain->getLeftLinkedHit()
          && $nextInChain->getLeftLinkedHit() != $nextInChain )
  {
    $nextInChain = $nextInChain->getLeftLinkedHit();
    push @cluster, $nextInChain;
  }
  $nextInChain = $this;
  while (    $nextInChain->getRightLinkedHit()
          && $nextInChain->getRightLinkedHit() != $nextInChain )
  {
    $nextInChain = $nextInChain->getRightLinkedHit();
    push @cluster, $nextInChain;
  }

  # Sort cluster
  @cluster = sort { $a->comparePositionOrder( $b ) } ( @cluster );

  my $lastAnnot = undef;
  foreach my $annot ( @cluster ) {
    next if ( $annot == $lastAnnot );
    $annot->setLeftLinkedHit( $lastAnnot );
    $annot->setRightLinkedHit( undef );
    if ( $lastAnnot ) {
      $lastAnnot->setRightLinkedHit( $annot );
    }
    $lastAnnot = $annot;
  }

  if ( $DEBUG ) {
    print "Now look what we did with this thing:\n";
    $this->printLinks();
  }

}

# Intelligently merge $partner information with us.  Then remove
# our partner's links to ready it for removal.
sub mergeSimpleLow {
  my $this    = shift;
  my $partner = shift;

  my $DEBUG = 0;

  if ( $DEBUG ) {
    print "mergeSimpleLow(): THIS links:\n";
    $this->printLinks();
    print "mergeSimpleLow(): PARTNER links\n";
    $partner->printLinks();
  }

  # Pick appropriate substitution/insertion/deletion stats for
  # merged element.
  if ( $this->getHitName() eq $partner->getHitName() ) {
    $this->adjustSubstLevel( $partner, "noQueryOverlapAdj" );
    $this->setScore( $partner->getScore() )
        if ( $this->getScore() < $partner->getScore() );
  }
  else {
    if ( $this->getScore() < $partner->getScore() ) {
      $this->setPctDiverge( $partner->getPctDiverge() );
      $this->setPctKimuraDiverge( $partner->getPctKimuraDiverge() );
      $this->setPctDelete( $partner->getPctDelete() );
      $this->setPctInsert( $partner->getPctInsert() );
      $this->setScore( $partner->getScore() );
      $this->setSubjName( $partner->getSubjName() );
    }    # else....keep ours
  }

  # Pick the lower query begin position
  $this->setQueryStart( $partner->getQueryStart() )
      if ( $this->getQueryStart() > $partner->getQueryStart() );

  # Pick the higher query end position
  if ( $this->getQueryEnd() < $partner->getQueryEnd() ) {
    $this->setQueryEnd( $partner->getQueryEnd() );
    $this->setQueryRemaining( $partner->getQueryRemaining() );
  }

  my $newQuerySize = $this->getQueryEnd() - $this->getQueryStart() + 1;
  my $newSeq2End   =
      $newQuerySize + ( $this->getPctDelete() - $this->getPctInsert() ) *
      $newQuerySize / 100;
  $this->setSubjStart( 1 );
  $this->setSubjEnd( sprintf( "%d", $newSeq2End ) );

  ## TODO: What happens when our partner has joins that
  ##       should point to the new merger?
  $partner->removeFromJoins();

  if ( $DEBUG ) {
    print "mergeSimpleLow(): Final THIS links:\n";
    $this->printLinks();
    print "mergeSimpleLow(): Leaving\n";
  }

}

# Intelligently merge $partner information with us.  Then remove
# our partner's links to ready it for removal.
sub merge {
  my $this    = shift;
  my $partner = shift;

  my $DEBUG = 0;

  if ( $DEBUG ) {
    print "merge(): Entered\n";
    print "merge(): THIS links:\n";
    $this->printLinks();
    print "merge(): PARTNER links\n";
    $partner->printLinks();
  }

  # Pick appropriate substitution/insertion/deletion stats for
  # merged element.
  $this->adjustSubstLevel( $partner );

  # Pick the lower query begin position
  $this->setQueryStart( $partner->getQueryStart() )
      if ( $this->getQueryStart() > $partner->getQueryStart() );

  # Pick the higher query end position
  if ( $this->getQueryEnd() < $partner->getQueryEnd() ) {
    $this->setQueryEnd( $partner->getQueryEnd() );
    $this->setQueryRemaining( $partner->getQueryRemaining() );
  }

  # Find the lower consensus begin position
  my $minBeg = $partner->getSubjStart();
  $minBeg = $this->getSubjStart()
      if ( $this->getSubjStart() < $minBeg );

  # Find the higher consensus end position
  my $maxEnd = $partner->getSubjEnd();
  $maxEnd = $this->getSubjEnd()
      if ( $this->getSubjEnd() > $maxEnd );

  # Sanity check consensus position choices. If
  # the hybrid range ( lowest begin & highest end ) of
  # the two elements is > 2* query length then
  # use the original range of the higher scoring
  # element instead of the hybrid.
  if ( ( $maxEnd - $minBeg + 1 ) >
       ( 2 * ( $this->getQueryEnd() - $this->getQueryStart() + 1 ) ) )
  {

    # Use the highest scoring annotation's range
    if ( $this->getScore() < $partner->getScore() ) {
      $this->setSubjStart( $partner->getSubjStart() );
      $this->setSubjEnd( $partner->getSubjEnd() );
      $this->setSubjRemaining( $partner->getSubjRemaining() );
    }
  }
  else {

    # Use a hybrid range
    $this->setSubjStart( $minBeg );
    if ( $this->getSubjEnd() < $partner->getSubjEnd() ) {
      $this->setSubjEnd( $maxEnd );
      $this->setSubjRemaining( $partner->getSubjRemaining() );
    }
  }

  # Use the higher score and it's associated consensus name
  if ( $this->getScore() < $partner->getScore() ) {
    $this->setScore( $partner->getScore() );
    $this->setSubjName( $partner->getSubjName() );
  }

  $partner->removeFromJoins();

  if ( $DEBUG ) {
    print "merge(): Final THIS links:\n";
    $this->printLinks();
    print "merge(): Leaving\n";
  }

}

#### NEW TESTING
sub getAdjustedSubstLevel {
  my $this    = shift;
  my $partner = shift;
  my $method  = shift;

  my $DEBUG = 0;
  if ( $DEBUG ) {
    print "getAdjustedSubstLevel( this, partner, method = $method ): Entered\n";
    print "this: ";
    $this->print();
    print "partner: ";
    $partner->print();
  }

  warn "getAdjustedSubstLevel(): Unknown method $method"
      if ( defined $method && $method !~ /noQueryOverlapAdj/ );

  # For calculation, choose best matching consensus in overlapped region
  my $subBases;
  my $delBases;
  my $insBases;

  #
  # This change allows us to calculate these stats based on
  # the complete alignemnt returned by RM.  Previously this would
  # only account for one fragment.
  # TODO: Allow this change after a sanity check!
  #
  my $thisLength = $this->getQueryEnd() - $this->getQueryStart() + 1;

  my $partnerLength = $partner->getQueryEnd() - $partner->getQueryStart() + 1;

  my $qo = 0;

  $qo = $this->getQueryOverlap( $partner )
      unless ( $method eq "noQueryOverlapAdj" );
  my $co = $this->getConsensusOverlap( $partner );
  print "adjustSubstLevel(): qo = $qo, co = $co, "
      . "thisLength = $thisLength, "
      . "partnerLength = $partnerLength\n"
      if ( $DEBUG );

  if ( $this->getPctDiverge() <= $partner->getPctDiverge() ) {
    $subBases =
        $thisLength * $this->getPctDiverge() + ( $partnerLength - $qo ) *
        $partner->getPctDiverge();
    $delBases = $thisLength * $this->getPctDelete() + ( $partnerLength - $qo ) *
        $partner->getPctDelete();
    $insBases = $thisLength * $this->getPctInsert() + ( $partnerLength - $qo ) *
        $partner->getPctInsert();
  }
  else {
    $subBases =
        $partnerLength * $partner->getPctDiverge() + ( $thisLength - $qo ) *
        $this->getPctDiverge();
    $delBases =
        $partnerLength * $partner->getPctDelete() + ( $thisLength - $qo ) *
        $this->getPctDelete();
    $insBases =
        $partnerLength * $partner->getPctInsert() + ( $thisLength - $qo ) *
        $this->getPctInsert();
  }

  # Include the query or consensus gap in $PctInsert or $PctDelete
  # Expressed in percent (not bases)
  if ( $qo <= 10 && $method ne "noQueryOverlapAdj" ) {
    if ( $co - $qo > 0 ) {
      $insBases += ( $co - $qo ) * 100;
    }
    else {
      $delBases -= ( $co - $qo ) * 100;

    }
  }

  my $totalLength = $thisLength + $partnerLength - $qo;
  $totalLength = 1 if $totalLength < 1;

  print "getAdjustedSubstLevel(): subBases="
      . ( $subBases / 100 )
      . ", delBases="
      . ( $delBases / 100 )
      . ", insBases= "
      . ( $insBases / 100 )
      . ", totalLength = $totalLength\n"
      if ( $DEBUG );

  return (
           ( $subBases / 100 ),
           sprintf( "%0.2f", ( $subBases / $totalLength ) ),
           ( $delBases / 100 ),
           sprintf( "%0.2f", ( $delBases / $totalLength ) ),
           ( $insBases / 100 ),
           sprintf( "%0.2f", ( $insBases / $totalLength ) )
  );

}

# Adjust the current element's divergence characteristics based
# on the overlap between it and a past element ( j ). If $ignorediv
# is set simply set the current to the past element -- note this
# is possibly not needed anymore.
sub adjustSubstLevel {
  my $this    = shift;
  my $partner = shift;
  my $method  = shift;

  my $DEBUG = 0;
  if ( $DEBUG ) {
    print "adjustSubstLevel( $this, $partner, $method ): Entered\n";
    $this->print();
    $partner->print();
  }

  warn "adjustSubstLevel(): Unknown method $method"
      if ( defined $method && $method !~ /noQueryOverlapAdj/ );

  my ( $subBases, $subPct, $delBases, $delPct, $insBases, $insPct ) =
      getAdjustedSubstLevel( $this, $partner, $method );

  $this->setPctDiverge( $subPct );

  # TODO...we should also adjust Kimura here
  $this->setPctDelete( $delPct );
  $this->setPctInsert( $insPct );

}

sub comparePositionOrder {
  my $this       = shift;
  my $refElement = shift;

  # Am I on the left or right of the reference?
  return ( $this->getQueryName() ) cmp( $refElement->getQueryName() )
      || ( $this->getQueryStart() ) <=>     ( $refElement->getQueryStart() )
      || ( $refElement->getQueryEnd() ) <=> ( $this->getQueryEnd() )
      || ( $refElement->getScore() ) <=>    ( $this->getScore() );
}

#
# Tests if this current element is joined
# so that it "contains" the refElement.
#
#
#      +-------------------------------------+
#      |                                     |
#   -------  ..  --refElement---  ..   ---this----
#
#  or
#
#      +-------------------------------------+
#      |                                     |
#   --this---  ..  --refElement---  ..   --------
#
sub containsElement {
  my $this       = shift;
  my $refElement = shift;

  # Am I on the left or right of the reference?
  if ( $this->comparePositionOrder( $refElement ) > 0 ) {

    # I am on the right
    my $leftLink = $this->getLeftLinkedHit();
    return 0 if ( $leftLink == undef || $leftLink == $this );
    return 1 if ( $leftLink->comparePositionOrder( $refElement ) < 0 );
  }
  else {

    # I am on the left
    my $rightLink = $this->getRightLinkedHit();
    return 0 if ( $rightLink == undef || $rightLink == $this );
    return 1 if ( $rightLink->comparePositionOrder( $refElement ) > 0 );
  }
}

sub getConsensusGap {
  my $this    = shift;
  my $partner = shift;
  return ( -$this->getConsensusOverlap( $partner ) );
}

sub getConsensusOverlap {
  my $this    = shift;
  my $partner = shift;

  warn "getConsensusOverlap() mismatched orientations!\n"
      if ( $this->getOrientation() != $partner->getOrientation() );

  my $left  = $this;
  my $right = $partner;
  if ( $left->comparePositionOrder( $right ) > 0 ) {

    # I am on the right
    $left  = $partner;
    $right = $this;
  }

  if ( $this->getOrientation() eq "C" ) {

    #                      0
    # <--left--, <--right---
    #
    return ( $right->getSubjEnd() - $left->getSubjStart() + 1 );
  }
  else {

    # 0
    # --left-->, --right-->
    #
    return ( $left->getSubjEnd() - $right->getSubjStart() + 1 );
  }
}

sub getQueryGap {
  my $this    = shift;
  my $partner = shift;
  return ( -$this->getQueryOverlap( $partner ) );
}

sub getQueryOverlap {
  my $this    = shift;
  my $partner = shift;

  my $sign = 1;
  $sign = -1
      if (    $this->getQueryStart() > $partner->getQueryEnd()
           || $partner->getQueryStart > $this->getQueryEnd() );

  my $left = $this->getQueryStart();
  $left = $partner->getQueryStart()
      if ( $partner->getQueryStart() > $left );
  my $right = $this->getQueryEnd();
  $right = $partner->getQueryEnd()
      if ( $partner->getQueryEnd() < $right );
  return ( abs( $right - $left + 1 ) * $sign );

}

1;
