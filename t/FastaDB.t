use strict;
use SeqDBI;
use FastaDB;
use Data::Dumper;
use Test::More;

#------------------------------------------------------------------------##
# Tests
#------------------------------------------------------------------------##
##
##
##
my $seqdb = FastaDB->new();
ok( defined $seqdb, "Testing object creation ->new()" );
$seqdb = undef;


#
# A place for our stuff
#
my $wrkDir = "wrkdir";
if (! -d $wrkDir ) { 
  system( "mkdir $wrkDir" );
}


##
##
##
system( "cp seqs/fastaformat/simple-good-small.fa $wrkDir" );
$seqdb = FastaDB->new( fileName=>"$wrkDir/simple-good-small.fa", 
                          openMode=>SeqDBI::ReadWrite );
ok( defined $seqdb, "Testing object creation - ReadWrite mode" );
                     
##
##
##
my $seq = $seqdb->getSubstr( "seq1", 0 , 3 );
my $desc = $seqdb->getDescription( "seq1" );
ok( ( $desc eq "Good format, small, one line sequence " .
               "DNA bases-only (21 bp)" &&
      $seq eq "ACC" ), "Testing sequence extraction 3-bp begin" );

##
##
##
my $seq1 = $seqdb->getSubstr( "seq1", 18 );
my $seq2 = $seqdb->getSubstr( "seq1", 18, 3 );
ok( $seq1 eq "TGA" && $seq2 eq "TGA",
    "Testing sequence extraction 3-bp end" );

##
##
##
$seq = $seqdb->getSubstr( "seq1", 4, 4 );
ok( $seq eq "TGTG",
    "Testing sequence extraction 4-bp middle" );

##
##
##
$seq1 = $seqdb->getSequence( "seq1" );
$seq2 = $seqdb->getSubstr( "seq1", 0, 21 ); 
ok( $seq1 eq $seq2 && $seq1 eq "ACCGTGTGTAGCTGTCGATGA", 
    "Testing sequence extraction full first seq" );

##
##
##
$seq1 = $seqdb->getSequence( "seq2" );
$seq2 = $seqdb->getSubstr("seq2", 0, 60 );
ok( $seq1 eq $seq2 &&
    $seq1 eq "TCCGACAGCATCGTACGTAGGCGAGTTATTCGGGAATTATACGTGATGCGCAGCTGATCG",
    "Testing sequence extraction full second seq" );

##
##
##
$seq = $seqdb->getSubstr( "seq6", 185, 120 );
ok( $seq eq "TATATCTTCAAATTTTTTTAAGATGAATGTTTAAAGATGTTTTTGACAGAATGAATACCTTTCAAAAACTATCCCTTTTGGAAGTCAGTCCGTAACTTTCGTTACGATTCCTTATTTTCA", "Testing sequence extraction partial sixth seq");
$seqdb = undef;

##
##
##
system( "cp seqs/fastaformat/simple-good-large.fa $wrkDir" );
$seqdb = FastaDB->new( fileName=>"$wrkDir/simple-good-large.fa",
                       openMode=>SeqDBI::ReadWrite );
ok( defined $seqdb, "Testing object creation - large file" );

##
##
##
$seq = $seqdb->getSubstr( "seq1", 35838, 20 );
ok( $seq eq "TGCTAGCCGTCACGTAAGTT", "Testing sequence extraction large file" );

##
##
##
$seqdb->setSubstr( "seq1", 605500, 140, "X" x 140 );
$seq = $seqdb->getSubstr( "seq1", 605635, 5 );
ok( $seq eq "XXXXX", "Testing sequence substitution large file" );



##
##
##
system( "cp seqs/fastaformat/simple-good-small.fa $wrkDir" );
$seqdb = FastaDB->new( fileName=>"$wrkDir/simple-good-small.fa",
                       openMode=>SeqDBI::ReadWrite );
my @ids = $seqdb->getIDs();
my @lengths = $seqdb->getSeqLengths();
ok ( defined $seqdb &&
     @ids ==  6 &&  @lengths == 6 &&
     $lengths[0] == 21 &&
     $lengths[1] == 60 &&
     $lengths[2] == 83 &&
     $lengths[3] == 120 &&
     $lengths[4] == 1330 &&
     $lengths[5] == 21741,
     "Testing counting sequence lengths and sequence count" );


##
## 
##
my $val1;
my $val2;
my $val3;
$val1 = $seqdb->getGCLength( "seq1" );
$val2 = $seqdb->getGCLength( "seq2" );
$val3 = $seqdb->getGCLength( "seq3" );
ok ( $val1 == 11  &&
  $val2 == 31 &&
  $val3 == 47, "Testing getGCLength" );

##
## 
##
$val1 = $seqdb->getSubtLength( "seq1" );
$val2 = $seqdb->getSubtLength( "seq2" );
$val3 = $seqdb->getSubtLength( "seq3" );
ok( 
  $val1 == 21 &&
  $val2 == 60 &&
  $val3 == 83, "Testing getSubtLength" );

##
##
##
$val1 = $seqdb->setSubstr( "seq1", 0, 4, "XXXX" );
$seq1 = $seqdb->getSequence( "seq1" );
my $subtVal1 = $seqdb->getSubtLength( "seq1" );
my $gcVal1 = $seqdb->getGCLength( "seq1" );

$seqdb->setSubstr( "seq2", 17, 6, "XXXXXXXXXX" ); #TAGGCG
$seq2 = $seqdb->getSequence( "seq2" );
my $subtVal2 = $seqdb->getSubtLength( "seq2" );
my $gcVal2 = $seqdb->getGCLength( "seq2" );

$seqdb->setSubstr( "seq3", 17, 23, "" );
my $seq3 = $seqdb->getSequence( "seq3" );
my $subtVal3 = $seqdb->getSubtLength( "seq3" );
my $gcVal3 = $seqdb->getGCLength( "seq3" );

$seqdb->setSubstr( "seq4", 119, 1, "X" );
my $seq4 = $seqdb->getSequence( "seq4" );
ok( 
        $seq1 eq "XXXXTGTGTAGCTGTCGATGA" && $val1 eq "ACCG" &&
        $subtVal1 == 17 && $gcVal1 == 8 &&
        $seq2 eq "TCCGACAGCATCGTACGXXXXXXXXXXAGTTATTCGGGAATTATACGTGATGCGCAGCTGATCG" &&
        $subtVal2 == 54 && $gcVal2 == 27 &&
        $seq3 eq "GCACGCGATTTTTGGGGGCATGATCGATCGTGCGGACGGTAGATGATGACGCATGCTAGG" && 
        $seq4 eq "CGATTTTGTAGCTACGACGCGCGGGCGAGCGCTACGTATATATCAGACGATCGTGCGACTACG" .
                 "ACTGTTGTGATGACTTATATGCGGCGCGGCGGCATCGATGACTACGGCTACGGGACX",
        "Testing numerous translations" );


##
##
##
$val1 = $seqdb->setSubstr( "seq1", 0, $seqdb->getSeqLength( "seq1"), "XXXX" );
$seq1 = $seqdb->getSequence( "seq1" );
$subtVal1 = $seqdb->getSubtLength( "seq1" );
$gcVal1 = $seqdb->getGCLength( "seq1" );
ok( 
  $seq1 eq "XXXX" && $subtVal1 == 0 &&
  $gcVal1 == 0, "Testing setSubstr -- completly masked" );


##
##
##
system( "cp seqs/fastaformat/missingids.fa $wrkDir" );
$seqdb = FastaDB->new( fileName=>"$wrkDir/missingids.fa",
                       openMode=>SeqDBI::ReadWrite );
@ids = $seqdb->getIDs();
@lengths = $seqdb->getSeqLengths();
ok ( defined $seqdb &&
     @ids ==  3 &&  @lengths == 3 &&
     $ids[0] eq "UnnamedSeq" &&
     $ids[1] eq "UnnamedSeq_1" &&
     $ids[2] eq "UnnamedSeq_2",
     "Testing missing ids"
   );

##
##
##
system( "cp seqs/fastaformat/emptyrecord.fa $wrkDir" );
$seqdb = FastaDB->new( fileName=>"$wrkDir/emptyrecord.fa",
                       openMode=>SeqDBI::ReadWrite );
@ids = $seqdb->getIDs();
@lengths = $seqdb->getSeqLengths();
ok( defined $seqdb &&
    @ids ==  6 &&  @lengths == 6 &&
    $ids[0] eq "seq1" &&
    $seqdb->getSubtLength( "seq2" ) == 0 &&
    $ids[5] eq "seq6" ,
    "Empty Records Test" );


##
##
##
$seqdb->setSubstr( "seq4", 0, 4, "" );
ok( defined $seqdb &&
    $seqdb->getSeqLength( "seq4" ) == "0",
    "Null a complete record" );
$seqdb = undef;


##
##
##
system( "cp seqs/fastaformat/simple-good-small.fa $wrkDir" );
$seqdb = FastaDB->new( fileName=>"$wrkDir/simple-good-small.fa",
                          openMode=>SeqDBI::ReadWrite );
my @seqIDs = $seqdb->getIDs();
my $failure = 0;
foreach my $seqID ( @seqIDs ) {
  print "\#   - Testing sequence $seqID\n";   
  for ( my $i = 0; $i < 2; $i++ ) {
    print "\#       - Test #$i\n";
    unless ( &ShiftTest( $seqdb, $seqID ) == 1 ) {
      print "\#          Sequence $seqID failed the ShiftTest\n";
      $failure = 1;
    }
    unless ( &SwapTest( $seqdb, $seqID ) == 1 ) {
      print "\#          Sequence $seqID failed the SwapTest\n";
      $failure = 1;
    }
  }
}
ok( $failure == 0, "Swap and Shift testing" );






$seqdb = undef;

######################## READ ONLY TESTS ############################

##
##
##
system( "cp seqs/fastaformat/simple-good-small.fa $wrkDir" );
$seqdb = FastaDB->new( fileName=>"$wrkDir/simple-good-small.fa",
                          openMode=>SeqDBI::ReadOnly );
ok( defined $seqdb, "Testing object creation - ReadOnly" );
                     
##
##
##
$seq = $seqdb->getSubstr( "seq1", 0 , 3 );
$desc = $seqdb->getDescription( "seq1" );
ok( $desc eq
    "Good format, small, one line sequence " .
    "DNA bases-only (21 bp)" &&
    $seq eq "ACC",
    "Testing sequence extraction 3-bp begin" );

##
##
##
$seq1 = $seqdb->getSubstr( "seq1", 18 );
$seq2 = $seqdb->getSubstr( "seq1", 18, 3 );
ok( $seq1 eq "TGA" &&
    $seq2 eq "TGA", 
    "Testing sequence extraction 3-bp end" );

##
##
##
$seq = $seqdb->getSubstr( "seq1", 4, 4 );
ok( $seq eq "TGTG",
    "Testing sequence extraction 4-bp middle" );

##
##
##
$seq1 = $seqdb->getSequence( "seq1" );
$seq2 = $seqdb->getSubstr( "seq1", 0, 21 ); 
ok( $seq1 eq $seq2 && 
    $seq1 eq "ACCGTGTGTAGCTGTCGATGA",
    "Testing sequence extraction full first seq" );

##
##
##
$seq1 = $seqdb->getSequence( "seq2" );
$seq2 = $seqdb->getSubstr("seq2", 0, 60 );
ok( $seq1 eq $seq2 &&
    $seq1 eq "TCCGACAGCATCGTACGTAGGCGAGTTATTCGGGAATTATACGTGATGCGCAGCTGATCG",
    "Testing sequence extraction full second seq" );

##
##
##
$seq = $seqdb->getSubstr( "seq6", 185, 120 );
ok( $seq eq "TATATCTTCAAATTTTTTTAAGATGAATGTTTAAAGATGTTTTTGACAGAATGAATACCTTTCAAAAACTATCCCTTTTGGAAGTCAGTCCGTAACTTTCGTTACGATTCCTTATTTTCA", "Testing sequence extraction partial sixth seq" );
$seqdb = undef;

##
##
##
system( "cp seqs/fastaformat/simple-good-large.fa $wrkDir" );
$seqdb = FastaDB->new( fileName=>"$wrkDir/simple-good-large.fa",
                       openMode=>SeqDBI::ReadOnly );
ok( defined $seqdb, "Testing object creation - large ReadOnly example" );

##
##
##
$seq = $seqdb->getSubstr( "seq1", 35838, 20 );
ok( $seq eq "TGCTAGCCGTCACGTAAGTT", "Testing sequence extraction large file" );

##
##
##
system( "cp seqs/fastaformat/simple-good-small.fa $wrkDir" );
$seqdb = FastaDB->new( fileName=>"$wrkDir/simple-good-small.fa",
                       openMode=>SeqDBI::ReadOnly );
@ids = $seqdb->getIDs();
@lengths = $seqdb->getSeqLengths();
ok( defined $seqdb &&
    @ids ==  6 &&  @lengths == 6 &&
    $lengths[0] == 21 &&
    $lengths[1] == 60 &&
    $lengths[2] == 83 &&
    $lengths[3] == 120 &&
    $lengths[4] == 1330 &&
    $lengths[5] == 21741,
    "Testing short ReadOnly sequence count and lengths"
             );
$seqdb = undef;


done_testing();

###

sub SwapTest {
  my $seqdb = shift;
  my $seqID = shift;
                                                                                
  my $origSeq = $seqdb->getSequence( $seqID );
  my $seqLen = $seqdb->getSeqLength( $seqID );
  my $chunkSize = 0;
  my $startPos = 0;
  my $index = 0;
  my @chunks = ();
  do {
    $chunkSize = sprintf("%1.0f",rand($seqLen - 1)) + 1;
    #print "Pushing $startPos , $chunkSize\n";
    push @chunks, { 'index' => $index,
                    'start' => $startPos,
                    'length' => $chunkSize };
    $index++;
    $startPos += $chunkSize;
    $seqLen -= $chunkSize;
  } while ( $seqLen > 1 );
  if ( $seqLen == 1 ) {
    #print "Pushing $startPos, 1\n";
    push @chunks, { 'index' => $index,
                    'start' => $startPos,
                    'length' => 1 };
  }
                                                                                
  &fisherYatesShuffle( \@chunks );
                                                                                
  foreach my $chunk ( @chunks ) {
    #print "Obtaining sequence for " . $chunk->{'start'} . " ( " . $chunk->{'length'} . " )\n";
    $chunk->{'seq'} = $seqdb->setSubstr( $seqID, $chunk->{'start'},
                                         $chunk->{'length'}, "x" );
    $seqdb->setSubstr( $seqID, $chunk->{'start'},
                       1, "x" x $chunk->{'length'} );
  }
                                                                                
  &fisherYatesShuffle( \@chunks );
                                                                                
  foreach my $chunk ( @chunks ) {
    #print "Setting sequence for " . $chunk->{'start'} . " ( " . $chunk->{'length'} . " )\n";
    $seqdb->setSubstr( $seqID, $chunk->{'start'},
                       $chunk->{'length'}, $chunk->{'seq'}  );
  }
                                                                                
  my $newSeq = $seqdb->getSequence( $seqID );
  if ( $newSeq ne $origSeq ) {
    return 0;
  }
  return 1;
}
                                                                                
sub fisherYatesShuffle {
  my $array = shift;
  my $i;
  for( $i = @$array-1; $i >= 0; --$i ){
    my $j = int rand ( $i + 1 );
    next if $i == $j;
    @$array[ $i, $j ] = @$array[ $j, $i ];
  }
}
                                                                                
sub ShiftTest {
  my $seqdb = shift;
  my $seqID = shift;
                                                                                
  my $origSeq = $seqdb->getSequence( $seqID );
  my $seqLen = $seqdb->getSeqLength( $seqID );
  my $chunkSize = 0;
  my @chunks = ();
  do {
    $chunkSize = sprintf("%1.0f",rand($seqLen - 1)) + 1;
    #print "Pushing $chunkSize\n";
    push @chunks, $chunkSize;
    $seqLen -= $chunkSize;
  } while ( $seqLen > 1 );
  if ( $seqLen == 1 ) {
    #print "Pushing 1\n";
    push @chunks, 1;
  }
  $seqLen = $seqdb->getSeqLength( $seqID );
  foreach $chunkSize ( @chunks ) {
    my $seq = $seqdb->setSubstr( $seqID, 0, $chunkSize, "" );
    $seqdb->setSubstr( $seqID, $seqLen - $chunkSize, 0, $seq );
  }
  my $newSeq = $seqdb->getSequence( $seqID );
  if ( $newSeq ne $origSeq ) {
    return 0;
  }
  return 1;
}


