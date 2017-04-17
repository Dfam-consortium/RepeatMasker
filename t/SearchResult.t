use strict;
use CrossmatchSearchEngine;
use SearchResult;
use SearchResultCollection;
use Data::Dumper;
use Test::More;


my $wrkDir = "wrkdir";
if (! -d $wrkDir ) { 
  system( "mkdir $wrkDir" );
}


#------------------------------------------------------------------------##
# Tests
#------------------------------------------------------------------------##

# Crossmatch 1.080812
# Unadjusted scoring -raw run
#/usr/local/bin/cross_match -raw -alignments -gap_init -30 -ins_gap_ext -6 -del_gap_ext -5 -minmatch 8 -minscore 225 -bandwidth 14 -masklevel 101 -matrix /home/rhubley/projects/RepeatMasker/Matrices/crossmatch/20p43g.matrix chr7-25kb.fa /home/rhubley/projects/RepeatMasker/Libraries/20130422/homo_sapiens/longlib
# 353 28.71 0.00 1.98  seq      969  1069 (23931)  C L2-1_AMi#LINE/L2   (680)   555   457
# 
#C seq                  1069 ACGTCTCCCCTTGGATGTCCTGCCAAGACCTCAAACTCAACATGTCTTAA 1020
#                              i     i v        iv   ivv v     i         v i    
#  L2-1_AMi#LINE/L       457 ACATCTCCTCGTGGATGTCTAGCCGTCAGCTCAAGCTCAACATGGCCTAA 506
#  
#C seq                  1019 AACCAAACTCATCTTCTTCCTCCTAAACCTCTGCCACCTACCTCCTCCTC 970
#                               vi i   v vv      i  i  i   --  iv  i iii   ii   
#  L2-1_AMi#LINE/L       507 AACAGAGCTCTTAATCTTCCCCCCAAGCCT--GCTCCCCATTCCCTTTTC 554
#
#C seq                   969 T 969
#
#  L2-1_AMi#LINE/L       555 T 555
#
#Transitions / transversions = 1.64 (18 / 11); Gap_init rate = 0.01 (1 / 101), avg. gap size = 2.00 (2 / 1)
# 
# Again with complexity adjustment
#
# 253 28.71 0.00 1.98  seq      969  1069 (23931)  C L2-1_AMi#LINE/L2   (680)   555   457
#
#C seq                  1069 ACGTCTCCCCTTGGATGTCCTGCCAAGACCTCAAACTCAACATGTCTTAA 1020
#                              i     i v        iv   ivv v     i         v i    
#  L2-1_AMi#LINE/L       457 ACATCTCCTCGTGGATGTCTAGCCGTCAGCTCAAGCTCAACATGGCCTAA 506
#
#C seq                  1019 AACCAAACTCATCTTCTTCCTCCTAAACCTCTGCCACCTACCTCCTCCTC 970
#                               vi i   v vv      i  i  i   --  iv  i iii   ii   
#  L2-1_AMi#LINE/L       507 AACAGAGCTCTTAATCTTCCCCCCAAGCCT--GCTCCCCATTCCCTTTTC 554
#
#C seq                   969 T 969
#                              
#  L2-1_AMi#LINE/L       555 T 555
#
my $searchResult = new SearchResult( 
      score => 353, pctDiverge => 28.71,  pctDelete => 0.00, pctInsert => 1.98,
      queryName => "seq",  queryStart => 969,  queryEnd => 1069,  queryRemaining => 23931, 
      orientation => "C",
      subjName => "L2-1_AMi#LINE/L2", subjStart => 457, subjEnd => 555, subjRemaining => 680,
      id => undef,  lineageId => undef,  overlap => undef,
      matrixName => "20p43g.matrix", queryString => "ACGTCTCCCCTTGGATGTCCTGCCAAGACCTCAAACTCAACATGTCTTAAAACCAAACTCATCTTCTTCCTCCTAAACCTCTGCCACCTACCTCCTCCTCT", subjString => "ACATCTCCTCGTGGATGTCTAGCCGTCAGCTCAAGCTCAACATGGCCTAAAACAGAGCTCTTAATCTTCCCCCCAAGCCT--GCTCCCCATTCCCTTTTCT" );

# Basic object accessor test
ok ( $searchResult->getScore() == 353, "Basic accessor test" ) or 
  diag( "Object accessor returned incorrect result! Returned " . $searchResult->getScore() . " instead of 353." );


# Rescore alignment ( raw )
use Matrix;
my $matrix = Matrix->new( fileName => "/home/rhubley/projects/RepeatMasker/Matrices/crossmatch/20p43g.matrix");

my ( $score, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions )  = $searchResult->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -30,
                             gapExtPenalty => -6
                            );

ok ( $score == $searchResult->getScore(), "Basic rescoring (raw) test" ) or
  diag("Failed to recalculate the original score " . $searchResult->getScore() . ", returned $score\n");
ok ( $divergence - 0.01 < 38.4766672846808 && $divergence + 0.01 > 38.4766672846808,
     "Divergence test" ) or
  diag("Failed to calculate the standard Kimura divergence!  Returned $divergence, expected 38.4766672846808\n");

# Rescore alignment ( complexity adjusted )
( $score, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions )  = $searchResult->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -30,
                             gapExtPenalty => -6,
                             complexityAdjust => 1
                            );

ok ( $score == 253, "Complexity adjusted scoring test" ) or
  diag("Failed to calculate complexity adjusted score. Returned $score expected 253\n");
ok ( $transitions == 18, "Transition calculation test" ) or
  diag("Transition calculation failed! Returned $transitions expected 18.\n");
ok ( $transversions == 11, "Transversion calculation test" ) or
  diag("Transversion calculation failed! Returned $transversions expected 11\n");
ok ( $cpgsites == 2, "CpG calculation test" ) or
  diag("Calculation of CpG sites failed!  Returned $cpgsites, expected 2\n");


# Artificial example 
# Raw score:
#  seq                  969 AGTTTATTTAT--AATT-T 984
#                            vv   i--i --i  i- 
#  L2-1_AMi#LINE        457 ACGTTAC--GTTCGATCGT 473
#
#  Manual analysis:
#     AGTTTATTTAT--AATT-TGGC
#      xx   x  x  xx  xx
#     ACGTTAC--GTTCGATCGTGGC
#           A/A = 3 * 9   =  27 
#           C/C = 3 * 10  =  30
#           T/T = 5 * 9   =  45 
#           G/G = 4 * 10  =  40
#           C/G = 1 * -15 = -15
#           G/T = 1 * -15 = -15
#      gap_init = 3 * -30 = -90
#       gap_ext = 2 * -6  = -12
#    ---------------------------
#                            10 
#    4 CpG Sites
#    2 Tansversions
#    1.2 Transitions ( using divCpGMod )
#    14 Well Characterized Bases    
#    27.2765174726037 % div ( Kimura using divCpGMod )
#                  
$searchResult = new SearchResult( 
      score => 30, pctDiverge => 0.00,  pctDelete => 0.00, pctInsert => 0.00,
      queryName => "seq",  queryStart => 969,  queryEnd => 987,  queryRemaining => 0, 
      orientation => "",
      subjName => "L2-1_AMi#LINE/L2", subjStart => 457, subjEnd => 476, subjRemaining => 0,
      id => undef,  lineageId => undef,  overlap => undef,
      matrixName => "20p43g.matrix", queryString => "AGTTTATTTAT--AATT-TGGC", subjString => "ACGTTAC--GTTCGATCGTGGC" );
( $score, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions )  = $searchResult->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -30,
                             gapExtPenalty => -6,
                             scoreCpGMod => 1,
                             divCpGMod => 1
                            );

my ( $div, $transi, $transv, $wellCharBases, $numCpGs )
        = $searchResult->calcKimuraDivergence( divCpGMod => 1 );

ok ( $score == 10 ) or
  diag("Failed to calculate raw score. Returned $score expected 10\n");
ok ( $transitions == $transi && ( $transi > (1.2 - 0.000001) && $transi < (1.2 + 0.000001) ) ) or
  diag("Failed to calculate transitions correctly using divCpGMod. Returned $transi expected 1.2\n");
ok ( $transversions == $transv && $transv == 2  ) or
  diag("Failed to calculate transversions correctly using divCpGMod. Returned $transv expected 2\n");
ok ( $wellCharBases == $well_characterized_bases && $wellCharBases == 17  ) or
  diag("Failed to calculate well characterized bases correctly using divCpGMod. Returned $wellCharBases expected 17\n");
#print "Num CpGs = $numCpGs, transi/transv = $transi/$transv, kimura = $div wcb = $wellCharBases other divergence = $divergence and $well_characterized_bases i/v = $transitions/$transversions\n";

my $qSeq = $searchResult->getQueryString();
ok ( $qSeq eq "AGTTTATTTAT--AATT-TGGC", "Query String Retrieval" );
substr( $qSeq, 0, 1 ) = lc( substr( $qSeq, 0, 1 ) );
$searchResult->setQueryString( $qSeq );
my $qSeq2 = $searchResult->getQueryString();
ok ( $qSeq2 eq "aGTTTATTTAT--AATT-TGGC", "Query String After Modification" ) or
  diag("Returned: $qSeq2 rather than aGTTTATTTAT--AATT-TGGC\n");
my $printStr = $searchResult->toStringFormatted( SearchResult::AlignWithQuerySeq );
ok ( $printStr =~ /aGTTTATTTAT--AATT-TGGC/, "Alignment toStringFormatted Contains Modification" ) or
  diag("Returned:\n$printStr\n" );

# Another Artificial example 
#      30 0 0 0 seq 969 984 (0) L2-1_AMi#LINE/L2 457 473 (0)
#
#         seq                  969 ATAG--------AA 974
#                                     v--------i  
#         L2-1_AMi#LINE        457 ATACGTTCGTTCGA 470
#
#  3 CpG Sites
#  0.1 Transitions
#  1 Transversion
#
$searchResult = new SearchResult( 
      score => 30, pctDiverge => 0.00,  pctDelete => 0.00, pctInsert => 0.00,
      queryName => "seq",  queryStart => 969,  queryEnd => 974,  queryRemaining => 0, 
      orientation => "",
      subjName => "L2-1_AMi#LINE/L2", subjStart => 457, subjEnd => 470, subjRemaining => 0,
      id => undef,  lineageId => undef,  overlap => undef,
      matrixName => "20p43g.matrix", queryString => "ATAG--------AA", subjString => "ATACGTTCGTTCGA" );
( $score, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions )  = $searchResult->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -30,
                             gapExtPenalty => -6,
                             scoreCpGMod => 1,
                             divCpGMod => 1
                            );

( $div, $transi, $transv, $wellCharBases, $numCpGs )
        = $searchResult->calcKimuraDivergence( divCpGMod => 1 );

ok ( $transitions == $transi && ( $transi > (0.1 - 0.000001) && $transi < (0.1 + 0.000001) ) ) or
  diag("Failed to calculate transitions correctly using divCpGMod. Returned $transi expected 0.1\n");
ok ( $transversions == $transv && $transv == 1  ) or
  diag("Failed to calculate transversions correctly using divCpGMod. Returned $transv expected 1\n");
ok ( $wellCharBases == $well_characterized_bases && $wellCharBases == 6  ) or
  diag("Failed to calculate well characterized bases correctly using divCpGMod. Returned $wellCharBases expected 6\n");
#

( $score, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions )  = $searchResult->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -30,
                             gapExtPenalty => -6
                            );

( $div, $transi, $transv, $wellCharBases, $numCpGs )
        = $searchResult->calcKimuraDivergence();

ok ( $transitions == $transi && $transi == 1) or
  diag("Failed to calculate transitions correctly. Returned $transi expected 1\n");
ok ( $divergence == $div &&  ( $div > (44.7939867307014 - 0.0001) && $div < (44.7939867307014 + 0.0001) ) ) or
  diag("Failed to calculate divergence correctly. Returned $div expected 44.7939867307014\n");

#print "Num CpGs = $numCpGs, transi/transv = $transi/$transv, kimura = $div wcb = $wellCharBases other divergence = $divergence and $well_characterized_bases i/v = $transitions/$transversions cpg=$cpgsites\n";


# 
# Fragmented Alignment Artefact
#   - RepeatMasker can generate these weird alignments during fragmentation.
#     We need to make sure we sensibly handle them in the rescoreAlignment 
#     routines.
#
#      300 10.0 10.0 0.0 seq 509 524 (300) L1_Mur2_orf2#LINE/L1 1771 1771 (100)
#
#       seq                  509 TGGAAGGAAGGAAGGA 524
#                                ----------------
#       L1_Mur2_orf2#       1771 ---------------- 1772
#
#    score = -30 + ( 15 * -6 ) = -120
#
$searchResult = new SearchResult( 
      score => 300, pctDiverge => 10.00,  pctDelete => 10.00, pctInsert => 10.00,
      queryName => "seq",  queryStart => 509,  queryEnd => 524,  queryRemaining => 300, 
      orientation => "",
      subjName => "L1_Mur2_orf2#LINE/L1", subjStart => 1771, subjEnd => 1771, subjRemaining => 100,
      id => undef,  lineageId => undef,  overlap => undef,
      matrixName => "20p43g.matrix", queryString => "TGGAAGGAAGGAAGGA", subjString => "----------------" );
( $score, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions )  = $searchResult->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -30,
                             gapExtPenalty => -6,
                             scoreCpGMod => 1,
                             divCpGMod => 1
                            );

( $div, $transi, $transv, $wellCharBases, $numCpGs )
        = $searchResult->calcKimuraDivergence( divCpGMod => 1 );

ok ( $divergence == 100 ) or
  diag("Failed to calculate divergence. Returned $divergence expected 100\n");
ok ( $percIns == 100 ) or
  diag("Failed to calculate percent insertion. Returned $percIns expected 100\n");
ok ( $percDel == 0 ) or
  diag("Failed to calculate percent deletion. Returned $percDel expected 0\n");
ok ( $score == -120 ) or
  diag("Failed to calculate raw score. Returned $score expected -120\n");
ok ( $transitions == $transi && ( $transi > (0 - 0.000001) && $transi < (0 + 0.000001) ) ) or
  diag("Failed to calculate transitions correctly using divCpGMod. Returned $transi expected 0\n");
ok ( $transversions == $transv && $transv == 0  ) or
  diag("Failed to calculate transversions correctly using divCpGMod. Returned $transv expected 0\n");
ok ( $wellCharBases == $well_characterized_bases && $wellCharBases == 0  ) or
  diag("Failed to calculate well characterized bases correctly using divCpGMod. Returned $wellCharBases expected 0\n");
#

#
#238 26.55 23.13 0.00 chr1 4162976 4163088 (191308883) C MLTR31FA_MM#LTR/ERVK (2) 574 428 1468 m_b3s501i17
#
#  chr1             4162976 TTGTGGGAAATATTAAAAAACGAACCACCAAG-TTCCCACATT-CCATCC 4163023
#                              i                i i   iv  i -     i ii - iiv  
#C MLTR31FA_MM#L        574 TTGCGGGAAATATTAAAAAATGGACCGGCAGGCTTCCCGCGCTACTGACC 525
#
#  chr1             4163024 AGCTACTCTTTGCCCA-------CCAAC-------TCTC----TGGTG-- 4163053
#                           i  v  ivviv v   -------   i -------    ----i  ? --
#C MLTR31FA_MM#L        524 GGCAACCGGCGGACCACATGGCTCCAGCCGGGTGATCTCAACTCGGYGGT 475
#
#  chr1             4163054 ------AGGCTGGAGCTCCCAAGTTCCCT------TCACTACCATGC 4163088
#                           ------i   i         i   i  i ------  i  i   i  
#C MLTR31FA_MM#L        474 CTCTCTGGGCCGGAGCTCCCGAGTCCCTTGGCGGATCGCTGCCACGC 428
#
#Matrix = 25p39g.matrix
#Transitions / transversions = 3.29 (23 / 7)
#Gap_init rate = 0.06 (7 / 112), avg. gap size = 4.86 (34 / 7)
$searchResult = new SearchResult( 
      score => 238, pctDiverge => 26.55,  pctDelete => 23.13, pctInsert => 0.00,
      queryName => "chr1",  queryStart => 4162976,  queryEnd => 4163088,  queryRemaining => 191308883, 
      orientation => "C",
      subjName => "MLTR31FA_MM#LTR/ERVK", subjStart => 428, subjEnd => 574, subjRemaining => 2,
      id => undef,  lineageId => undef,  overlap => undef,
      matrixName => "25p39g.matrix", queryString => "TTGTGGGAAATATTAAAAAACGAACCACCAAG-TTCCCACATT-CCATCCAGCTACTCTTTGCCCA-------CCAAC-------TCTC----TGGTG--------AGGCTGGAGCTCCCAAGTTCCCT------TCACTACCATGC", subjString => "TTGCGGGAAATATTAAAAAATGGACCGGCAGGCTTCCCGCGCTACTGACCGGCAACCGGCGGACCACATGGCTCCAGCCGGGTGATCTCAACTCGGYGGTCTCTCTGGGCCGGAGCTCCCGAGTCCCTTGGCGGATCGCTGCCACGC" );

$matrix = Matrix->new( fileName => "/home/rhubley/projects/RepeatMasker/Matrices/crossmatch/25p39g.matrix");
( $score, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions )  = $searchResult->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -27,
                             gapExtPenalty => -5,
                             complexityAdjust => 1
                            );

ok ( $score == $searchResult->getScore() ) or
  diag("Failed to recalculate the original score " . $searchResult->getScore() . ", returned $score\n");

#2331 25.75 8.90 6.51  Human     6606  8065 (62835)    L1MDa_5end#LINE/L1     1152  2646 (4)
#
#  Human                6606 ATAGCAGGGATGTTGGAATTATCAGACCAGGATTTTAAAACAACTATGAT 6655
#                              ii   i      i        v   v i         i  i i  ?   
#  L1MDa_5end#LINE      1152 ATGACAGAGATGTTAGAATTATCTGACAAAGATTTTAAAGCAGCCATNAT 1201
#  
#  Human                6656 TAATAT-------------ATTAAG-----GCCTGTAATG---GAAATAG 6684
#                            v  v  -------------    v -----  i  v  ii---    - i 
#  L1MDa_5end#LINE      1202 AAAAATGCTTCAATGAGCAATTACGAACACGCTTGAAACAAATGAAA-AA 1250
#  
#  Human                6685 TAGA------CGACATACAAGAACAAATGGGTAAGCAGAGGGATGGAGAT 6728
#                                ------ ii  v v  iv v   - vv v    i  ii  i  iiv 
#  L1MDa_5end#LINE      1251 TAGAAAGTCTCAGCAAAGAAATAGAAA-GTCTCAGCAAAGAAATAGAAGA 1299
#  
#  Human                6729 TCTAAGAAAGAATTTTAAT---AATGCTAAAGATAAAAAA-ACACTGCAA 6774
#                             v   ii     iiv-   ---   vi  i iv i     -   v --   
#  L1MDa_5end#LINE      1300 TATAAAGAAGAACCA-AATGGAAATTTTAGAACTGAAAAATACAAT--AA 1346
#  
#  Human                6775 CAAAAAGGAAGAATGCCTTTGATGGGCTCATTAGCACACTAGACGTGGCT 6824
#                             vi   vi  i iiv vv v          vi    v v i  v v i v 
#  L1MDa_5end#LINE      1347 CCGAAATAAAAAGCTCAGTGGATGGGCTCAACAGCAGAATGGAGGGGACA 1396
#
#  Human                6825 GAGGAATCTCTCAG---ATTG-AGGATATGGCAATAGAAACT-TCCAA-- 6867
#                                  vvvv    --- i  - i    vii         i -i    --
#  L1MDa_5end#LINE      1397 GAGGAAAGAATCAGTGAACTGGAAGATAGAACAATAGAAATTACCCAATC 1446
#
#  Human                6868 ----------AGAGAAAAGAGACTGAAAAAAAAAAAGGAAGAGAATTTCA 6907
#                            ----------        v                -v   v   iii
#  L1MDa_5end#LINE      1447 TGAACAACAGAGAGAAAATAGACTGAAAAAAAAAA-TGAACAGAGCCTCA 1495
#
#  Human                6908 AGAAC-TGTGGGACAACTACAAAAGGTGGAAAATAGGCATAATGGAAATA 6956
#                            i i  -        v iv       i vv  v  vviii v  v i i v
#  L1MDa_5end#LINE      1496 GGGACCTGTGGGACTATAACAAAAGATCTAACATTTATGTCATCGGAGTC 1545
#
#  Human                6957 TCAAAAGGAGGAGAAAAAGAGAAAGGAACAGAAGAAATACTTGAAGCAAT 7006
#                            i  i      ii  i      iiv  ii v   i  i    i    v
#  L1MDa_5end#LINE      1546 CCAGAAGGAGAGGAGAAAGAGGGCGGGGCTGAAAAAGTACTCGAAGAAAT 1595
#
#  Human                7007 AATGACTGAGAATTTCTCCAAATTAATGGCAGTTAAC-CAAACCT-CAGA 7054
#                                i    i  i   -       --     ivvi  -i      -
#  L1MDa_5end#LINE      1596 AATGGCTGAAAACTTC-CCAAATT--TGGCAAGAGACATAAACCTACAGA 1642
#
#  Human                7055 TCCAGGAAGCTCAGAGAACGCAAAGCAGAATAAATGCCAAAAACAACAAC 7104
#                             i  i      v  v    v v  i   i     i-     i -------
#  L1MDa_5end#LINE      1643 TTCAAGAAGCTGAGCGAACCCCAAACAGGATAAAC-CCAAAGA------- 1684
#
#  Human                7105 AAAGCCCTACACCTAGTCATATTATATTCAAACCACAGAAAATCAAAGAT 7154
#                            -  v  --     v  v  i  i   v      iv v     ii     i
#  L1MDa_5end#LINE      1685 -AATCC--ACACCAAGACACATCATAATCAAACTTCTGAAAACTAAAGAC 1731
#
#  Human                7155 GGAGAAG-------------TCAGGAAGTGGGAAAAATCCCATTACCTAT 7191
#                            ii    i-------------v   v?  v i    viv v v
#  L1MDa_5end#LINE      1732 AAAGAAAAAAAATCTTGAAAGCAGCNAGAGAGAAACGACACCTTACCTAT 1781
#
#  Human                7192 AGAGGAGGAAAGGATAAAACTACATTAGACTTCTCA---GAAACCATGCA 7238
#                              i   ii  viiv vi  iv   vii  i      ---         v
#  L1MDa_5end#LINE      1782 AGGGGAAAAACAATTCGAATGACAGCGGATTTCTCATCAGAAACCATGGA 1831
#
#  Human                7239 AGCAAGAGAGAGGAGAGTGAAATAAAGTGTT--AAG-ACAGAAAAAAAAA 7285
#                            i  v   ---    -    iv i  vi v  --   -i v    i    -
#  L1MDa_5end#LINE      1832 GGCCAGA---AGGA-AGTGGCACAACATTTTTCAAGTGCTGAAAGAAAA- 1876
#
#  Human                7286 GCAC---CAACCTAAAATTCT-------GTGAAATTGTCCTTCAAAAGTG 7325
#                             v  ---     i i   i  ------- i    v i       ii i
#  L1MDa_5end#LINE      1877 GAACTGTCAACCCAGAATCCTATATCCAGCGAAAATATCCTTCAGGAATG 1926
#
#  Human                7326 GAAGAGAAATAA-----TTACTCAGACAAACAAAAACTGAGAGAATTTGT 7370
#                            i i i     v -----  -      ii  vi      i
#  L1MDa_5end#LINE      1927 AAGGGGAAATCAAGACATT-CTCAGATGAAGGAAAACTAAGAGAATTTGT 1975
#
#  Human                7371 TGCC--------TGCCTTGAAAGAAATGCTAA-GACAGTTTTTTCAGAGA 7411
#                            i   -------- i  i i      vv     - iv    i i - i v
#  L1MDa_5end#LINE      1976 CGCCAGCAGACCTACCCTAAAAGAATGGCTAAAGGAAGTTCTCT-AAACA 2024
#
#  Human                7412 CAAGGAAAATGATATAGGTCAGAAACTCTGTTCTAT-----ATAAAGAA- 7455
#                            v  i i        v i -v  i  iii  vvvi  ----- vi     -
#  L1MDa_5end#LINE      2025 GAAAGGAAATGATAAAAG-AAGGAATCTTGGAACATCAGGAAGGAAGAAA 2073
#
#  Human                7456 GAGCATCAG--AGAAGAAATAAGTGAAGGTAAAATAGAACTTATTTTTCT 7503
#                              i  -   --  v i     v v v ii  v     ivi  vi  i  i
#  L1MDa_5end#LINE      2074 GAACA-CAGTAAGCAAAAATATGGGTAAATACAATAGGCTTTCCTTCTCC 2122
#
#  Human                7504 TATACTTAATTGATCTAACATAATAGTTTGCTCAAAATAATTATAGCAAC 7553
#                             ---   i v vv     -  v  -     ------ iii  iv     -
#  L1MDa_5end#LINE      2123 T---CTTGAGTTTTCTAA-ATTAT-GTTTG------ACGGTTGAAGCAA- 2160
#
#  Human                7554 AAAGTATTATATTATGGATACTTTGTGTGTCTATACATGCTTATATA--- 7600
#                               v   v i i i ------i vi    vi v  v   --   i  ---
#  L1MDa_5end#LINE      2161 AAATTATAACACTGT------CTGATGTGGTTCTAAATG--TATGTAGAG 2202
#
#  Human                7601 ----TGTTTA-------TGTATTAAGTAAAATGAATGAGTGATTAAAAGG 7639
#                            ---- i    ------- i     --    -  iiv   v --    i
#  L1MDa_5end#LINE      2203 GAAATATTTAAGACAATTATATTA--TAAA-TGGGGGAGGG--TAAAGGG 2247
#
#  Human                7640 ATG-AGAGGGAGGAATTAGGAATATTTTGTTATTATAAGGTACTTTCACC 7688
#                             ii- i       ---  i --------   vi    ------i     -
#  L1MDa_5end#LINE      2248 ACATAAAGGGAGG---TAAG--------GTTTCTATA------CTTCAC- 2279
#
#  Human                7689 ATCCCTGAAGTGTTATAGTGTTATTTGAAA---GTAGACTTGGATTAATT 7735
#                            -  ---i  v  v  v i  -- i--  v ---       vv   v i
#  L1MDa_5end#LINE      2280 -TC---AAACTGGTAAAATG--AC--GACACCAGTAGACTGTGATAAGTT 2321
#
#  Human                7736 GTAAACATATATTGCAAATTCTAGGGCAACCTCTAAAAAACATTTTTTAA 7785
#                            i iv i     v  i  vvi    i      v        -ii v vi
#  L1MDa_5end#LINE      2322 ATGTATATATAATGTAATACCTAGAGCAACCACTAAAAAA-GCTATACAA 2370
#
#  Human                7786 AAAGTATATATTAATTAATATGCTA-AGAAAGGAGAGAAAATGGAATAAT 7834
#                             i  -   i i v --  v ii   -   v -- ivv          vv
#  L1MDa_5end#LINE      2371 AGAG-ATACACTCA--AAAACACTATAGATA--AATCAAAATGGAATTCT 2415
#
#  Human                7835 ACAAAATGCTAAATTAAAGCCACAA-AAGGCAGAAAAAGAATAGAAGACA 7883
#                             v      i v  v   v-     i-       i    ----    i
#  L1MDa_5end#LINE      2416 AAAAAATGTTCAAGTAAC-CCACAGGAAGGCAGGAAAA----AGAAAACA 2460
#
#  Human                7884 AAAATAGGAACAA-AGAATAATGAATAGAAAATAGTAA-TAACTATGGTA 7931
#                            i i v v   v  -    i  --- i      i iv  -   v-    i
#  L1MDa_5end#LINE      2461 GAGAAATGAAAAACAGAACAA---ACAGAAAACAAAAAATAAA-ATGGCA 2506
#
#  Human                7932 GATATTAATCCAAATGTATCAATAATAATTTTAAACATCAGTTGTCTAA- 7980
#                              iv v ii  v  ii          v iv     ii v i v      -
#  L1MDa_5end#LINE      2507 GACTTAAGCCCTAACATATCAATAATTACATTAAATGTAAATGGTCTAAA 2556
#
#  Human                7981 TATACCAATGAAAAGACAGA--TTGTCAGAGTAGATTAAAAAAAAT---C 8025
#                              i      v          --   v      i          v  ---
#  L1MDa_5end#LINE      2557 TACACCAATTAAAAGACAGAGATTGGCAGAGTGGATTAAAAAACATGACC 2606
#  Human                8026 CAACTATATGTTGTATACAAGAAACCCACTTTTAATATAA 8065
#                                      i   v          i     iv
#  L1MDa_5end#LINE      2607 CAACTATATGCTGTCTACAAGAAACTCACTTCAAATATAA 2646
#
#Transitions / transversions = 1.34 (214 / 160); Gap_init rate = 0.05 (78 / 1460), avg. gap size = 2.88 (225
# / 78)
#
# And with complexity adjusted scoring:
# 1530 25.75 8.90 6.51  Human     6606  8065 (62835)    L1MDa_5end#LINE/L1     1152  2646 (4)
# ..
#Transitions / transversions = 1.34 (214 / 160); Gap_init rate = 0.05 (78 / 1460), avg. gap size = 2.88 (225 / 78)
#

$searchResult = new SearchResult( 
      score => 2331, pctDiverge => 25.75,  pctDelete => 8.90, pctInsert => 6.51,
      queryName => "Human",  queryStart => 6606, queryEnd => 8065,  queryRemaining => 62835, 
      orientation => "",
      subjName => "L1MDa_5end#LINE/L1", subjStart => 1152, subjEnd => 2646, subjRemaining => 4,
      id => undef,  lineageId => undef,  overlap => undef,
      matrixName => "20p43g.matrix",  
queryString =>
"ATAGCAGGGATGTTGGAATTATCAGACCAGGATTTTAAAACAACTATGAT" .
"TAATAT-------------ATTAAG-----GCCTGTAATG---GAAATAG" .
"TAGA------CGACATACAAGAACAAATGGGTAAGCAGAGGGATGGAGAT" .
"TCTAAGAAAGAATTTTAAT---AATGCTAAAGATAAAAAA-ACACTGCAA" .
"CAAAAAGGAAGAATGCCTTTGATGGGCTCATTAGCACACTAGACGTGGCT" .
"GAGGAATCTCTCAG---ATTG-AGGATATGGCAATAGAAACT-TCCAA--" .
"----------AGAGAAAAGAGACTGAAAAAAAAAAAGGAAGAGAATTTCA" .
"AGAAC-TGTGGGACAACTACAAAAGGTGGAAAATAGGCATAATGGAAATA" .
"TCAAAAGGAGGAGAAAAAGAGAAAGGAACAGAAGAAATACTTGAAGCAAT" .
"AATGACTGAGAATTTCTCCAAATTAATGGCAGTTAAC-CAAACCT-CAGA" .
"TCCAGGAAGCTCAGAGAACGCAAAGCAGAATAAATGCCAAAAACAACAAC" .
"AAAGCCCTACACCTAGTCATATTATATTCAAACCACAGAAAATCAAAGAT" .
"GGAGAAG-------------TCAGGAAGTGGGAAAAATCCCATTACCTAT" .
"AGAGGAGGAAAGGATAAAACTACATTAGACTTCTCA---GAAACCATGCA" .
"AGCAAGAGAGAGGAGAGTGAAATAAAGTGTT--AAG-ACAGAAAAAAAAA" .
"GCAC---CAACCTAAAATTCT-------GTGAAATTGTCCTTCAAAAGTG" .
"GAAGAGAAATAA-----TTACTCAGACAAACAAAAACTGAGAGAATTTGT" .
"TGCC--------TGCCTTGAAAGAAATGCTAA-GACAGTTTTTTCAGAGA" .
"CAAGGAAAATGATATAGGTCAGAAACTCTGTTCTAT-----ATAAAGAA-" .
"GAGCATCAG--AGAAGAAATAAGTGAAGGTAAAATAGAACTTATTTTTCT" .
"TATACTTAATTGATCTAACATAATAGTTTGCTCAAAATAATTATAGCAAC" .
"AAAGTATTATATTATGGATACTTTGTGTGTCTATACATGCTTATATA---" .
"----TGTTTA-------TGTATTAAGTAAAATGAATGAGTGATTAAAAGG" .
"ATG-AGAGGGAGGAATTAGGAATATTTTGTTATTATAAGGTACTTTCACC" .
"ATCCCTGAAGTGTTATAGTGTTATTTGAAA---GTAGACTTGGATTAATT" .
"GTAAACATATATTGCAAATTCTAGGGCAACCTCTAAAAAACATTTTTTAA" .
"AAAGTATATATTAATTAATATGCTA-AGAAAGGAGAGAAAATGGAATAAT" .
"ACAAAATGCTAAATTAAAGCCACAA-AAGGCAGAAAAAGAATAGAAGACA" .
"AAAATAGGAACAA-AGAATAATGAATAGAAAATAGTAA-TAACTATGGTA" .
"GATATTAATCCAAATGTATCAATAATAATTTTAAACATCAGTTGTCTAA-" .
"TATACCAATGAAAAGACAGA--TTGTCAGAGTAGATTAAAAAAAAT---C" .
"CAACTATATGTTGTATACAAGAAACCCACTTTTAATATAA",
subjString =>
"ATGACAGAGATGTTAGAATTATCTGACAAAGATTTTAAAGCAGCCATNAT" .
"AAAAATGCTTCAATGAGCAATTACGAACACGCTTGAAACAAATGAAA-AA" .
"TAGAAAGTCTCAGCAAAGAAATAGAAA-GTCTCAGCAAAGAAATAGAAGA" .
"TATAAAGAAGAACCA-AATGGAAATTTTAGAACTGAAAAATACAAT--AA" .
"CCGAAATAAAAAGCTCAGTGGATGGGCTCAACAGCAGAATGGAGGGGACA" .
"GAGGAAAGAATCAGTGAACTGGAAGATAGAACAATAGAAATTACCCAATC" .
"TGAACAACAGAGAGAAAATAGACTGAAAAAAAAAA-TGAACAGAGCCTCA" .
"GGGACCTGTGGGACTATAACAAAAGATCTAACATTTATGTCATCGGAGTC" .
"CCAGAAGGAGAGGAGAAAGAGGGCGGGGCTGAAAAAGTACTCGAAGAAAT" .
"AATGGCTGAAAACTTC-CCAAATT--TGGCAAGAGACATAAACCTACAGA" .
"TTCAAGAAGCTGAGCGAACCCCAAACAGGATAAAC-CCAAAGA-------" .
"-AATCC--ACACCAAGACACATCATAATCAAACTTCTGAAAACTAAAGAC" .
"AAAGAAAAAAAATCTTGAAAGCAGCNAGAGAGAAACGACACCTTACCTAT" .
"AGGGGAAAAACAATTCGAATGACAGCGGATTTCTCATCAGAAACCATGGA" .
"GGCCAGA---AGGA-AGTGGCACAACATTTTTCAAGTGCTGAAAGAAAA-" .
"GAACTGTCAACCCAGAATCCTATATCCAGCGAAAATATCCTTCAGGAATG" .
"AAGGGGAAATCAAGACATT-CTCAGATGAAGGAAAACTAAGAGAATTTGT" .
"CGCCAGCAGACCTACCCTAAAAGAATGGCTAAAGGAAGTTCTCT-AAACA" .
"GAAAGGAAATGATAAAAG-AAGGAATCTTGGAACATCAGGAAGGAAGAAA" .
"GAACA-CAGTAAGCAAAAATATGGGTAAATACAATAGGCTTTCCTTCTCC" .
"T---CTTGAGTTTTCTAA-ATTAT-GTTTG------ACGGTTGAAGCAA-" .
"AAATTATAACACTGT------CTGATGTGGTTCTAAATG--TATGTAGAG" .
"GAAATATTTAAGACAATTATATTA--TAAA-TGGGGGAGGG--TAAAGGG" .
"ACATAAAGGGAGG---TAAG--------GTTTCTATA------CTTCAC-" .
"-TC---AAACTGGTAAAATG--AC--GACACCAGTAGACTGTGATAAGTT" .
"ATGTATATATAATGTAATACCTAGAGCAACCACTAAAAAA-GCTATACAA" .
"AGAG-ATACACTCA--AAAACACTATAGATA--AATCAAAATGGAATTCT" .
"AAAAAATGTTCAAGTAAC-CCACAGGAAGGCAGGAAAA----AGAAAACA" .
"GAGAAATGAAAAACAGAACAA---ACAGAAAACAAAAAATAAA-ATGGCA" .
"GACTTAAGCCCTAACATATCAATAATTACATTAAATGTAAATGGTCTAAA" .
"TACACCAATTAAAAGACAGAGATTGGCAGAGTGGATTAAAAAACATGACC" .
"CAACTATATGCTGTCTACAAGAAACTCACTTCAAATATAA" );

#v*, **,  
#vi, v*,
#v*, i*,
#v*, vi,
#vi, ii,
#i*, i*,
#ii, i*,
#
#Normal: 
#trans = 11, tranv = 6 
#CpG:
#trans = 1/10 + 1/10 + 1/10 + 1/10 + 1 + 1/10 + 1/10 + 1 + 1/10 
#      = 2 - 7/10 = 2.7
#
#
#CpG transitions = normal - ( 8.3 )
#

## doesn't matter for this one
( $score, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions )  = $searchResult->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -30,
                             gapExtPenalty => -6,
                             complexityAdjust => 1,
                            );
                             #scoreCpGMod => 1

diag("Now score = $score\n");
diag("div = $divergence\n");
diag("transitions = $transitions\n");
diag("transversions = $transversions\n");


my $startTime = time();
for ( my $i = 0; $i < 1000; $i++ )
{
 ( $score, $divergence, $cpgsites, $percIns, $percDel,
     $positionScores, $xdrop_fragments, $well_characterized_bases,
     $transitions, $transversions )  = $searchResult->rescoreAlignment(
                             scoreMatrix => $matrix,
                             gapOpenPenalty => -30,
                             gapExtPenalty => -6,
                             scoreCpGMod => 1
                            );
}
diag("Time to run 1000 rescores: " . ( time() - $startTime ) . " secs ( typically 8 seconds )\n");
  


#diag( "" . $searchResult->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n" );

done_testing();

1;
