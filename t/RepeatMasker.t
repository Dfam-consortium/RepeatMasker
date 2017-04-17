print "1..5\n";

use strict;

my $testFile = "";
my $testDataDir = "";
my $resultComparisonDir = "";


#------------------------------------------------------------------------##
# Test to see if we broke something
#------------------------------------------------------------------------##
$testDataDir = "seqs/general/";
$resultComparisonDir = "seqs/general/small-1-rcmp-1/";
$testFile = "small-1.fa";
print "\n# Testing - Simple Run -gc 43 $testDataDir$testFile\n";
system("cp $testDataDir$testFile wrkdir/$testFile");
system("../RepeatMasker -engine crossmatch -gc 43 wrkdir/$testFile");
system("cat wrkdir/$testFile.tbl | egrep -v \"RepeatMasker|cross_match|RepBase\" > wrkdir/$testFile.tbl.cmp");
print "not " unless( 
  `diff -b wrkdir/$testFile.tbl.cmp $resultComparisonDir/$testFile.tbl.cmp` eq "" &&
  `diff wrkdir/$testFile.masked $resultComparisonDir/$testFile.masked` eq "" &&
  `diff wrkdir/$testFile.out $resultComparisonDir/$testFile.out` eq "" );
print "ok 1\n";

$testDataDir = "seqs/general/";
$resultComparisonDir = "seqs/general/small-1-rcmp-2/";
$testFile = "small-1.fa";
print "\n# Testing - SoftMasking RepeatMasker -gc 43 -xsmall $testDataDir$testFile\n";
system("cp $testDataDir$testFile wrkdir/$testFile");
system("../RepeatMasker -engine crossmatch -gc 43 -xsmall wrkdir/$testFile");
system("cat wrkdir/$testFile.tbl | egrep -v \"RepeatMasker|cross_match|RepBase\" > wrkdir/$testFile.tbl.cmp");
print 
  "diff -b wrkdir/$testFile.tbl.cmp $resultComparisonDir/$testFile.tbl.cmp\n". 
  "diff wrkdir/$testFile.masked $resultComparisonDir/$testFile.masked\n".
  "diff wrkdir/$testFile.out $resultComparisonDir/$testFile.out\n";
print "not " unless( 
  `diff -b wrkdir/$testFile.tbl.cmp $resultComparisonDir/$testFile.tbl.cmp` eq "" &&
  `diff wrkdir/$testFile.masked $resultComparisonDir/$testFile.masked` eq "" &&
  `diff wrkdir/$testFile.out $resultComparisonDir/$testFile.out` eq "" );
print "ok 2\n";

$testDataDir = "seqs/general/";
$resultComparisonDir = "seqs/general/small-1-rcmp-3/";
$testFile = "small-1.fa";
print "\n# Testing - Mask with Xs RepeatMasker -gc 43 -x $testDataDir$testFile\n";
system("cp $testDataDir$testFile wrkdir/$testFile");
system("../RepeatMasker -engine crossmatch -gc 43 -x wrkdir/$testFile");
system("cat wrkdir/$testFile.tbl | egrep -v \"RepeatMasker|cross_match|RepBase\" > wrkdir/$testFile.tbl.cmp");
print 
  "diff -b wrkdir/$testFile.tbl.cmp $resultComparisonDir/$testFile.tbl.cmp\n".
  "diff wrkdir/$testFile.masked $resultComparisonDir/$testFile.masked\n".
  "diff wrkdir/$testFile.out $resultComparisonDir/$testFile.out\n";
print "not " unless( 
  `diff -b wrkdir/$testFile.tbl.cmp $resultComparisonDir/$testFile.tbl.cmp` eq "" &&
  `diff wrkdir/$testFile.masked $resultComparisonDir/$testFile.masked` eq "" &&
  `diff wrkdir/$testFile.out $resultComparisonDir/$testFile.out` eq "" );
print "ok 3\n";

$testDataDir = "seqs/ISElements/";
$resultComparisonDir = "seqs/ISElements/is1-rcmp-1/";
$testFile = "is1.fa";
print "\n# Testing - ISClipping RepeatMasker -gc 43 -is_clip $testDataDir$testFile\n";
system("cp $testDataDir$testFile wrkdir/$testFile");
system("../RepeatMasker -engine crossmatch -gc 43 -is_clip wrkdir/$testFile");
system("cat wrkdir/$testFile.tbl | egrep -v \"RepeatMasker|cross_match|RepBase\" > wrkdir/$testFile.tbl.cmp");
print
  "diff -b wrkdir/$testFile.tbl.cmp $resultComparisonDir$testFile.tbl.cmp\n".
  "diff wrkdir/$testFile.masked $resultComparisonDir$testFile.masked\n".
  "diff wrkdir/$testFile.out $resultComparisonDir$testFile.out\n";
print "not " unless( 
  `diff -b wrkdir/$testFile.tbl.cmp $resultComparisonDir$testFile.tbl.cmp` eq "" &&
  `diff wrkdir/$testFile.masked $resultComparisonDir$testFile.masked` eq "" &&
  `diff wrkdir/$testFile.out $resultComparisonDir$testFile.out` eq "" );
print "ok 4\n";


$testDataDir = "seqs/general/";
$resultComparisonDir = "seqs/general/hum-1-rcmp-1/";
$testFile = "hum-1.fa";
print "\n# Testing - Med Human Sequence with RepeatMasker -gc 43 -qq $testDataDir$testFile\n";
system("cp $testDataDir$testFile wrkdir/$testFile");
system("../RepeatMasker -engine crossmatch -gc 43 -qq wrkdir/$testFile");
system("cat wrkdir/$testFile.tbl | egrep -v \"RepeatMasker|cross_match|RepBase\" > wrkdir/$testFile.tbl.cmp");
unless (
  `diff -b wrkdir/$testFile.tbl.cmp $resultComparisonDir/$testFile.tbl.cmp` eq "" &&
  `diff --ignore-blank-lines wrkdir/$testFile.masked $resultComparisonDir/$testFile.masked` eq "" &&
  `diff wrkdir/$testFile.out $resultComparisonDir/$testFile.out` eq "" ) {
  # Bad!
  print "diff wrkdir/$testFile.tbl.cmp $resultComparisonDir/$testFile.tbl.cmp\n";
  print "diff --ignore-blank-lines wrkdir/$testFile.masked $resultComparisonDir/$testFile.masked\n";
  print "diff wrkdir/$testFile.out $resultComparisonDir/$testFile.out\n";
  print "not ";
}
print "ok 5\n";


