#!/usr/local/bin/perl
use strict;
use strict;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/../";
use Data::Dumper;
use lib "$FindBin::RealBin/../Libraries";
use RepeatAnnotationData;


my %classHash = ();
open IN,"../famdb.py -i ../Libraries/RepeatMaskerLib.h5 families --descendants 1 --curated --format fasta_name --include-class-in-name |" or die;
while ( <IN> ) {
  if ( />(\S+)\#(\S+)/ ) {
    $classHash{lc($1)} = $2;
  }
}
close IN;


my $repeatDB = \%RepeatAnnotationData::repeatDB;

foreach my $id ( keys %classHash ) {
  if ( exists $repeatDB->{$id} ) {
    my $type = $classHash{$id};
    my $subtype = "";
    if ( $classHash{$id} =~ /(\S+)\/(\S+)/ ) {
      $type = $1;
      $subtype = $2;
    }
    if ( $repeatDB->{$id}->{'type'} ne "" && $repeatDB->{$id}->{'type'} ne $type ) {
       print "$id  - Type Change: Dfam=$type and RepeatAnnotationData=" . $repeatDB->{$id}->{'type'} . "\n";
       $repeatDB->{$id}->{'type'} = $type;
    }
    if ( $repeatDB->{$id}->{'subtype'} ne "" && $repeatDB->{$id}->{'subtype'} ne $subtype ) {
       print "$id  - Subtype Change: Dfam=$subtype and RepeatAnnotationData=" . $repeatDB->{$id}->{'subtype'} . "\n";
       $repeatDB->{$id}->{'subtype'} = $subtype;
    }
  }
}
##
## Output results
##
my %newDB;
open OUT, ">RepeatAnnotationDataModified.pm";
print OUT "package RepeatAnnotationData;\n";
print OUT "require Exporter;\n";
print OUT "\@EXPORT_OK = qw( \%repeatDB \%lineHash \%preProcData );\n";
print OUT "\%EXPORT_TAGS = ( all => [ \@EXPORT_OK ] );\n";
print OUT "\@ISA         = qw(Exporter);\n";
print OUT "\n";
print OUT "BEGIN {\n";
print OUT "\n";
print OUT "  \%repeatDB = (\n";
my $output = Dumper( \%newDB );
$output =~ s/^\$VAR1 = \{//;
$output =~ s/\};//;
print OUT $output;
print OUT "\n";
print OUT "  );\n";
print OUT "  \%lineHash = (\n";
$output = Dumper( \%RepeatAnnotationData::lineHash );
$output =~ s/^\$VAR1 = \{//;
$output =~ s/\};//;
print OUT $output;
print OUT "\n";
print OUT "  );\n";
print OUT "  \%preProcData = (\n";
$output = Dumper( \%RepeatAnnotationData::preProcData );
$output =~ s/^\$VAR1 = \{//;
$output =~ s/\};//;
print OUT $output;
print OUT "\n";
print OUT "  );\n";
print OUT "}\n";
close OUT;

## Done
exit;


