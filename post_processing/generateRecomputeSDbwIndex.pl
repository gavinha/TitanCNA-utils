#!/usr/bin/perl -w 

=head1 NAME

generateClusterIndexEval.pl

-head1 SYNOPSIS

=head1 OPTIONS

    -inDir|i <string>
    -outDir|o <string>
    -datatype|d <string>  {'LogRatio','AllelicRatio'}; default: LogRatio
	-scale <integer> Default: 25
    -rscript|r <string>  default: ~/software/code/scripts/titan/analysis/recomputeSDbwIndex.R
    -script|s <string>

=head1 DESCRIPTION

=head1 CONTACT

Gavin Ha <gha@bccrc.ca>

=cut

use strict;
use DBI;
use File::Basename;
use Getopt::Long;

sub usage () {
    exec('perldoc', $0);
    exit;
}

my ($inDir,$outDir,$datatype,$scale,$rscript,$script,$Rbinary,$help);
$Rbinary = "/xchip/cga_home/gavinha/software/R-3.1.1/bin/R";
$datatype = "Both";
$scale = 25;
$rscript = "~/software/code/scripts/titan/analysis/recomputeSDbwIndex.R";
GetOptions (
    'inDir|i=s' => \$inDir,
    'outDir|o=s' => \$outDir,
    'datatype|d=s' => \$datatype,
    'scale=i' => \$scale,
    'script|s=s' => \$script,
    'rscript|r=s' => \$rscript,
    'help|?' => \$help
    );

if($help) {
    &usage();
}
print "Parameters:\ninDir=$inDir\ndatatype=$datatype\nscale=$scale\noutDir=$outDir\nrscript=$rscript\nscript=$script\n";

my $ls = `ls -1 $inDir/*titan.txt`;
my @files = split(/\n/,$ls);

open OUT, ">$script" || die("Can't write to $script\n");
foreach my $file (@files){
	my ($name,$path,$suffix) = fileparse($file);
	my $id = $name;
	$id =~ s/_titan.txt//g;
	my $objFile = $file; $objFile =~ s/_titan.txt/.RData/g;
	my $outfile = "$outDir/$id\_params.txt";
	print OUT "$Rbinary --vanilla --args $file $objFile $scale $datatype $outfile < $rscript\n";

}
close OUT;
