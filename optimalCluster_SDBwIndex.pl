#!/usr/bin/perl -w 

=head1 NAME

optimalCluster_SDBwIndex.pl

-head1 SYNOPSIS

=head1 OPTIONS

    -inDir|i <string>
    -scale|s <float>
    -outFile|o <string>

=head1 DESCRIPTION

=head1 CONTACT

Gavin Ha <gha@bccrc.ca>

=cut

use strict;
use DBI;
use File::Basename;
use Getopt::Long;
use List::Util qw(min);

sub usage () {
    exec('perldoc', $0);
    exit;
}

my ($inDir,$outFile,$scale,$help);
GetOptions (
    'inDir|i=s' => \$inDir,
    'outFile|o=s' => \$outFile,
    'scale|s=f' => \$scale,
    'help|?' => \$help
    );

if($help) {
    &usage();
}
print "Parameters:\ninDir=$inDir\noutFile=$outFile\nscale=$scale\n";

my $ls = `ls -1 $inDir/*params.txt`;
#my $grep = `grep "S_Dbw validity index" $inDir/*/titan/*params.txt`;
#my @indexArray = split(/\n/,$grep);
my @files = split(/\n/,$ls);

my %sdbw;
foreach my $file (@files){
	my ($name,$path,$suffix) = fileparse($file); my $id = $name;
	$id =~ s/\_params.txt//g;
	my @token = split(/\_/,$id); 
	my $clust = $token[-1];
	my $sampleid = $id; $sampleid =~ s/\_$clust//g;
	$clust =~ s/cluster0//g;
	#print $name . "\t" . $id . "\t" . $sampleid . "\t" . $clust . "\n";
	my ($value,$norm,$ploidy,$cellPrev,$dens,$scat);
	open IN, $file || die("Can't open $file\n");
	while (<IN>){
		chomp;
		my ($type,$val) = split(/\t/,$_);
		
		if ($type =~ m/S_Dbw validity index/){
			$value = $val;
			#print $sampleid . "\t" . $clust . "\t" . $value . "\n";
		}elsif ($type =~ m/S_Dbw dens.bw/){
			$dens = $val;
		}elsif ($type =~ m/S_Dbw scat/){
			$scat = $val;
		}elsif ($type =~ m/Normal contamination estimate/){
			$norm = $val;
		}elsif ($type =~ m/Average tumour ploidy estimate/){
			$ploidy = $val;
		}elsif ($type =~ m/Clonal cluster cellular prevalence/){
			$val =~ s/\ /\,/g;
			$cellPrev = $val;
		}
	}
	#print "$scale\t$dens\t$scat\n";
	my $validity = $value;
	$validity = $scale*$dens + $scat if (defined $scale);
	$sdbw{$sampleid}{$validity}{'cluster'} = $clust;
	$sdbw{$sampleid}{$validity}{'normal'} = $norm;
	$sdbw{$sampleid}{$validity}{'ploidy'} = $ploidy;
	$sdbw{$sampleid}{$validity}{'cellPrev'} = $cellPrev;
	close IN;

}

open OUT, ">$outFile" || die("Can't write to $outFile\n");
foreach my $sid (sort keys %sdbw){
	my @k = keys %{$sdbw{$sid}};
	my @sdbwStr;
	foreach my $v (@k){
		push(@sdbwStr, $sdbw{$sid}{$v}{'cluster'}."=".$v);
	}
	#print join(", ",sort(@sdbwStr)) . "\n";
	my $minVal = min @k;
	#print $minVal . "\n";
	print OUT $sid . "_cluster0" . $sdbw{$sid}{$minVal}{'cluster'} . "\t" . $sid . "\t" . $sdbw{$sid}{$minVal}{'cluster'} . "\t" . $sdbw{$sid}{$minVal}{'cellPrev'} . "\t" . $sdbw{$sid}{$minVal}{'normal'} . "\t" . $sdbw{$sid}{$minVal}{'ploidy'} . "\t" . join(", ",sort(@sdbwStr)) . "\n";
}
close OUT;



