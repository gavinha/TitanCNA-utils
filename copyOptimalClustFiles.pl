#!/usr/bin/perl -w 

=head1 NAME

.pl

-head1 SYNOPSIS

=head1 OPTIONS

	-inDir|i <string>
    -inFile|f <string>     #optimal clusters file: uses first column
    -clustNum|c <integer>  #copy only TITAN results for cluster [$clustNum]
    -keyFile|k <string>    #2 column file: 1) ID to use, 2) File ID
    -mvSeg|seg <boolean>   Use if want to move .seg and segs.txt files
    -mvParam|param <boolean>
    -mvTitan|titan <boolean>
    -mvPlot|plot <boolean>
    -outDir|o <string>

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

my ($inDir,$inFile,$keyFile,$clustNum,$mvSeg,$mvParam,$mvTitan,$mvPlot,$outDir,$help);
GetOptions (
	'inDir|i=s' => \$inDir,
    'inFile|f=s' => \$inFile,
    'keyFile|k=s' => \$keyFile,
    'clustNum|c=i' => \$clustNum,
    'mvSeg|seg' => \$mvSeg,
    'mvParam|param' => \$mvParam,
    'mvTitan|titan' => \$mvTitan,
    'mvPlot|plot' => \$mvPlot,
    'outDir|o=s' => \$outDir,
    'help|?' => \$help
    );

if($help) {
    &usage();
}

#$clustNum = undef if (!defined($clustNum));
#$mvSeg = 0 if (!defined($mvSeg));
#$mvParam = 0 if (!defined($mvParam));
#$mvTitan = 0 if (!defined($mvTitan));

print "Parameters:\ninDir=$inDir\ninFile=$inFile\nkeyFile=$keyFile\nclustNum=$clustNum\noutDir=$outDir\nmvSeg=$mvSeg\nmvParam=$mvParam\nmvTitan=$mvTitan\nmvPlot=$mvPlot\n";

my %idmap;
if (defined $keyFile){
open KEY, $keyFile || die("Can't open $keyFile.\n");
while (<KEY>){
	chomp;
	my ($idtoUse,$fileid) = split(/\t/,$_);
	$idmap{$fileid} = $idtoUse;
}
close KEY;
}

my $ls = `ls -1 $inDir/*params.txt`;
$ls = 
my @files = split(/\n/,$ls);

open IN, $inFile || die("Can't open $inFile.\n");
while(<IN>){
	chomp;
	my ($fileid,$norm,$ploidy) = split(/\t/,$_);
	my $filepath = "";
	foreach my $file (@files){
		my ($name,$path,$suffix) = fileparse($file);
		$name =~ s/_params.txt//g;
		if ($name =~ m/$fileid/){
			my ($tmpId, @jnk) = split(/\_/,$fileid); 
			
			my $idToUse = $fileid; 
			if (defined $keyFile){
				my $newId = $idmap{$tmpId};
				$idToUse =~ s/$tmpId/$newId/g;
			}
			$filepath = "$path/$fileid";
			
			## use only specified cluster if given ####
			$filepath =~ s/cluster\_$jnk[-1]/cluster\_$clustNum/g if (defined $clustNum);
			$idToUse =~ s/cluster\_$jnk[-1]/cluster\_$clustNum/g if (defined $clustNum);  
			#print "filepath=$filepath\tfileid=$fileid\ttmpid=$tmpId\tnewId=$newId\tidToUse=$idToUse\n";
			#### params.txt ####
			my $tmpFile = "$filepath\_params.txt";
			my $outFile = "$outDir/$idToUse\_params.txt";
			if (-e $tmpFile && defined $mvParam){
				#print "cp $tmpFile $outFile\n";
				system("cp $tmpFile $outFile");
			}
			#### segs.txt ####
			$tmpFile = "$filepath\_segs.txt";
			$outFile = "$outDir/$idToUse\_segs.txt";
			if (-e $tmpFile && defined $mvSeg){
				#print "cp $tmpFile $outFile\n";
				system("cp $tmpFile $outFile");
			}
			#### .segs ####
			$tmpFile = "$filepath\.seg";
			$outFile = "$outDir/$idToUse\.seg";
			if (-e $tmpFile && defined $mvSeg){
				#print "cp $tmpFile $outFile\n";
				system("cp $tmpFile $outFile");
			}
			#### titan.txt ####
			$tmpFile = "$filepath\_titan.txt";
			$outFile = "$outDir/$idToUse\_titan.txt";
			if (-e $tmpFile && defined $mvTitan){
				#print "cp $tmpFile $outFile\n";
				system("cp $tmpFile $outFile");
			}
			#### plots ####
			$tmpFile = "$filepath\/";
			$outFile = "$outDir/$idToUse\/";
			if (-e $tmpFile && defined $mvPlot){
				#print "cp $tmpFile $outFile\n";
				system("cp -r $tmpFile $outFile");
			}
		}
	}
	print "Can't find file for $fileid\n" if ($filepath eq "");
	
}
close IN;

