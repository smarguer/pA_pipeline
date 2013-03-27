#!/usr/bin/perl

use strict;


use GD::Graph::linespoints;
use GD::Graph::lines;
use GD::Graph::bars;
use GD::Graph::pie;
use List::MoreUtils qw(uniq);
#use Math::CDF qw(:all);
use Math::Trig ':pi';


my @datasets=("0304_10");
#insert which datasets should be included: "UCLstrandspecific" (=proliferating cells, strand specific)
#					"fasta" (=proliferating cells, non strand specific)	
#					"meiotic" (=meiosis arrested cells, strand specific)	
#					"quiescent7days" (=quiescence arrested cell, 7 days of nitrogen depletion, strand specific)	
#					"quiescentpt2" (=quiescence arrested cell, 24 hours of nitrogen depletion, strand specific)	
#					"long" (=proliferatinf cells, strand specific, ~70nt reads)
my $difference=6; #difference between CS to be summarized into one (individual separation)
my $minmaps=0;	#minimal number of RNA-seq hits for identified CS to be considered
my $maxutrlength=1000; #maximal distance from stop codon for CS to be considered
my $insideoroutside="outside"; #for Motifs: CS inside or outside ORf (input "inside", "outside", "both")
my $cleavagemultiplicity=0; #0 stands for "all", for Motifs: UTRs with how many CS should be considered?
my $cleavageorder=0;	#0 stands for "all", for Motifs: which CS in UTR should be considered
my $tandemorconvergent="all"; #"all", "tand", or "convergent" which genes should be considered?
my $leftboundary=-50;	#for Motifs: which region should be scanned? left boundary of Motif position (regarding to start of motif, 0 is the first nt overlapping poly(A) tail)
my $rightboundary=-6;	#for Motifs: which region should be scanned? right boundary of Motif position (regarding to start of motif, 0 is the first nt overlapping poly(A) tail)

my $maxmotifnumber=500;
my $leftwindowsize=-50; #for Motif profile plotting, left boundary of window
my $rightwindowsize=10; #for Motif profile plotting, right boundary of window
my $numberofgraphs=10;	#how manyy top motifs should be plotted?
my $binsize=20;
my $minmotiflength=6;	#motifs of which length should be considered? (minimum)
my $maxmotiflength=6;	#motifs of which length should be considered? (maximum)
my $maxproportions=9;
my $maxgeneproportions=10;
my $basesbeforestop=1500;
my $bases=3000;
my $basesplotright=(-200); 	#for nt-distribution profile around CSs: left boundary of window
my $basesplotleft=200;		#for nt-distribution profile around CSs: right boundary of window
my @db;
my @db1;
#============================read in coords===========================================================
open(UTRsandGenesnewpmid,"<", "usedfiles/pombe3UTRsgenesnewpmid")or die("Unable to open 3UTRsandGenesnewpmid");
my $n=0;
my $m;
my @features;
my @u;
while(my $line=<UTRsandGenesnewpmid>)
{
	my @feat=(trim(substr($line,84,17)),substr($line,35,3),int(substr($line,53,10))-1,int(substr($line,62,14))-1);

	$u[$n]=int(substr($line,10,1))-1;
	if($n>0 && $u[$n] != $u[$n-1])
	{
	$m=$n;
	}
          push@{$features[$u[$n]][$n-$m]}, @feat;
	$n++;
}

# First index of @features: chromosome number
# Second index of @features: numbering of elements for each chromosome
# Third index of @features: each feature corresponding to the one of interes. In order: strand, region, start, stop
#============================read in coords for introns===============================================
open(Intronsnewpmid,"<", "usedfiles/pombeintronsnewpmid")or die("Unable to open Intronsnewpmid");
$n=0;
my @introns;

while(my $line=<Intronsnewpmid>)
{
	my @intr=(int(substr($line,10,1)),trim(substr($line,82,10)),int(substr($line,46,10))-1,int(substr($line,59,14))-1);
          push@{$introns[$n]}, @intr;
	$n++;
}

for(my $u=0; $u<=$#introns; $u++)
{
@{$db[$introns[$u][0]-1][$introns[$u][2]]}=("Intron",$u,$introns[$u][1],"start");
@{$db[$introns[$u][0]-1][$introns[$u][3]]}=("Intron",$u,$introns[$u][1],"stop");
push @{$db1[$introns[$u][0]-1]},[$introns[$u][2],"Intron",$u,$introns[$u][1],"start"];
push @{$db1[$introns[$u][0]-1]},[$introns[$u][3],"Intron",$u,$introns[$u][1],"stop"];
}
#======================================systematic ids===============================================
open(contig1,"<", "usedfiles/contig1")or die("Unable to open contig1");
open(contig2,"<", "usedfiles/contig2")or die("Unable to open contig2");
open(contig3,"<", "usedfiles/contig3")or die("Unable to open contig3");

my %contig;
for (my $u=1; $u<=3; $u++)
{
my $string="contig"."$u";
while(my $line=<$string>)
{
my @array=split(/ /,$line);
if($array[2]ne"c")
{
$contig{($u-1,$array[0]-1)}="contig=".'"'.$array[2].'";';

}
elsif($array[2]eq"c")
{
$contig{($u+2,$array[0]-1)}="contig=".'"'.$array[3].'";';

}
}
}
close(contig1);
close(contig2);
close(contig3);


#===================================Gene names======================================================
open(Geneschr1,"<", "usedfiles/geneschr1")or die("Unable to open Genesechr1");
open(Geneschr2,"<", "usedfiles/geneschr2")or die("Unable to open Genesechr2");
open(Geneschr3,"<", "usedfiles/geneschr3")or die("Unable to open Genesechr3");
my @genenames1;
my @genenames2;
my @genenames3;
my @genenames=(\@genenames1, \@genenames2, \@genenames3);
my @dbgenenames;
my %names;
my %dbnames;

my $h=0;
for (my $u=1; $u<=3; $u++)
{
my $r=-1;
my $string="Geneschr"."$u";
my $k=0;

while(my $line=<$string>)
{

	if($line=~/CDS/)
	{
	$h++;

		my $x=1;

		for(my $k=$r+1; $k<=$#{$features[($u-1)]}; $k++)
		{
			my $start=$features[($u-1)][$k][2]+1;
			if($x==1 && $line=~$start&&$features[($u-1)][$k][1] eq "CDS")
			{

			$x=0;
			$r=$k;
			@{$db[$u-1][$features[($u-1)][$k][2]]}=("CDS",$h,$features[($u-1)][$k][0],"start");
			@{$db[$u-1][$features[($u-1)][$k][3]]}=("CDS",$h,$features[($u-1)][$k][0],"stop");
			push @{$db1[$u-1]},[$features[($u-1)][$k][2],"CDS",$h,$features[($u-1)][$k][0],"start"];
			push @{$db1[$u-1]},[$features[($u-1)][$k][3],"CDS",$h,$features[($u-1)][$k][0],"stop"];
			$features[($u-1)][$k][5]=$h;
			}
		}
	if ($features[($u-1)][$r][0] eq "plus")
	{
	push(@{$names{(($u-1),$features[($u-1)][$r][2])}}, $contig{(($u-1),$features[($u-1)][$r][2])});
	push(@{$dbgenenames[($u-1)][$h]}, $contig{(($u-1),$features[($u-1)][$r][2])});
	}
	elsif($features[($u-1)][$r][0] eq "minus")
	{
	push(@{$names{(($u-1),$features[($u-1)][$r][2])}}, $contig{(($u+2),$features[($u-1)][$r][2])});
	push(@{$dbgenenames[($u-1)][$h]}, $contig{(($u+2),$features[($u-1)][$r][2])});
	}
	}
	elsif (($line=~/gene=/) || ($line=~/synonym/) || ($line=~/obsolete/))
	{
	$genenames[($u-1)][$r]= substr($line,0,28);
	$genenames[($u-1)][$r] =~ s/\R//g; #removes linebreak at end of the line
	push(@{$names{(($u-1),$features[($u-1)][$r][2])}}, $genenames[($u-1)][$r]);
	push(@{$dbgenenames[($u-1)][$h]}, $genenames[($u-1)][$r]);
	}

}
}
close(Geneschr1);
close(Geneschr2);
close(Geneschr3);

my @dbnames;
for(my $u=0; $u<=2; $u++)
{
	for(my $v=0; $v<=$#{$dbgenenames[$u]}; $v++)
	{
	my @geneq;
	my @synq;
	my @contigq;
		for(my $w=0; $w<=$#{$dbgenenames[$u][$v]}; $w++)
		{
		my @array1=split(/"/,$dbgenenames[$u][$v][$w]);
			if ($array1[0]eq"contig=")
			{
			push @contigq, $array1[1];
			}
			elsif ($array1[0]eq"gene=")
			{
			push @geneq, $array1[1];
			}
			else
			{
			push @synq, $array1[1];
			}
		}

		my @sorted = sort { $a cmp $b } @geneq;	

	my @array1=(@contigq,@sorted, @synq);
		for(my $w=0; $w<=$#{$dbgenenames[$u][$v]}; $w++)
		{
	 $dbnames[$u][$v][$w]=$array1[$w];

		}
	}
}



#======================================================================================================
open(Chromosome1,"<", "usedfiles/wholeseqchr1") or die("Unable to open chromosome1");
open(Chromosome2,"<", "usedfiles/wholeseqchr2") or die("Unable to open chromosome2");
open(Chromosome3,"<", "usedfiles/wholeseqchr3") or die("Unable to open chromosome3");
open(Chromosome4,"<", "usedfiles/wholechr1rev") or die("Unable to open chromosome1rev");
open(Chromosome5,"<", "usedfiles/wholechr2rev") or die("Unable to open chromosome2rev");
open(Chromosome6,"<", "usedfiles/wholechr3rev") or die("Unable to open chromosome3rev");

my @chromosome1;
my @chromosome2;
my @chromosome3;
my @chromosome1rev;
my @chromosome2rev;
my @chromosome3rev;
my @chromosomes=(\@chromosome1,\@chromosome2,\@chromosome3,\@chromosome1rev,\@chromosome2rev,\@chromosome3rev);

for (my $u=1; $u<=6; $u++)
{
my $string="Chromosome"."$u";
	while(my $line=<$string>)
	{
	chomp($line);
	my @intermediate = split(//,$line);


	push @{$chromosomes[$u-1]},@intermediate;

	}
}

close(Chromosome1);
close(Chromosome2);
close(Chromosome3);
close(Chromosome4);
close(Chromosome5);
close(Chromosome6);



print "Chromosomes read\n";

#======================================================================================================

my $dataset;
my @n;
my @loadcs1fwd;
my @loadcs2fwd;
my @loadcs3fwd;
my @loadcs1rev;
my @loadcs2rev;
my @loadcs3rev;
my @loadcs1opp;
my @loadcs2opp;
my @loadcs3opp;

my @loadcs=(\@loadcs1fwd,\@loadcs2fwd,\@loadcs3fwd,\@loadcs1rev, \@loadcs2rev, \@loadcs3rev);




my @loadocc1fwd;
my @loadocc2fwd;
my @loadocc3fwd;
my @loadocc1rev;
my @loadocc2rev;
my @loadocc3rev;
my @loadocc1opp;
my @loadocc2opp;
my @loadocc3opp;
my @loadocc=(\@loadocc1fwd,\@loadocc2fwd,\@loadocc3fwd,\@loadocc1rev, \@loadocc2rev, \@loadocc3rev);
for (my $u=1; $u<=6; $u++)
{
$n[$u]=0;
}
foreach (@datasets)
{
$dataset=$_;
print $dataset."\n";
open(csocc1,"<", 'usedfiles/'.$dataset.'/'.$dataset.'cs1fwdocc') or die("Unable to open csocc1");
open(csocc2,"<", 'usedfiles/'.$dataset.'/'.$dataset.'cs2fwdocc') or die("Unable to open csocc2");
open(csocc3,"<", 'usedfiles/'.$dataset.'/'.$dataset.'cs3fwdocc') or die("Unable to open csocc3");
open(csocc4,"<", 'usedfiles/'.$dataset.'/'.$dataset.'cs1revocc') or die("Unable to open csocc4");
open(csocc5,"<", 'usedfiles/'.$dataset.'/'.$dataset.'cs2revocc') or die("Unable to open csocc5");
open(csocc6,"<", 'usedfiles/'.$dataset.'/'.$dataset.'cs3revocc') or die("Unable to open csocc6");



my @loadocc=(\@loadocc1fwd,\@loadocc2fwd,\@loadocc3fwd,\@loadocc1rev, \@loadocc2rev, \@loadocc3rev);

for (my $u=1; $u<=6; $u++)
{
my $string="csocc"."$u";
	while(my $line=<$string>)
	{
my @ary=split("\t",$line);
	$loadcs[$u-1][$n[$u]]=int($ary[0]);
	$loadocc[$u-1][$n[$u]]=int($ary[1]);
	$n[$u]++;
	}
}

}
for (my $u=0; $u<=5; $u++)
{
my @array;
	for (my $v=0; $v<=$#{$loadcs[$u]};$v++)
	{
		$array[$v][0]=$loadcs[$u][$v];
		$array[$v][1]=$loadocc[$u][$v];
		
	}
my @array1=sort{$a->[0] <=> $b->[0]} @array;
	for (my $v=0; $v<=$#{$loadcs[$u]};$v++)
	{
	$loadcs[$u][$v]=$array1[$v][0];
	$loadocc[$u][$v]=$array1[$v][1];
	}
}

for (my $u=0; $u<=5; $u++)
{
	for (my $v=0; $v<=$#{$loadcs[$u]};$v++)
	{
#	$loadcs[$u][$v]-=1;
	}
}

my @loadcsfwd=(\@loadcs1fwd, \@loadcs2fwd, \@loadcs3fwd);
my @loadcsrev=(\@loadcs1rev, \@loadcs2rev, \@loadcs3rev);
my @loadcsopp=(\@loadcs1opp, \@loadcs2opp, \@loadcs3opp);
my @loadoccfwd=(\@loadocc1fwd, \@loadocc2fwd, \@loadocc3fwd);
my @loadoccrev=(\@loadocc1rev, \@loadocc2rev, \@loadocc3rev);
my @loadoccopp=(\@loadocc1opp, \@loadocc2opp, \@loadocc3opp);

for (my $u=0; $u<=2; $u++)
{
	for(my $v=0; $v<=$#{$loadcsrev[$u]}; $v++)
	{
	$loadcsopp[$u][$v]=$#{$chromosomes[$u]}-$loadcsrev[$u][$v]-1;
	$loadoccopp[$u][$v]=$loadoccrev[$u][$v];
	}

@{$loadcsopp[$u]}=reverse @{$loadcsopp[$u]};
@{$loadoccopp[$u]}=reverse @{$loadoccopp[$u]};
}


close(csocc1);
close(csocc2);
close(csocc3);
close(csocc4);
close(csocc5);
close(csocc6);

print $#{$loadcsopp[0]}." ".$#{$loadcsopp[1]}." ".$#{$loadcsopp[2]}." ".$#{$loadcsfwd[0]}." ".$#{$loadcsfwd[1]}." ".$#{$loadcsfwd[2]}."\n";
print "cleavage sites read\n";

print "cleavage sites read\n";


#==================Maximal distance between cleavage====================================================

my @csmed;

my @cs1fwd;
my @cs2fwd;
my @cs3fwd;
my @cs1rev;
my @cs2rev;
my @cs3rev;

my @cs1opp;
my @cs2opp;
my @cs3opp;

my @cs=(\@cs1fwd,\@cs2fwd,\@cs3fwd,\@cs1rev, \@cs2rev, \@cs3rev);


my @occmed;

my @occ1fwd;
my @occ2fwd;
my @occ3fwd;
my @occ1rev;
my @occ2rev;
my @occ3rev;

my @occ1opp;
my @occ2opp;
my @occ3opp;

my @occ=(\@occ1fwd,\@occ2fwd,\@occ3fwd,\@occ1rev, \@occ2rev, \@occ3rev);
my @occfwd=(\@occ1fwd, \@occ2fwd, \@occ3fwd);
my @occrev=(\@occ1rev, \@occ2rev, \@occ3rev);
my @occopp=(\@occ1opp, \@occ2opp, \@occ3opp);

for (my $u=0; $u<=5; $u++)
{
my $r=0;
my $i=0;
	while( $i<=$#{$loadcs[$u]})
	{
	my $k=1;
		if ($loadcs[$u][$i+1]-$loadcs[$u][$i] <=$difference)
		{	
			while ($loadcs[$u][$i+$k]-$loadcs[$u][$i+$k-1] <=$difference && $i+$k<=$#{$loadcs[$u]})
			{
			$k++;
			}
			if(maxloc(@{$loadocc[$u]}[$i..($i+$k-1)])!=0)
			{
			$csmed[$u][$r]=$loadcs[$u][$i+maxloc(@{$loadocc[$u]}[$i..($i+$k-1)])];
			$occmed[$u][$r]=$loadocc[$u][$i+maxloc(@{$loadocc[$u]}[$i..($i+$k-1)])];
			}
			else
			{	
			if($u<=2)
			{
			$csmed[$u][$r]=$loadcs[$u][$i]; #select shortest one
			$occmed[$u][$r]=$loadocc[$u][$i];
			}
			else
			{
			$csmed[$u][$r]=$loadcs[$u][$i+$k-1]; #select shortest one
			$occmed[$u][$r]=$loadocc[$u][$i+$k-1];
			}
			}
		
		$i=$i+$k-1;

		}
		else 
		{
		$csmed[$u][$r]=$loadcs[$u][$i];

		$occmed[$u][$r]=$loadocc[$u][$i];	

		}	
	$i++;		
	$r++;
	}

	my $s=0;
	for (my $t=0; $t<=$#{$csmed[$u]}; $t++)
	{
	if($occmed[$u][$t]>=$minmaps)
		{
		 $cs[$u][$s]=$csmed[$u][$t];

		$occ[$u][$s]=$occmed[$u][$t];
		$s++;
		
		}
	}
}


my @csfwd=(\@cs1fwd, \@cs2fwd, \@cs3fwd);
my @csrev=(\@cs1rev, \@cs2rev, \@cs3rev);
my @csopp=(\@cs1opp, \@cs2opp, \@cs3opp);

for (my $u=0; $u<=2; $u++)
{
	for(my $v=0; $v<=$#{$csrev[$u]}; $v++)
	{
	$csopp[$u][$v]=$#{$chromosomes[$u]}-$csrev[$u][$v]-1;
	$occopp[$u][$v]=$occrev[$u][$v];
	}

@{$csopp[$u]}=reverse @{$csopp[$u]};
@{$occopp[$u]}=reverse @{$occopp[$u]};
}

for (my $u=3; $u<=5;$u++)
{
for(my $v=0; $v<=$#{$cs[$u]}; $v++)
{
for (my $f=0; $f<=10;$f++)
{
#rint $chromosomes[$u][$csopp[$u-3][$v]-10+$f];
}
#print "\n";
}
}

print "Maximal distance between cleavage done\n";

#===========================================Orientation for each cleavage site==========================================================
my @direction;
my @plusfeatures;
my @minusfeatures;
my @cds;
my $g=0;
for (my $u=0; $u<=2; $u++)
{
	for (my $v=0; $v<=$#{$features[$u]}; $v++)
	{
		if ($features[$u][$v][1] eq "CDS")
		{
		@{$cds[$g]}=($u, $v);
		push(@{$cds[$g]}, @{$features[$u][$v]});		
		$g++;
		}
	}
}



for (my $u=0; $u<=$#cds; $u++)
{
if ($u<$#cds && $cds[$u][0] eq $cds[$u+1][0] && $cds[$u][2] eq "plus" && $cds[$u+1][2] eq "plus")
{
$direction[$cds[$u][0]][$cds[$u][1]]="tand";
$features[$cds[$u][0]][$cds[$u][1]][4]=$direction[$cds[$u][0]][$cds[$u][1]];
}
elsif($u<$#cds && $cds[$u][0]==$cds[$u+1][0] && $cds[$u][2] eq "plus" && $cds[$u+1][2] eq "minus")
{
$direction[$cds[$u][0]][$cds[$u][1]]="cvg";
$direction[$cds[$u+1][0]][$cds[$u+1][1]]="cvg";
$features[$cds[$u][0]][$cds[$u][1]][4]=$direction[$cds[$u][0]][$cds[$u][1]];
$features[$cds[$u+1][0]][$cds[$u+1][1]][4]=$direction[$cds[$u+1][0]][$cds[$u+1][1]];
}
elsif($u>0 && $cds[$u][0]==$cds[$u-1][0] && $cds[$u][2] eq "minus" && $cds[$u-1][2] eq "minus")
{
$direction[$cds[$u][0]][$cds[$u][1]]="tand";
$features[$cds[$u][0]][$cds[$u][1]][4]=$direction[$cds[$u][0]][$cds[$u][1]];
}
elsif($features[$cds[$u][0]][$cds[$u][1]][4] ne "tand" && $features[$cds[$u][0]][$cds[$u][1]][4] ne "cvg")
{
$direction[$cds[$u][0]][$cds[$u][1]]="none";
$features[$cds[$u][0]][$cds[$u][1]][4]=$direction[$cds[$u][0]][$cds[$u][1]];
}
}


for (my $u=0; $u<=2; $u++)
{
	for (my $v=0; $v<=$#{$features[$u]}; $v++)
	{

		if ($features[$u][$v][1] eq "UTR")
		{
			for(my $h=0; $h<=10; $h++)
			{
			if ($v>$h && $features[$u][$v][0] eq "plus"&& $features[$u][($v-$h-1)][1] eq "CDS"  && $features[$u][$v][2] == ($features[$u][($v-$h-1)][3]+1))
			{
			$direction[$u][$v]=$features[$u][$v-$h-1][4];
			$features[$u][$v][4]=$direction[$u][$v];
			
			}
			elsif ($v<($#{$features[$u]}-$h) && $features[$u][$v][0] eq "minus"&& $features[$u][($v+$h+1)][1] eq "CDS"  && ($features[$u][$v][3]+1) ==$features[$u][$v+1+$h][2])
			{
			$direction[$u][$v]=$features[$u][$v+1+$h][4];
			$features[$u][$v][4]=$direction[$u][$v];
			}
			
			}
			if($features[$u][$v][4] ne "tand" &&$features[$u][$v][4] ne "cvg"&&$features[$u][$v][4] ne "none")
			{
			$features[$u][$v][4]="randomUTR";
			}
		}
	}
}
print "Tnd and Cvg determined\n";
#=====================================================================================================

my @cdsfeatures;
my %cleavagefeatures; 				#chromosome, genenumber, tand/cvg, occurence, inside or outside, UTRlength,cleavagemultiplicity,cleavageorder
my %oppcleavagefeatures; 				#chromosome, genenumber, tand/cvg, occurence, inside or outside, UTRlength
my %cdsnames;
my @cleavages;
my @utrlengths;
my @independentgene;
my @outside;
my @inside;
my @cleavagenames;
my %genestop;


for (my $u=0; $u<=2; $u++)
{	
	$n=0;
	$m=0;
	for (my $v=0; $v<=$#{$features[$u]}; $v++)
	{
		if($features[$u][$v][1] eq "CDS" && $features[$u][$v][0] eq "plus")
		{
		$cdsfeatures[$u][$n]=$features[$u][$v];
		$cdsnames{(($u),$features[$u][$v][2])}=$names{(($u),$features[$u][$v][2])};
		$n++;
		}
		elsif($features[$u][$v][1] eq "CDS" && $features[$u][$v][0] eq "minus")
		{
		$cdsfeatures[$u+3][$m]=$features[$u][$v];
		$cdsnames{(($u+3),$features[$u][$v][2])}=$names{(($u),$features[$u][$v][2])};	
		$m++;
		}
	}
}


my $totalcds=0;
my $tandems=0;
my $convergents=0;
my $cleavagenumber;

for (my $u=0; $u<=5; $u++)
{
	 for(my $v=0; $v<=$#{$cdsfeatures[$u]}; $v++)
	{
	$totalcds++;
		if($cdsfeatures[$u][$v][4] eq "tand")
		{
		$tandems++;
		}
		elsif($cdsfeatures[$u][$v][4] eq "cvg")
		{
		$convergents++;
		}
	}
}
print "total number of genes is: ".$totalcds."\n";
print "out of which tandem: ".$tandems." and convergent: ".$convergents."\n";

my $cleavagenumber;

for (my $u=0; $u<=2; $u++)
{
	 for(my $v=0; $v<$#{$cdsfeatures[$u]}; $v++)
	{
		my $k=0;
		while($k<=$#{$csfwd[$u]})
		{
			if($csfwd[$u][$k]>=$cdsfeatures[$u][$v][3] && $csfwd[$u][$k]<$cdsfeatures[$u][$v+1][2] && ($csfwd[$u][$k]-$cdsfeatures[$u][$v][3])<=$maxutrlength)
			{
			my $str="$u"."_"."$v";
				push @independentgene, $str;
				push @outside, $str;
				$cleavagenumber=0;

				while($csfwd[$u][$k]>=$cdsfeatures[$u][$v][3] && $csfwd[$u][$k]<$cdsfeatures[$u][$v+1][2] && ($csfwd[$u][$k]-$cdsfeatures[$u][$v][3])<=$maxutrlength)
				{

				$cleavagenumber++;
				push(@{$cleavages[$u][$v]}, $csfwd[$u][$k]);
				push(@{$utrlengths[$u][$v]}, ($csfwd[$u][$k]-$cdsfeatures[$u][$v][3]));
				$k++;
				}
				for( my $n=$k-$cleavagenumber; $n<$k; $n++)
				{
				@{$cleavagefeatures{($u,$csfwd[$u][$n])}}=($u, $v, $cdsfeatures[$u][$v][4], $occfwd[$u][$n],"outside",($csfwd[$u][$n]-$cdsfeatures[$u][$v][3]),$cleavagenumber,($cleavagenumber-($k-$n)+1),$n);
				push @{$db1[$u]},[$csfwd[$u][$n],"CS",$cdsfeatures[$u][$v][5],"plus", $cleavagenumber,$n,"outside"];
				}

			}
			if($csfwd[$u][$k]>=$cdsfeatures[$u][$v][2] && $csfwd[$u][$k]<$cdsfeatures[$u][$v][3])
			{
				my $str="$u"."_"."$v";
				push @independentgene, $str;
				push @inside, $str;		
			while($csfwd[$u][$k]>=$cdsfeatures[$u][$v][2] && $csfwd[$u][$k]<$cdsfeatures[$u][$v][3])
				{
				@{$cleavagefeatures{($u,$csfwd[$u][$k])}}=($u, $v, $cdsfeatures[$u][$v][4], $occfwd[$u][$k],"inside","na","na","na",$k);
				push @{$db1[$u]},[$csfwd[$u][$k],"CS",$cdsfeatures[$u][$v][5],"plus","na",$k,"inside"];
				$k++;
				}
			}
			else
			{
			$k++;
			}
			
		}
	}
}


for (my $u=3; $u<=5; $u++)
{
	 for(my $v=0; $v<$#{$cdsfeatures[$u]}; $v++)
	{
		my $k=0;
		while($k<=$#{$csrev[$u-3]})
		{
			if($csrev[$u-3][$k]>$cdsfeatures[$u][$v][3] && $csrev[$u-3][$k]<=$cdsfeatures[$u][$v+1][2] && ($cdsfeatures[$u][$v+1][2]-$csrev[$u-3][$k])<=$maxutrlength)
			{
				my $str="$u"."_"."$v";
				push @independentgene, $str;
				push @outside, $str;
				$cleavagenumber=0;
				while($csrev[$u-3][$k]>$cdsfeatures[$u][$v][3] && $csrev[$u-3][$k]<=$cdsfeatures[$u][$v+1][2] && ($cdsfeatures[$u][$v+1][2]-$csrev[$u-3][$k])<=$maxutrlength)
				{
				my $str="$u"."_"."$v";
				$cleavagenumber++;
				push(@{$cleavages[$u][($v+1)]}, $csrev[$u-3][$k]);
				push(@{$utrlengths[$u][($v+1)]}, ($cdsfeatures[$u][$v+1][2]-$csrev[$u-3][$k]));
				$k++;
				}	
				for( my $n=$k-$cleavagenumber; $n<$k; $n++)
				{
				@{$cleavagefeatures{($u,$csrev[$u-3][$n])}}=($u, ($v+1), $cdsfeatures[$u][($v+1)][4], $occrev[$u-3][($n)],"outside",($cdsfeatures[$u][$v+1][2]-$csrev[$u-3][$n]),$cleavagenumber,($k-$n),$n);
				push @{$db1[$u-3]},[$csrev[$u-3][$n],"CS",$cdsfeatures[$u][$v+1][5],"minus", $cleavagenumber,$n,"outside"];
				@{$oppcleavagefeatures{($u,($#{$chromosomes[$u-3]}-$csrev[$u-3][$n]-1))}}=($u, ($v+1), $cdsfeatures[$u][($v+1)][4], $occrev[$u-3][($n)],"outside",($cdsfeatures[$u][$v+1][2]-$csrev[$u-3][$n]), $cleavagenumber, ($k-$n),$n);


				}
			}
			if($csrev[$u-3][$k]>$cdsfeatures[$u][$v][2] && $csrev[$u-3][$k]<=$cdsfeatures[$u][$v][3])
			{
				my $str="$u"."_"."$v";
				push @independentgene, $str;
				push @inside, $str;
				while($csrev[$u-3][$k]>$cdsfeatures[$u][$v][2] && $csrev[$u-3][$k]<=$cdsfeatures[$u][$v][3])
				{
				@{$cleavagefeatures{($u,$csrev[$u-3][$k])}}=($u, $v, $cdsfeatures[$u][$v][4], $occrev[$u-3][$k],"inside","na","na","na",$k);
				push @{$db1[$u-3]},[$csrev[$u-3][$k],"CS",$cdsfeatures[$u][$v][5],"minus", "na",$k,"inside"];

				@{$oppcleavagefeatures{($u,($#{$chromosomes[$u-3]}-$csrev[$u-3][$k]-1))}}=($u, $v, $cdsfeatures[$u][$v][4], $occrev[$u-3][$k],"inside","na","na","na",$k);


				$k++;
				}
			}
			else
			{
			$k++;
			}
			
		}
	}
}





print "genes and cleavage sites written\n";


my @unindependentgene=uniq @independentgene;
my @unoutside=uniq @outside;
my @uninside=uniq @inside;

my $numberofindependentgenes=($#unindependentgene+1);
my $numberofoutsidegenes=($#unoutside+1);
my $numberofinsidegenes=($#uninside+1);

#=================================Tandem and convergent depending on number of cleavage sites=====================================================================
open(GENESPLITTING,">", 'Results/'.$dataset.'/GENESPLITTING.dat') or die("Unable to open GENESPLITTING");
my @tandem;
my @convergent;
my @tandemcleavages;
my @convergentcleavages;
my @cleavagecounter;

for (my $u=0; $u<=5; $u++)
{
	 for(my $v=0; $v<=$#{$cs[$u]}; $v++)
	{
		if ($cleavagefeatures{($u,$cs[$u][$v])}[4] eq "outside")
		{
		$cleavagecounter[($cleavagefeatures{($u,$cs[$u][$v])}[6]-1)]++;
		
			if ($cleavagefeatures{($u,$cs[$u][$v])}[2] eq "tand")
			{
			$tandem[($cleavagefeatures{($u,$cs[$u][$v])}[6]-1)]++;
			$tandemcleavages[$cleavagefeatures{($u,$cs[$u][$v])}[0]][$cleavagefeatures{($u,$cs[$u][$v])}[1]]=$cleavages[$u][$v];
			}
			if ($cleavagefeatures{($u,$cs[$u][$v])}[2] eq "cvg")
			{
			$convergent[($cleavagefeatures{($u,$cs[$u][$v])}[6]-1)]++;
			$convergentcleavages[$cleavagefeatures{($u,$cs[$u][$v])}[0]][$cleavagefeatures{($u,$cs[$u][$v])}[1]]=$cleavages[$u][$v];
			}
		}
	}
}

	 for(my $v=0; $v<=$#tandem; $v++)
	{
	$tandem[$v]/=($v+1);
	print GENESPLITTING "tandem ".($v+1)." ".$tandem[$v]."\n";
	}
	print GENESPLITTING "tandem "."all"." ".sum(@tandem)."\n";

	 for(my $v=0; $v<=$#convergent; $v++)
	{
	$convergent[$v]/=($v+1);
	print GENESPLITTING "convergent ".($v+1)." ".$convergent[$v]."\n";
	}
	print GENESPLITTING "convergent "."all"." ".sum(@convergent)."\n";

	 for(my $v=0; $v<=$#cleavagecounter; $v++)
	{
	$cleavagecounter[$v]/=($v+1);
	print GENESPLITTING "all ".($v+1)." ".$cleavagecounter[$v]."\n";
	}
	print "all"." ".sum(@cleavagecounter)."\n";
	

	print GENESPLITTING "All genes considered:".$numberofindependentgenes." genes with cleavages outside:".$numberofoutsidegenes." genes with cleavages inside:".$numberofinsidegenes."\n";
close(GENESPLITTING);
print "Genesplitting done\n";


#=================================Select kind of cleavage sites to be considered=======================================================
my @cleavage;
my @oppcleavage;

for(my $u=0; $u<=5; $u++)
{
my $n=0;
	for(my $v=0; $v<=$#{$cs[$u]}; $v++)
	{
	if ($cleavagefeatures{($u,$cs[$u][$v])}[5]>=0 || $cleavagefeatures{($u,$cs[$u][$v])}[5] eq "na")
	{
		if($insideoroutside ne 0 && $cleavagefeatures{($u,$cs[$u][$v])}[4] eq $insideoroutside )
		{
			if($cleavagemultiplicity ne 0 && $cleavagefeatures{($u,$cs[$u][$v])}[6] eq $cleavagemultiplicity)
			{

				if($cleavageorder ne 0 && $cleavagefeatures{($u,$cs[$u][$v])}[7] eq $cleavageorder)
				{

					if($tandemorconvergent ne "all" && $cleavagefeatures{($u,$cs[$u][$v])}[2] eq $tandemorconvergent)
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}
					elsif($tandemorconvergent eq "all")
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}

				}
				elsif($cleavageorder == 0)
				{
					if($tandemorconvergent ne "all" && $cleavagefeatures{($u,$cs[$u][$v])}[2] eq $tandemorconvergent)
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}
					elsif($tandemorconvergent == "all")
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}

				}

			}
			elsif($cleavagemultiplicity == 0)
			{
				if($cleavageorder ne 0 && $cleavagefeatures{($u,$cs[$u][$v])}[7] eq $cleavageorder)
				{
					if($tandemorconvergent ne "all" && $cleavagefeatures{($u,$cs[$u][$v])}[2] eq $tandemorconvergent)
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}
					elsif($tandemorconvergent eq "all")
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}

				}
				elsif($cleavageorder == 0)
				{
					if($tandemorconvergent ne "all" && $cleavagefeatures{($u,$cs[$u][$v])}[2] eq $tandemorconvergent)
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}
					elsif($tandemorconvergent eq "all")
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}

				}

			}


		}
		elsif($insideoroutside eq 0)
		{
			if($cleavagemultiplicity ne 0 && $cleavagefeatures{($u,$cs[$u][$v])}[6] eq $cleavagemultiplicity)
			{
				if($cleavageorder ne 0 && $cleavagefeatures{($u,$cs[$u][$v])}[7] eq $cleavageorder)
				{
					if($tandemorconvergent ne "all" && $cleavagefeatures{($u,$cs[$u][$v])}[2] eq $tandemorconvergent)
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}
					elsif($tandemorconvergent == "all")
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}

				}
				elsif($cleavageorder == 0)
				{
					if($tandemorconvergent ne "all" && $cleavagefeatures{($u,$cs[$u][$v])}[2] eq $tandemorconvergent)
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}
					elsif($tandemorconvergent eq "all")
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}

				}

			}
			elsif($cleavagemultiplicity == 0)
			{
				if($cleavageorder ne 0 && $cleavagefeatures{($u,$cs[$u][$v])}[7] eq $cleavageorder)
				{
					if($tandemorconvergent ne "all" && $cleavagefeatures{($u,$cs[$u][$v])}[2] eq $tandemorconvergent)
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}
					elsif($tandemorconvergent eq "all")
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}

				}
				elsif($cleavageorder == 0)
				{
					if($tandemorconvergent ne "all" && $cleavagefeatures{($u,$cs[$u][$v])}[2] eq $tandemorconvergent)
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}
					elsif($tandemorconvergent eq "all")
					{
					$cleavage[$u][$n]=$cs[$u][$v];
					$n++;
					}

				}

			}

		}
	}	
	}
}

for(my $u=3; $u<=5; $u++)
{
	for(my $v=0; $v<=$#{$cleavage[$u]}; $v++)
	{
	$oppcleavage[$u-3][$v]=$#{$chromosomes[$u-3]}-$cleavage[$u][$v]-1;
#print $cleavage[$u][$v]." ".$#{$chromosomes[$u-3]}." ".$u." ".$v."\n";
	}

@{$oppcleavage[$u-3]}=reverse@{$oppcleavage[$u-3]};
}

print "cleavage type selected\n";

#=================================Sequences from chromosome============================================================================

my @utrcleavage=(\@{$cleavage[0]},\@{$cleavage[1]},\@{$cleavage[2]},\@{$oppcleavage[0]},\@{$oppcleavage[1]},\@{$oppcleavage[2]});




my $totalcleavage;
for (my $u=0; $u<=$#utrcleavage; $u++)
{
$totalcleavage+=($#{$utrcleavage[$u]}+1);
}
print "number of sequences is ".$totalcleavage."\n";

my $newindependentgene=0;

for(my $u=0; $u<=2; $u++)
{
	for(my $v=1; $v<=$#{$utrcleavage[$u]}; $v++)
	{
		if($cleavagefeatures{($u,$utrcleavage[$u][$v])}[1]ne$cleavagefeatures{($u,$utrcleavage[$u][$v-1])}[1])
		{
		$newindependentgene++;
		}
	}
}

for(my $u=3; $u<=5; $u++)
{
	for(my $v=1; $v<=$#{$utrcleavage[$u]}; $v++)
	{
		if($oppcleavagefeatures{($u,$utrcleavage[$u][$v])}[1]ne$oppcleavagefeatures{($u,$utrcleavage[$u][$v-1])}[1])
		{
		$newindependentgene++;
		}
	}
}

$newindependentgene+=6;


print "New number of independent genes:".$newindependentgene."\n";
#=================================basecounter==========================================================================================
#open(BASECOUNTER,">", 'Results/'.$dataset.'/basecounter.dat') or die("Unable to open BASECOUNTER");
my @zeros;
my @zeros1;
my @zeros2;
my @zeros3;
		
for(my $f=0; $f<=$bases; $f++)
{
$zeros[$f]=0;
$zeros1[$f]=0;
$zeros2[$f]=0;
$zeros3[$f]=0;
}
my %basecounter=("A"=>[@zeros],"T"=>[@zeros1],"G"=>[@zeros2],"C"=>[@zeros3]);



for(my $u=0; $u<=5; $u++)
{
	for(my $v=0; $v<=$#{$utrcleavage[$u]}; $v++)
	{
		for(my $f=0; $f<=$bases; $f++)
		{
			if($basesbeforestop<=$utrcleavage[$u][$v] && ($utrcleavage[$u][$v]-$basesbeforestop+$bases)<=$#{$chromosomes[$u]})
			{
			$basecounter{$chromosomes[$u][($utrcleavage[$u][$v]-$basesbeforestop+$f)]}[$f]++;
			
			} 
		}
	}
} 

my $key;
my $value;
while(($key,$value)=each(%basecounter))
{
for(my $f=0; $f<=$bases; $f++)
{
#print BASECOUNTER $basecounter{$key}[$f]." ";
}
#print BASECOUNTER "\n";
}


my @xs;
my @As;
my @Gs;
my @Cs;
my @Ts;
for(my $u=$basesplotright; $u<=$basesplotleft; $u++)
{
$xs[$u-$basesplotright]="$u";
}
for(my $u=($bases-$basesbeforestop+$basesplotright); $u<=($bases-$basesbeforestop+$basesplotleft); $u++)
{
$As[$u-($bases-$basesbeforestop+$basesplotright)]=$basecounter{A}[$u];
$Gs[$u-($bases-$basesbeforestop+$basesplotright)]=$basecounter{G}[$u];
$Cs[$u-($bases-$basesbeforestop+$basesplotright)]=$basecounter{C}[$u];
$Ts[$u-($bases-$basesbeforestop+$basesplotright)]=$basecounter{T}[$u];
}

my @data=(\@xs,\@As,\@Gs,\@Cs,\@Ts);
my $maxa=maxval(@As);
my $maxc=maxval(@Cs);
my $maxg=maxval(@Gs);
my $maxt=maxval(@Ts);
my @bb=($maxa, $maxc, $maxg, $maxt);
my $maxb=maxval(@bb);

  my $graph = GD::Graph::lines->new(1500, 1000);

  $graph->set( 
      x_label           => 'Distance from cleavage site',
      y_label           => 'bases at x-position out of'.$totalcleavage.'sequences',
      title             => 'Base distribution around cleavage site. location: '.$insideoroutside.'; Multiplicity: '.$cleavagemultiplicity.'; order:'.$cleavageorder.'; tnd or cvg'.$tandemorconvergent,
      y_max_value       => ($maxb+5),
      y_tick_number     => 15,
      y_label_skip      => 5,
      x_label_skip      => 5,
      transparent	=> 0
  ) or die $graph->error;

$graph->set_legend(("A","G","C","T"));
$graph->set( dclrs => [ qw(lred lgreen lblue cyan lpurple lorange) ] );

  my $gd = $graph->plot(\@data) or die $graph->error;

open(Basespic, '>Results/'.$dataset.'/basecounter.png') or die $!;
 binmode Basespic;
 print Basespic $gd->png;




#close(BASECOUNTER);
print"basecounter done\n";

#============================Calculate frequencies of nucleotides etc.=============================================
my @unicspergene;
my %frequency;
my $markovorder=1;
my @genestop; #indeces correspond to @unicspergene
for(my $u=0; $u<=2; $u++)
{
my $r=0;
my $v=0;
while($v<=$#{$utrcleavage[$u]})
{
my $k=1;
	if ($cleavagefeatures{($u, $utrcleavage[$u][$v])}[1] == $cleavagefeatures{($u, $utrcleavage[$u][$v+$k])}[1])
	{
		while($cleavagefeatures{($u, $utrcleavage[$u][$v])}[1] == $cleavagefeatures{($u, $utrcleavage[$u][$v+$k])}[1])
		{
		$k++;
		}
	$unicspergene[$u][$r]=$utrcleavage[$u][$v+$k-1];
	$genestop[$u][$r]=$cdsfeatures[$u][$cleavagefeatures{($u,$utrcleavage[$u][$v+$k-1])}[1]][3];
	$r++;
	}
	else
	{
	$unicspergene[$u][$r]=$utrcleavage[$u][$v];
	$genestop[$u][$r]=$cdsfeatures[$u][$cleavagefeatures{($u,$utrcleavage[$u][$v])}[1]][3];
	$r++
	}
$v+=$k;
}
}


for(my $u=3; $u<=5; $u++)
{
my $r=0;
my $v=0;
while($v<=$#{$utrcleavage[$u]})
{
my $k=1;
	if ($oppcleavagefeatures{($u, $utrcleavage[$u][$v])}[1] == $oppcleavagefeatures{($u, $utrcleavage[$u][$v+$k])}[1])
	{
		while($oppcleavagefeatures{($u, $utrcleavage[$u][$v])}[1] == $oppcleavagefeatures{($u, $utrcleavage[$u][$v+$k])}[1])
		{my @genestop;
		$k++;
		}
	$unicspergene[$u][$r]=$utrcleavage[$u][$v+$k-1];
	$genestop[$u][$r]=$#{$chromosomes[$u-3]}-$cdsfeatures[$u][$oppcleavagefeatures{($u,$utrcleavage[$u][$v+$k-1])}[1]][2]-1;
	$r++;
	}
	else
	{
	$unicspergene[$u][$r]=$utrcleavage[$u][$v];
	$genestop[$u][$r]=$#{$chromosomes[$u-3]}-$cdsfeatures[$u][$oppcleavagefeatures{($u,$utrcleavage[$u][$v])}[1]][2]-1;
	$r++
	}

$v+=$k;
}
}


my $totalnumberofnt;

for(my $g=1; $g<=($markovorder+1); $g++)
{
for(my $u=0; $u<=5; $u++)
{
	for(my $v=0; $v<=$#{$unicspergene[$u]}; $v++)
	{
	my $c=$genestop[$u][$v];
	my $d=$unicspergene[$u][$v];
	$totalnumberofnt+=(abs($c-$d)+1);
	my @array=@{$chromosomes[$u]}[$c..$d];
	my $intstring="@array";
	$intstring =~ s/(.)\s/$1/seg; #removes white spaces
		for(my $i=0; $i<=($#array-$g+1); $i++)
		{	
		my $intkey=substr($intstring,$i,$g);
		$frequency{$intkey}++;
		}
	}
}
}

my $key;
foreach $key(keys %frequency)
{
$frequency{$key}=$frequency{$key}/($totalnumberofnt-length($key)+1);
}


#================================Motifproportion===================================================================





my @maxmotifs;

for (my $maxlength=$minmotiflength; $maxlength<=$maxmotiflength; $maxlength++)
{
open(PVAL,">", 'db/'.$dataset.'/pvals_motiflength'.$maxlength."_leftbd".$leftboundary."_righbd".$rightboundary."_".$cleavagemultiplicity.'_order'.$cleavageorder.'_tndcvg'."$tandemorconvergent"."$insideoroutside".'.dat') or die("Unable to open PVAL");
print "Total number of sequences considered is ".$totalcleavage."\n";
print $maxlength."\n";
my %propcleavage;
my @seqfrequency;
for(my $u=0; $u<=5; $u++)
{
	for (my $v=0; $v<=$#{$utrcleavage[$u]}; $v++)
	{
	@{$propcleavage{($u,$v)}}=($u,$v,$utrcleavage[$u][$v]);
	}
}


	my $r=0;
	my $bound=$totalcleavage;

#while($r<=0)#$maxproportions && $bound!=0)
while($bound>0)		
	{
	my %motifprop;
	my %sequences;
	my %genes;
		my $key;
		foreach $key(keys %propcleavage)
			{
			my $cut=1;
			my $c=($propcleavage{$key}[2]+$leftboundary);
			my $d=($propcleavage{$key}[2]+$rightboundary);
			my @array=@{$chromosomes[$propcleavage{$key}[0]]}[$c..$d];
			my $intstring="@array";
			$intstring =~ s/(.)\s/$1/seg; #removes white spaces


				for(my $i=0; $i<=($#array-$maxlength+1); $i++)
				{	
					my $intkey=substr($intstring,$i,$maxlength);
					$motifprop{$intkey}[0][0]++;
					$motifprop{$intkey}[1][($motifprop{$intkey}[0][0]-1)][0]=$propcleavage{$key}[0];
					$motifprop{$intkey}[1][($motifprop{$intkey}[0][0]-1)][1]=$propcleavage{$key}[1];
					$motifprop{$intkey}[1][($motifprop{$intkey}[0][0]-1)][2]=$i;
					$motifprop{$intkey}[1][($motifprop{$intkey}[0][0]-1)][2]=$c+$i;
					push @{$sequences{$intkey}},"$propcleavage{$key}[0]"."_"."$propcleavage{$key}[1]";
					if ($propcleavage{$key}[0]<=2)
					{
					push @{$genes{$intkey}},"$cleavagefeatures{($propcleavage{$key}[0],$propcleavage{$key}[2])}[0]"."_"."$cleavagefeatures{($propcleavage{$key}[0],$propcleavage{$key}[2])}[1]";
					}
					else
					{
					push @{$genes{$intkey}},"$oppcleavagefeatures{($propcleavage{$key}[0],$propcleavage{$key}[2])}[0]"."_"."$oppcleavagefeatures{($propcleavage{$key}[0],$propcleavage{$key}[2])}[1]";
					}
				}

			
		}

my $key2;
my %uni;
my %rank;
my $gh=0;

		foreach $key2(keys %sequences)
		{
#$gh++;
#print $key2." ".$gh." \n";
		my @un = @{$sequences{$key2}};
		my @unique= uniq @un;
		$uni{$key2}=($#unique+1);


		}


my $key2;
my @maxes;
foreach $key2 (sort {$uni{$a} <=> $uni{$b} }
           keys %uni)
{
push (@maxes, $key2);
}
@maxes=reverse @maxes;

my $size;
if($#maxes>100)
{
$size=100;
}
else
{
$size=$#maxes;
}

for (my $kj=0; $kj<=$size; $kj++)
{
print $kj."\n";
my @motifpos;
for(my $g=0; $g<=$bases; $g++)
{
$motifpos[$g]=0;
}
		my @positions = ();
		my $KEY;
		foreach $KEY (keys %propcleavage)
		{	
			my $U=$propcleavage{$KEY}[0];
			my $CLS=$propcleavage{$KEY}[2];
			my $tempstring=join("",@{$chromosomes[$U]}[($CLS-$basesbeforestop)..($CLS-$basesbeforestop+$bases)]);
			my $pos = 0;



			while (1) 
			{
			  $pos = index($tempstring, $maxes[$kj], $pos);
			  last if($pos < 0);
			  push(@positions, $pos++);
			}
		}

			for (my $k=0; $k<=$#positions; $k++)
			{
			$motifpos[$positions[$k]]++;				
			}


		$rank{$maxes[$kj]}=(maxval(@motifpos[(int($bases/2)+1+$leftboundary)..(int($bases/2)+1+$rightboundary)])-med(@motifpos));

		}
my @unisort;

my $key4;
foreach $key4 (sort {$rank{$a} <=> $rank{$b} }
           keys %uni)
{
push (@unisort, $key4);
}
@unisort=reverse @unisort;



my %pvalue;
my $e=0;
my $N=$bound;
print "N ".$N."\n";
		for (my $i=0; $i<=9; $i++)
		{
		my $key3=$unisort[$i];
		print $key3."\n";
			
			my $K=$uni{$key3};
			my $h=$K-1;

			if ($K != 0)
			{	

			my $a=substr($key3,0,2);
			my $numerator=$frequency{$a};
			my $denominator=1;
			for (my $t=1; $t<=(length($key3)-2); $t++)
			{
			my $e=substr($key3,$t,2);
			my $f=substr($key3,$t,1);
			$numerator*=$frequency{$e};
			$denominator*=$frequency{$f};
			}
			for(my $t=$K;$t<=$N; $t++)
			{
			my $P;
 			my $probability=($numerator/$denominator);
			my $qprobability=(1-($numerator/$denominator));
			my $n=abs($leftboundary-$rightboundary)+1;
			$P=(1-($qprobability**($n)));
			$pvalue{$key3}+=exp(logsum($N)-logsum($t)-logsum($N-$t)+logzum($t,$P)+logzum(($N-$t),(1-$P)));
			}
			}

		}
	

my @orderedprops;	


my $key1;
foreach $key1 (sort {$pvalue{$a} <=> $pvalue{$b} }
           keys %pvalue)
{
push (@orderedprops, $key1);
}

my $stddev;
my $avg;
my $sig=0;
my $x=1;


my %cl;
my %unig;
			
		my @ung = @{$genes{$orderedprops[$sig]}};
		my @uniqueg= uniq @ung;
		$unig{$orderedprops[$sig]}=($#uniqueg+1);


		for(my $jj=0; $jj<=$uni{$orderedprops[$sig]}-1; $jj++)
		{
		my @un = @{$sequences{$orderedprops[$sig]}};
		my @unique= uniq @un;
		@{$cl{$orderedprops[$sig]}[$jj]}=split(/_/, $unique[$jj]);
		$cl{$orderedprops[$sig]}[$jj][0]=int($cl{$orderedprops[$sig]}[$jj][0]);
		$cl{$orderedprops[$sig]}[$jj][1]=int($cl{$orderedprops[$sig]}[$jj][1]);
		}

my @positions = ();
		for(my $k=0; $k<=$#{$cl{$orderedprops[$sig]}}; $k++)
		{
		my $tempstring=join("",@{$chromosomes[$cl{$orderedprops[$sig]}[$k][0]]}[($propcleavage{($cl{$orderedprops[$sig]}[$k][0],$cl{$orderedprops[$sig]}[$k][1])}[2]+$leftboundary)..($propcleavage{($cl{$orderedprops[$sig]}[$k][0],$cl{$orderedprops[$sig]}[$k][1])}[2]+$rightboundary)]);
			my $pos = 0;



			while (1) 
			{
			  $pos = index($tempstring, $orderedprops[$sig], $pos);
			  last if($pos < 0);
			  push(@positions, $pos++);
			}
		}
		
	
	$stddev=stdev(\@positions);
	$avg=average(\@positions);


my @motifpositions=();
my @positions = ();
		for(my $k=0; $k<=$#{$cl{$orderedprops[0]}}; $k++)
		{
		my $tempstring=join("",@{$chromosomes[$cl{$orderedprops[0]}[$k][0]]}[($propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2]+$leftboundary)..($propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2]+$rightboundary)]);
			my $pos = 0;



			while (1) 
			{
			  $pos = index($tempstring, $orderedprops[0], $pos);
			  last if($pos < 0);
			my $dum=(($propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2]+$leftboundary)+$pos);
			if($cl{$orderedprops[0]}[$k][0]<=2)
			{
			my $t1=$cleavagefeatures{($cl{$orderedprops[0]}[$k][0], $propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2])}[0];
			my $t2=$cleavagefeatures{($cl{$orderedprops[0]}[$k][0], $propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2])}[1];
			my $t3=$cleavagefeatures{($cl{$orderedprops[0]}[$k][0], $propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2])}[6];
			my $t4=$cleavagefeatures{($cl{$orderedprops[0]}[$k][0], $propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2])}[8];
			my $r1=$r+1;
			push @{$db1[$cl{$orderedprops[0]}[$k][0]]},[$dum,$orderedprops[0],$cdsfeatures[$t1][$t2][5],"plus",$t3,$t4,"outside",$r1];
			}	
			elsif($cl{$orderedprops[0]}[$k][0]>2)
			{
			my $t1=$oppcleavagefeatures{($cl{$orderedprops[0]}[$k][0], $propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2])}[0];
			my $t2=$oppcleavagefeatures{($cl{$orderedprops[0]}[$k][0], $propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2])}[1];
			my $t3=$oppcleavagefeatures{($cl{$orderedprops[0]}[$k][0], $propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2])}[6];
			my $t4=$oppcleavagefeatures{($cl{$orderedprops[0]}[$k][0], $propcleavage{($cl{$orderedprops[0]}[$k][0],$cl{$orderedprops[0]}[$k][1])}[2])}[8];
			my $cr=$#{$chromosomes[$cl{$orderedprops[0]}[$k][0]]}-$dum-1;
			my $r1=$r+1;
			push @{$db1[$cl{$orderedprops[0]}[$k][0]-3]},[$cr,$orderedprops[0],$cdsfeatures[$t1][$t2][5],"minus",$t3,$t4,"outside",$r1];
			}	
			  push(@positions, $pos++);

			}
		}
		
	
my $fufi=($leftboundary+$avg);
		print PVAL $orderedprops[0]." ".$pvalue{$orderedprops[0]}." ".$uni{$orderedprops[0]}." ".$fufi." ".$stddev." ".$unig{$orderedprops[0]}."\n";
		print $orderedprops[0]." ".$pvalue{$orderedprops[0]}." ".$uni{$orderedprops[0]}." ".$fufi." ".$stddev." ".$unig{$orderedprops[0]}."\n";

push @{$maxmotifs[$maxlength]}, $orderedprops[0];
push @seqfrequency, $uni{$orderedprops[0]}."\n";

		for(my $k=0; $k<=($motifprop{$orderedprops[0]}[0][0]-1); $k++)
		{
		delete($propcleavage{($motifprop{$orderedprops[0]}[1][$k][0],$motifprop{$orderedprops[0]}[1][$k][1])});
		}
$bound=($totalcleavage-sum(@seqfrequency));
my $frum=keys ( %propcleavage );
print "bound ".$bound." ".$frum."\n";
$r++;
	}
close(PVAL);

}


print "Motifproportion done\n";
for(my$u=0; $u<=2;$u++)
{
my @db2=sort {$a->[0] <=> $b->[0];} @{$db1[$u]};
my $oi=($u+1);
open(DB,">", 'db/'.$dataset.'/db'.$oi.'_'.$dataset) or die("Unable to open db".$u);
	for(my$v=0; $v<=$#db2;$v++)
	{
		for(my $w=0; $w<=$#{$db2[$v]};$w++)
		{
		print DB $db2[$v][$w]."\t";
		}
		if($db2[$v][1]eq"CDS")
		{
		for (my $er=0; $er<=$#{$dbnames[$u][$db2[$v][2]]}; $er++)
		{
		print DB  $dbnames[$u][$db2[$v][2]][$er]."\t";
		}
		}
	print DB "\n";
	}
close(DB);
}

#========================================SUBROUTINES===================================================
#subroutine to remove whitespaces from strings=========================================================
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


#subroutine to get the maximum of an array=======================================
sub maxval
{

my @array=@_;
my $max_value = 0;

for(my $i = 0; $i<=$#array;$i++)
{
  if($array[$i] > $max_value)
	{
	$max_value = $array[$i];
	}
}
return $max_value;
}



#subroutine to get the maximum location of an array=======================================
sub maxloc
{
my @array=@_;
my $max_index = 0;
my $max_value = 0;
for(my $i = 0; $i <=$#array;$i++)
{
  if($array[$i] > $max_value)
	{
	$max_value = $array[$i];
	$max_index = $i;
	}
}
return $max_index;
}

#subroutine to get the minimum of an array=======================================
sub minval
{
my @array=@_;
my $min_value = $#chromosome1;

for(my $i = 0; $i <=$#array;$i++)
{
  if($array[$i] < $min_value)
	{
	$min_value = $array[$i];
	}
}
return $min_value;
}


#subroutine to get the minimum location of an array=======================================
sub minloc
{
my @array=@_;
my $min_index = 0;guenter grass
my $min_value = $#chromosome1;
for(my $i = 0; $i <=$#array;$i++)
{
  if($array[$i] < $min_value)
	{
	$min_value = $array[$i];
	$min_index = $i;
	}
}
return $min_index;
}

#subroutine to sum array-elements==========================================================
sub sum
{
my @array=@_;
my $sumarray;
for(my $row = 0; $row <=$#array; $row++)
	{
	$sumarray+=$array[$row];
	}
return $sumarray;
}
#subroutine to get complimentary strand in direction of transcription=====================
sub invert
{
my @array=@_;
my @invarray;
	for(my $a=0; $a<=$#array; $a++)
	{
		if($array[$a] eq "A")
		{
		@invarray[$a]="T";
		}
		elsif($array[$a] eq "C")
		{
		@invarray[$a]="G";
		}
		elsif($array[$a] eq "G")
		{
		@invarray[$a]="C";
		}
		elsif($array[$a] eq "T")
		{
		@invarray[$a]="A";
		}
		else
		{
		@invarray[$a]=@array[$a];
		}
	}
@invarray=reverse(@invarray);
return @invarray;
}
#==========================================================

sub med
{

my @array= @_;
my $count = $#array;
my @arraysort = sort { $a <=> $b } (@array);
	if ($count % 2)
	{
	return $arraysort[int($count/2)];
	} 
	else 
	{
	return ($arraysort[$count/2] + $arraysort[$count/2 - 1]) / 2;
	}
}
#=====================MEan===============================
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

#============standard deviation==============================
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

#==========================================================
sub fac {
  my ($n) = $_;

  if ($n < 2) {
    return 1;
  }
  else {
    return $n * fac($n-1);
  }
}
#==================================
  sub log10 {
 my $n = shift;
 return log($n)/log(10);
 }
#==========================================================
sub logsum {
 # my ($n) = $_;

  if ($_[0] == 1) {
    return log($_[0]);
  }
  elsif  ($_ [0]> 1){
    return log($_[0]) + logsum($_[0]-1);
  }
}
#==========================================================
sub logzum {
  my @n = @_;

    return $n[0]*(log($n[1]));
}
