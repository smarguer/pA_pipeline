#!/usr/bin/perl

use strict;

use List::MoreUtils qw(uniq);
use Math::Trig ':pi';

my $dataset="0304_10";


#============================read in coords===========================================================
open(UTRsandGenesnewpmid,"<", "../usedfiles/pombe3UTRsgenesnewpmid")or die("Unable to open 3UTRsandGenesnewpmid");
my $n=0;
my $m;
my @features;
my @plusfeatures;
my @minusfeatures;
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

for(my $u=0; $u<=2; $u++)
{
	for(my $t=0; $t<=$#features;$t++)
	{
		if($features[$u][$t][0]eq"plus"&&$features[$u][$t][1]eq "CDS")
		{
		push @{$plusfeatures[$u]}, \@{$features[$u][$t]};
		}
		elsif($features[$u][$t][0]eq"minus"&&$features[$u][$t][1]eq "CDS")
		{
		push @{$minusfeatures[$u]}, \@{$features[$u][$t]};
		}
	}
}

# First index of @features: chromosome number
# Second index of @features: numbering of elements for each chromosome
# Third index of @features: each feature corresponding to the one of interes. In order: strand, region, start, stop
#======================================================================================================
open(Chromosome1,"<", "../usedfiles/wholeseqchr1") or die("Unable to open chromosome1");
open(Chromosome2,"<", "../usedfiles/wholeseqchr2") or die("Unable to open chromosome2");
open(Chromosome3,"<", "../usedfiles/wholeseqchr3") or die("Unable to open chromosome3");
open(Chromosome4,"<", "../usedfiles/wholechr1rev") or die("Unable to open chromosome1rev");
open(Chromosome5,"<", "../usedfiles/wholechr2rev") or die("Unable to open chromosome2rev");
open(Chromosome6,"<", "../usedfiles/wholechr3rev") or die("Unable to open chromosome3rev");

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
open(polyAs1filtered,"<", '../forfilter3/'.$dataset.'1fwdforfilter3') or die("Unable to open polyAs1filtered");
open(polyAs2filtered,"<", '../forfilter3/'.$dataset.'2fwdforfilter3') or die("Unable to open polyAs2filtered");
open(polyAs3filtered,"<", '../forfilter3/'.$dataset.'3fwdforfilter3') or die("Unable to open polyAs3filtered");
open(polyAs4filtered,"<", '../forfilter3/'.$dataset.'1revforfilter3') or die("Unable to open polyAs1revfiltered");
open(polyAs5filtered,"<", '../forfilter3/'.$dataset.'2revforfilter3') or die("Unable to open polyAs2revfiltered");
open(polyAs6filtered,"<", '../forfilter3/'.$dataset.'3revforfilter3') or die("Unable to open polyAs3revfiltered");

my $div=0;
my @rnaseq;
my @seqno;
my @chromseq;
my @map;
my @cs;
my @csrev;
my @quality;
for(my $u=1; $u<=6; $u++)
{
my $string="polyAs".$u."filtered";
	while(my $line=<$string>)
	{
	$div++;
	chomp($line);
	my @intermediate = split(/ /,$line);

		if($div % 3==2)
		{
		push @{$rnaseq[$u-1]}, $intermediate[0];
		push @{$seqno[$u-1]}, $intermediate[1];
		}
		elsif($div % 3==0)
		{
		push @{$chromseq[$u-1]}, $intermediate[0];
		push @{$map[$u-1]}, int($intermediate[1])-1;
		}
		elsif($div % 3==1)
		{
			if($u<=2)
			{
			push @{$cs[$u-1]}, int($intermediate[0])-1;


			}
			if($u>2)
			{
			push @{$cs[$u-1]},int($intermediate[0])-1;
	
			}

		}
	}
}

close(polyAs1filtered);
close(polyAs2filtered);
close(polyAs3filtered);
close(polyAs4filtered);
close(polyAs5filtered);
close(polyAs6filtered);

for(my $u=0; $u<=5; $u++)
{
	for(my $v=0; $v<=$#{$cs[$u]}; $v++)
	{
	my $q=0;
	my $y=$cs[$u][$v]-$map[$u][$v]-1;
		for(my $i=0; $i<=$y; $i++)
		{

			if($chromosomes[$u][$map[$u][$v]+$i] eq substr($rnaseq[$u][$v],$i,1))
			{
			$q++;
			}
			
		}
		$quality[$u][$v]=$q/($cs[$u][$v]-$map[$u][$v]);
		#print $quality[$u][$v]."\n";
	}
}
print "qualITY DONE\n";
my %rnaseqhash;
my %maphash;
my %chromseqhash;
my %cshash;
my %qualityhash;
my %chromosomehash;
for(my $u=0; $u<=5; $u++)
{
	for(my $v=0; $v<=$#{$cs[$u]}; $v++)
	{
#print $seqno[$u][$v]."\n";
	push @{$rnaseqhash{$seqno[$u][$v]}},$rnaseq[$u][$v];
	push @{$maphash{$seqno[$u][$v]}},$map[$u][$v];
	push @{$chromseqhash{$seqno[$u][$v]}},$chromseq[$u][$v];
	push @{$cshash{$seqno[$u][$v]}},$cs[$u][$v];
	push @{$qualityhash{$seqno[$u][$v]}},$quality[$u][$v];
	push @{$chromosomehash{$seqno[$u][$v]}},$u;
	}
}


my @csfiltered;

my $key;
my $size=keys( %maphash );
my $count=0;
foreach $key(keys %maphash)
{
$count++;
my $done=$count/$size;
print $count." ". $done."\n";

	if($#{$maphash{$key}}==0)
	{
	push @{$csfiltered[$chromosomehash{$key}[0]]},$cshash{$key}[0];
	}
	elsif($#{$maphash{$key}}>0)
	{
	my @maxqs;
	my @maxchr;
	my $M=maxval(@{$qualityhash{$key}});
		for(my $i=0; $i<=$#{$maphash{$key}};$i++)
		{
			if($qualityhash{$key}[$i]==$M)
			{
			push @maxqs,$i;
			push @maxchr,$chromosomehash{$key}[$i];
			}

		}
	my @css;
	my @chr;
	my @utr;
		if($#maxqs==0)
		{
		push @{$csfiltered[$maxchr[0]]},$cshash{$key}[$maxqs[0]];
		}
		elsif($#maxqs>0)
		{
		for(my $j=0; $j<=$#maxqs;$j++)
		{
			if($maxchr[$j]<=2)
			{
			for(my $g=0; $g<=$#{$plusfeatures[$maxchr[$j]]}; $g++)
			{
				if($g<$#{$plusfeatures[$maxchr[$j]]})
				{
					if($cshash{$key}[$maxqs[$j]]>$plusfeatures[$maxchr[$j]][$g][3]&&$cshash{$key}[$maxqs[$j]]<$plusfeatures[$maxchr[$j]][$g+1][2])
					{
					push @css, $cshash{$key}[$maxqs[$j]];
					push @chr, $maxchr[$j];
					push @utr, $cshash{$key}[$maxqs[$j]]-$plusfeatures[$maxchr[$j]][$g][3];
					}

				}
				if($g==$#{$plusfeatures[$maxchr[$j]]})
				{
					if($cshash{$key}[$maxqs[$j]]>$plusfeatures[$maxchr[$j]][$g][3])
					{
					push @css, $cshash{$key}[$maxqs[$j]];
					push @chr, $maxchr[$j];
					push @utr, $cshash{$key}[$maxqs[$j]]-$plusfeatures[$maxchr[$j]][$g][3];
					}

				}
			}
			}
			if($maxchr[$j]>2)
			{
			for(my $g=0; $g<=$#{$minusfeatures[$maxchr[$j]-3]}; $g++)
			{
				if($g==0)
				{
				my $stop=$#{$chromosomes[$maxchr[$j]]}-$minusfeatures[$maxchr[$j]-3][$g][2]-1;
					if($cshash{$key}[$maxqs[$j]]>$stop)
					{
					push @css, $cshash{$key}[$maxqs[$j]];
					push @chr, $maxchr[$j];
					push @utr, $cshash{$key}[$maxqs[$j]]-$stop;
					}
				}
				if($g>0)
				{
				my $otherstart=$#{$chromosomes[$maxchr[$j]]}-$minusfeatures[$maxchr[$j]-3][$g-1][3]-1;
				my $stop=$#{$chromosomes[$maxchr[$j]]}-$minusfeatures[$maxchr[$j]-3][$g][2]-1;
					if($cshash{$key}[$maxqs[$j]]>$stop&& $cshash{$key}[$maxqs[$j]]<$otherstart)
					{
					push @css, $cshash{$key}[$maxqs[$j]];
					push @chr, $maxchr[$j];
					push @utr, $cshash{$key}[$maxqs[$j]]-$stop;
					}
				}
			}
			}
		}
		my @mins;
		my $Min;
		if($#utr>0)
		{       
		$Min=minval(@utr);
		}
		elsif($#utr==0)
		{
		$Min=$utr[0];
		}
		
		for(my $i=0; $i<=$#utr;$i++)
		{
			if($utr[$i]==$Min)
			{
			push @mins,$i;
			}

		}
		if($#mins==0)
		{
		##push @{$csfiltered[$chr[0]]},$css[$mins[0]];
		push @{$csfiltered[$chr[$mins[0]]]},$css[$mins[0]];
		}
		elsif($#mins<0)
		{
		print $cshash{$key}[$maxqs[0]]." ".$maxchr[0]."\n";
		}
		else
		{
		print $#mins."\n";
		for(my $i=0; $i<=$#mins;$i++)
		{
		print $mins[$i]." ".$chr[$mins[$i]]." ".$css[$mins[$i]]." ".$utr[$mins[$i]]."\n";
		}
		}
		}
	}
                                                                                                                                                                                                                                                                               
}
my @uncsfiltered;
for(my $u=0; $u<=2; $u++)
{
	@{$uncsfiltered[$u]}=uniq @{$csfiltered[$u]};
	@{$uncsfiltered[$u]}=sort {$a <=> $b} @{$uncsfiltered[$u]};
}
for(my $u=3; $u<=5; $u++)
{
	@{$uncsfiltered[$u]}=uniq @{$csfiltered[$u]};
	@{$uncsfiltered[$u]}=sort {$b <=> $a} @{$uncsfiltered[$u]};
}
my @occ;
for(my $u=0; $u<=5; $u++)
{
	for(my $v=0; $v<=$#{$csfiltered[$u]};$v++)
	{
	$occ[$u][$csfiltered[$u][$v]]=0;
	}
}
for(my $u=0; $u<=5; $u++)
{
	for(my $v=0; $v<=$#{$csfiltered[$u]};$v++)
	{
	$occ[$u][$csfiltered[$u][$v]]++;
	}
}

open(cs0,">", $dataset.'/'.$dataset.'cs1fwdocc') or die("Unable to open cs0");
open(cs1,">", $dataset.'/'.$dataset.'cs2fwdocc') or die("Unable to open cs1");
open(cs2,">", $dataset.'/'.$dataset.'cs3fwdocc') or die("Unable to open cs2");
open(cs3,">", $dataset.'/'.$dataset.'cs1revocc') or die("Unable to open cs3");
open(cs4,">", $dataset.'/'.$dataset.'cs2revocc') or die("Unable to open cs4");
open(cs5,">", $dataset.'/'.$dataset.'cs3revocc') or die("Unable to open cs5");

	for(my $v=0; $v<=$#{$uncsfiltered[0]}; $v++)
	{
	print cs0 $uncsfiltered[0][$v]."\t".$occ[0][$uncsfiltered[0][$v]]."\n";
	}
	for(my $v=0; $v<=$#{$uncsfiltered[1]}; $v++)
	{
	print cs1 $uncsfiltered[1][$v]."\t".$occ[1][$uncsfiltered[1][$v]]."\n";
	}
	for(my $v=0; $v<=$#{$uncsfiltered[2]}; $v++)
	{
	print cs2 $uncsfiltered[2][$v]."\t".$occ[2][$uncsfiltered[2][$v]]."\n";
	}
	for(my $v=0; $v<=$#{$uncsfiltered[3]}; $v++)
	{
	my $dudi=$#{$chromosomes[3]}-$uncsfiltered[3][$v]-1;
	print cs3 $dudi."\t".$occ[3][$uncsfiltered[3][$v]]."\n";
	}
	for(my $v=0; $v<=$#{$uncsfiltered[4]}; $v++)
	{
	my $dudi=$#{$chromosomes[4]}-$uncsfiltered[4][$v]-1;
	print cs4 $dudi."\t".$occ[4][$uncsfiltered[4][$v]]."\n";
	}
	for(my $v=0; $v<=$#{$uncsfiltered[5]}; $v++)
	{
	my $dudi=$#{$chromosomes[5]}-$uncsfiltered[5][$v]-1;
	print cs5 $dudi."\t".$occ[5][$uncsfiltered[5][$v]]."\n";
	}
close(cs0);
close(cs1);
close(cs2);
close(cs3);
close(cs4);
close(cs5);
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
	return ($arraysort[$count/2] + $arraysort[$count/2 + 1]) / 2;

	} 
	else 
	{
	return $arraysort[int($count/2)];
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
