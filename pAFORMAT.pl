#!/usr/bin/perl

use strict;
use warnings;

my $line;
my $line1;
my $line2;
my $count=0;
my $count1=0;
my $count2=0;
my $count3=0;
my @holder;
my $up;
my $dwn;
my $name;
my $lchr1=5579133;
my $lchr2=4539804;
my $lchr3=2452883;


my %chr;
my %pAreads;
my %pAcount;
my $holder;
my $read;
my $mapStart;
my $mapEnd;
my @seg;
my @gen=();
my $hitStart;
my $cs;

if (@ARGV != 4) {die "wrong number of files";}
(my $in,my $in1, my $in2, my $tag)=@ARGV;

open (IN, $in) or die 'could not find the input file';
open (IN1, $in1) or die 'could not find the input file';
open (SEQ, $in2) or die 'could not find fasta sequences';

open (OUT1,">","$tag.all") or die 'could not create "all" file';
open (OUT2,">","$tag.rita") or die 'could not create "rita" file';
open (OUT3,">","$tag.paCount") or die 'could not create "paCount" file';
open (OUT4,">","$tag.ALLpA") or die 'could not create "ALLpA" file';
open (OUT5,">","$tag"."1fwdforfilter3") or die 'could not create "Rita" file';
open (OUT6,">","$tag"."2fwdforfilter3") or die 'could not create "Rita" file';
open (OUT7,">","$tag"."3fwdforfilter3") or die 'could not create "Rita" file';
open (OUT8,">","$tag"."1revforfilter3") or die 'could not create "Rita" file';
open (OUT9,">","$tag"."2revforfilter3") or die 'could not create "Rita" file';
open (OUT10,">","$tag"."3revforfilter3") or die 'could not create "Rita" file';


#####################
#parse sequence file#
#####################
while ($line2=<SEQ>)
{
 $count++;
 chomp($line2);
 if ($line2=~/\>/)
 {
  $name=substr $line2,1;
  $chr{$name}="A";
 }
 else
 {
  $chr{$name}.=$line2;
 }
}

$chr{"CH1_bases"}=substr($chr{"CH1_bases"},1);
$chr{"CH2_bases"}=substr($chr{"CH2_bases"},1);
$chr{"CH3_bases"}=substr($chr{"CH3_bases"},1);
$chr{"CH4_bases"}=substr($chr{"CH4_bases"},1);
$chr{"CH5_bases"}=substr($chr{"CH5_bases"},1);
$chr{"CH6_bases"}=substr($chr{"CH6_bases"},1);

close SEQ;

#########################
#reads 'ALLtrs.left' file
#########################
while ($line=<IN>)
 {
  $count++;
  chomp ($line);
  if($line=~/>/)
  {
   $holder=substr $line,1;
  }
  elsif($line=~/(A{3}A+$)/)
  {
   $count1++;
   $pAreads{$holder}{SEQ}=$line;
   $pAreads{$holder}{PA}=length($1);
  }
 }
close IN;

##########################
#reads 'ALLchr' file
##########################

while ($line1=<IN1>)
{

  @holder = split (/ /, $line1);
  if ($holder[0] !~ /vulgar/){next;}
  unless($pAreads{$holder[1]}{SEQ}){next;}

  @seg=($holder[6],$holder[7]);
  @seg = sort {$a <=> $b} @seg;

  $read=$pAreads{$holder[1]}{SEQ};
  $mapStart=lc(substr $read,$holder[2],1);
  $mapEnd=lc(substr $read,($holder[3]-1),1);
  $up= substr $read,$holder[2],1,$mapStart;  
  $dwn= substr $read,($holder[3]-1),1,$mapEnd;
  
  $hitStart=$holder[6];
  $cs=$holder[7];

############
 
# if($pAreads{$holder[1]}{PA}<=51-($holder[3]-1))
  if((($pAreads{$holder[1]}{PA} + $holder[3]) == 51)||(($pAreads{$holder[1]}{PA} + $holder[3]) == 52)) 
  {
   $count3++; 
   $pAcount{$holder[5]}{$holder[8]}{$holder[7]}{COUNT}++;
   print OUT4 "$line1"; 
   
   print OUT1 ">$holder[1]\n$pAreads{$holder[1]}{SEQ}\n$read\n";
   @gen=(substr($chr{$holder[5]},$seg[0],($seg[1]-$seg[0])),substr($chr{$holder[5]},($seg[0]-$holder[2]),51));
   my $N="----------------------------------";
   
   if ($holder[8] eq '-')
   {
    $gen[1]=substr($chr{$holder[5]},$seg[0]-(51-($seg[1]-$seg[0])-$holder[2]),51);
    $gen[0]=~tr/natcg/ntagc/;
    $gen[0]=reverse $gen[0];
    $gen[1]=~tr/natcg/ntagc/;
    $gen[1]=reverse $gen[1];
   }
   
   if($holder[2] != 0)
   {
    $gen[0]=substr($N,1,$holder[2]).$gen[0];
   }
   
   print OUT1 "$gen[0]\n$gen[1]\n";

#####print files for Rita's pipeline

   $read=uc($read);
   $gen[1]=uc($gen[1]);

   if($holder[8] eq "+")
   {
    if($holder[5]=~/1/)
    {
     #print OUT5 "$hitStart\n$read $count3\n$gen[1] ".($cs+1)."\n";
     print OUT5 ($cs+1)."\n$read $count3\n$gen[1] ".($hitStart+1)."\n";
    }
    elsif($holder[5]=~/2/)
    {
     #print OUT6 "$hitStart\n$read $count3\n$gen[1] ".($cs+1)."\n";
     print OUT6 ($cs+1)."\n$read $count3\n$gen[1] ".($hitStart+1)."\n";
    }
    elsif($holder[5]=~/3/)
    {
     #print OUT7 "$hitStart\n$read $count3\n$gen[1] ".($cs+1)."\n";
     print OUT7 ($cs+1)."\n$read $count3\n$gen[1] ".($hitStart+1)."\n";
    }
   }
   elsif ($holder[8] eq "-")
   {
    if($holder[5]=~/1/)
    {
     #print OUT8 (($lchr1-$hitStart+1)+1)."\n$read $count3\n$gen[1] ".(($lchr1-$cs+1)+1)."\n";
     print OUT8 (($lchr1-$cs+1)+1)."\n$read $count3\n$gen[1] ".(($lchr1-$hitStart+1)+1)."\n";
    }
    elsif($holder[5]=~/2/)
    {
     #print OUT9 (($lchr2-$hitStart+1)+1)."\n$read $count3\n$gen[1] ".(($lchr2-$cs+1)+1)."\n";
     print OUT9 (($lchr2-$cs+1)+1)."\n$read $count3\n$gen[1] ".(($lchr2-$hitStart+1)+1)."\n";
    }
    elsif($holder[5]=~/3/)
    {
     #print OUT10 (($lchr3-$hitStart+1)+1)."\n$read $count3\n$gen[1] ".(($lchr3-$cs+1)+1)."\n";
     print OUT10 (($lchr3-$cs+1)+1)."\n$read $count3\n$gen[1] ".(($lchr3-$hitStart+1)+1)."\n";
    }
   }
###################################
  }
  
  $up=undef;
  $dwn=undef;
}
 
foreach my $out (keys %pAcount)
{
 foreach my $out1 (keys %{$pAcount{$out}})
 {
  foreach my $out2 (keys %{$pAcount{$out}{$out1}})
  {
   print OUT3 "vulgar: ".$pAcount{$out}{$out1}{$out2}{COUNT}." NA NA NA $out $out2 $out2 $out1 NA NA NA\n"; 
   $count2=$count2 + $pAcount{$out}{$out1}{$out2}{COUNT};
  }
 }
}

print "$count\n$count1\n$count2\n";

close IN1;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
close OUT7;
close OUT8;
close OUT9;
close OUT10;

#bahler-pc7956: ~ $perl -e '$test="AAAAA"; print "$test\n"; $test1=substr $test,3,1,"a"; print "$test1\n$test\n";'
#AAAAA
#A
#AAAaA
#>@HWUSI-EAS483_0008_FC:1:3:12141:17448#0/1;33.0   
#>@HWUSI-EAS483_0008_FC:1:1:15188:12046#0/1;52.8 
############

