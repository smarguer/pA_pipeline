#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV != 3) {die "wrong number of files";}
(my $in, my $gff, my $tag)=@ARGV;

open (IN, $in) or die 'could not find the input file';
open (IN1, $gff) or die 'could not find the input file';
open (OUT,'>',"$tag.pA.map") or die 'could not open the output file';
  

my $line=<IN>;
my $line1;
my $count=0;
my $count1=0;
my $hit;
my $utr;
my $name;
my @holder;
my @holder1;
my @match;
my %gene;
my %strand;
my %utr;
my %cs;
my %chr;

while ($line1=<IN1>)
 {
  chomp $line1;
  @holder1 = split (/\t/, $line1);

  unless($holder1[2] eq "gene"){next;}
  if($holder1[9]=~/LTR/)
  {
   $name='LTR.'.$holder1[3].'.'.$holder1[12];
   #print "$name\n";
  } 
  else
  {
   $name=$holder1[9];
  }
  $strand{$name}=$holder1[6]; 
  $chr{$name}=$holder1[12]; 
  if($strand{$name} eq '+')
  {
   $gene{$name}=$holder1[4];
  }
  elsif($strand{$name} eq '-')
  {
   $gene{$name}=$holder1[3];
  }
 }

close IN1;

while ($line=<IN>)
 {
  $count++;
  chomp ($line);

  @holder = split (/\t/, $line);
  @match=@holder[5..$#holder];
  #print $#match;
   foreach (@match)
   {
    $hit="NA";
    unless(($_=~/AS\.S/)||($_=~/SPNCRNA/)||($_=~/NEW/))
    {
     $_=~/\.[EI]{1}\d+/;
     $hit=$`;
    }
   }
  #print "$hit\n";
  if ($hit eq "NA")
  {
   #print "$holder[0]\t$holder[1]\t$holder[2]\t$holder[3]\t$holder[4]\t$hit\tNA\tNA\n";
  }
  else
  {
   unless(defined $gene{$hit})
   {
    print "$line\n";
    next;
   }
   if($holder[3] eq "+")
   {
    $utr=$holder[1]-$gene{$hit};
   }
   else
   {
    $utr=$gene{$hit}-$holder[1];
   }
   #print "$holder[0]\t$holder[1]\t$holder[2]\t$holder[3]\t$holder[4]\t$hit\t$gene{$hit}\t$utr\n";
   $utr{$hit}{NUM}++;
   if(!(defined $utr{$hit}{AWAY}) || ($utr>$utr{$hit}{AWAY}))
   {
    $utr{$hit}{AWAY}=$utr;
    $cs{$hit}{AWAY}=$holder[1];
   }
   if(!(defined $utr{$hit}{HIGH}) || ($holder[0]>$utr{$hit}{HIGH}))
   {
    $utr{$hit}{HIGH}=$utr;
    $cs{$hit}{HIGH}=$holder[1];
   }
   if(!(defined $utr{$hit}{CLOSE}) || ($holder[0]<$utr{$hit}{CLOSE}))
   {
    $utr{$hit}{CLOSE}=$utr;
    $cs{$hit}{CLOSE}=$holder[1];
   }
  }
 }

 print OUT "name\torf_end\tchr\tstrand\tmain_utr\tshortest_utr\tlongest_utr\tmain_cs\tproximal_cs\tdistal_cs\tcs_number\n";

 foreach (keys %gene)
 {
  if (defined $utr{$_}{NUM})
  {
   print OUT "$_\t$gene{$_}\t$chr{$_}\t$strand{$_}\t$utr{$_}{HIGH}\t$utr{$_}{CLOSE}\t$utr{$_}{AWAY}\t$cs{$_}{HIGH}\t$cs{$_}{CLOSE}\t$cs{$_}{AWAY}\t$utr{$_}{NUM}\n";
  }
  else
  {
   print OUT "$_\t$gene{$_}\t$chr{$_}\t$strand{$_}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
  }
 }

print "\n";

close IN;
close OUT;

