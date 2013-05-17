#!/usr/bin/perl

use strict;
use warnings;

my $line;
my $line1;
my $count=0;
my @holder;
my @holder1;
my $what=0;
my $chr;
my $strand;
my %CS;
my %DB;
my %strand;

if (@ARGV != 1) {die "wrong number of dataset";}
(my $in)=@ARGV;
my @suffix=("cs1fwdocc","cs1revocc","cs2fwdocc","cs2revocc","cs3fwdocc","cs3revocc");
foreach my $cs (@suffix)
{
 $what++;
 if(($what ==1)||($what ==3)||($what ==5))
  {
   $strand="plus";
  }
  else
  {
   $strand="minus";
  }

 open (CS, 'usedfiles/'.$in.'/'.$in.$cs) or die 'pffff, cannot find the CS count files...';
 $chr=$what;
 if($what < 3){$chr=1;}
 elsif($what < 5){$chr=2;}
 elsif($what < 7){$chr=3;}

 #print "#################$what###$chr#################\n";
 while ($line=<CS>)
 {
  $count++;
  #print "$count\n";

  chomp ($line);
    @holder = split (/\t/, $line);
  
  $CS{$holder[0]}{$chr}{$strand}=$holder[1];
  #print "$holder[0]\t$chr\t$strand\t$CS{$holder[0]}{$chr}{$strand}\n"; 
 }
 close CS;
}

$chr=0;

my @db=("db1","db2","db3");

$strand{"plus"}="+";
$strand{"minus"}="-";

foreach my $db (@db)
{
 $chr++;
 open (DB, 'db/'.$in.'/'.$db.'_'.$in) or die 'pffff, cannot find the CS count files...';
 #print "OK $chr\n"; 

 while ($line1=<DB>)
 {
  $count++;
  chomp ($line1);
  @holder1 = split (/\t/, $line1);
  unless($holder1[1]=~/CS/){next;}  
  #unless($CS{$holder1[0]}{$chr}{$holder1[3]})
  #{
  # print "$line1\n";
  #}
 
  #print "$holder1[0]\t";
  print 'vulgar: '.$CS{$holder1[0]}{$chr}{$holder1[3]}.' NA NA NA CH'.$chr.'_bases '.$holder1[0].' '.$holder1[0].' '.$strand{$holder1[3]}.' NA NA NA'."\n";
  #print "$CS{$holder1[0]}{$chr}{$strand}\n";
  #print "$CS{$holder1[0]}{$chr}\n";
  #print "$CS{$holder1[0]}{$chr}{$holder1[3]}\n";
 }


 close DB;
}



