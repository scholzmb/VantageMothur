#!/bin/env perl
# parse log file from mothur and rename files for input into analysis
#

open FILE, "<", "mothur.quality" or die "could not open mothur.quality: ".$!;
open OUTLOG, ">", "input.shared" or die "could not open input.shared for writing: ".$!;

my $log = 0;
my @outString;

while(<FILE>){
   if($log){
      $log=0 if ($_ =~ /^\n/);
      chomp;
      my ($variable,$value) = split/\=/, $_;
      next unless (defined $value && $variable);
      if($variable !~ /processors/){
        my @suffixes = split /\./, $value;
	$cmd ="ln -s $value quality.$suffixes[-1]";
	system $cmd;
	$value = "quality.$suffixes[-1]";
      }
      push @outString, "$variable=$value";
      #print OUTLOG "$variable=$value\n" if (defined $variable);
      next;
   }
   if ($_ =~ /Current\ files\ saved\ by\ mothur\:/){
	$log = 1;
   }
   
   next;
   
}

my $string = join ",", @outString;
my $cmd = "sed s/XXX/$string/g ../mothur/shared.template.m > batch_mothur_shared.inputs.m";
system($cmd);

