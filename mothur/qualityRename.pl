#!/bin/env perl
# parse log file from mothur and rename files for input into analysis
#

open FILE, "<", "mothur.quality" or die "could not open mothur.quality: ".$!;
open OUTLOG, ">", "input.shared" or die "could not open input.shared for writing: ".$!;

my $log = 0;
my $header = 1;
my $names = 0;

while(<FILE>){
   if($header == 1){
	print OUTLOG $_;
	if($_ =~ /Batch\ Mode/){
		$header=0;
		print "resetting $header\n";
	}
        next;
   }
   elsif($log){
      $log=0 if ($_ =~ /^\n/);
      chomp;
      my ($variable,$value) = split/\=/, $_;
      my @suffixes = split /\./, $value;
      $value = "quality.$suffixes[-1]";
      print OUTLOG "$variable=$value\n" if (defined $variable);
      next;
   }
   elsif($names){
	if ($_ =~ /^\w/){
		chomp;
		my @temp = split /\./, $_;
		print OUTLOG "quality.$temp[-1]\n";
	}
	else {
		print OUTLOG "\n";
		$names=0;
	}
	next;
   }

   if($_ =~ /Output\ File\ Names/i){
	$names=1;
   }
   if ($_ =~ /Current\ files\ saved\ by\ mothur\:/){
	$log = 1;
   }
   
   print OUTLOG $_;
   next;
   
}


