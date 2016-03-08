#!/bin/env perl 
#


open FILE, "$ARGV[0]" or die;

while(<FILE>){
  chomp;
  my @split = split /\t/, $_;
  if($split[1] > 0 ){
	print "$split[0]\n";
  }

}
