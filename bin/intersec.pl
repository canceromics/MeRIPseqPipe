#! /usr/bin/perl -w
use strict;
use Getopt::Long;


my $file1;
my $column1;
my $file2;
my $column2;
my $type;
GetOptions(
	'a=s' => \$file1,
	'na=i' => \$column1,
	'b=s' => \$file2,
	'nb=i' => \$column2,
	't=s' => \$type
);
open FH1, $file1 or die "can not open $file1: $!";

my %id1;
while(<FH1>){
	chomp;
	$_=~s/"//g;
	my @field=split /\s+/;
	if(!exists $id1{$field[$column1-1]}){
		$id1{$field[$column1-1]}=$_;
	}else{
		$id1{$field[$column1-1]}.="\n".$_;
	}
}
close(FH1);
open FH2,$file2 or die "can not open $file2:$!";
my %id2;
while(<FH2>){
	chomp;
	$_=~s/"//g;
	my @field=split /\s+/;
	if(!exists $id2{$field[$column2-1]}){
                $id2{$field[$column2-1]}=$_;
        }else{
                $id2{$field[$column2-1]}.="\n".$_;
        }
}

if($type eq "ua"){
	foreach my $a (keys %id1){
		if(!exists $id2{$a}){
			print $id1{$a}."\n";
		}
	}
}

if($type eq "ub"){
	foreach my $b (keys %id2){
                if(!exists $id1{$b}){
                        print $id2{$b}."\n";
                }
        }
}

if($type eq "d"){
        foreach my $b (keys %id2){
                if(exists $id1{$b}){
                        print $b."\n";
                }
        }

}

if($type eq "da"){
        foreach my $b (keys %id2){
                if(exists $id1{$b}){
                        print $id1{$b}."\n";
                }
        }
}
if($type eq "db"){
        foreach my $b (keys %id2){
                if(exists $id1{$b}){
                        print $id2{$b}."\n";
                }
        }
}
if($type eq "dab"){
	foreach my $b (keys %id2){
                if(exists $id1{$b}){
			my @tmp1=split "\n",$id1{$b};
			my @tmp2=split "\n",$id2{$b};
			foreach my $t1 (@tmp1){
				foreach my $t2 (@tmp2){
					print $t1."\t".$t2."\n";
				}
			}
                        #print $id1{$b}."\t".$id2{$b}."\n";
                }
        }
}
