#!/usr/bin/perl
#************************************************Only identified MISSENSE VARIANTS***********************************

#Input: Names of series (Eg: ADSP), Patient Info (FORMAT: PatientID\tGender\tAge\tAD/Control\tAPOE), BioR Annotation File, Total number of samples, Flat file and variant counts output file names.
#USAGE: $perl generate_flat_file_SG.pl ADSP ADSP_sample_sex_age_AD_APOE.txt SORCS2_variants.txt 1000 SORCS2_flat_file.txt SORCS2_allele_counts.txt

use strict;
use warnings;
use Data::Dumper qw(Dumper);
my $series=$ARGV[0];
my $sample_count=$ARGV[3];
my %patient_bio=();
open FILE1, "$ARGV[1]" or die $!;
my $patient_bio_header=<FILE1>;
#print "Header:$patient_bio_header";
while(<FILE1>)
{	
	$_=~s/\n|\r//g;	
	my @p_info=split(/\t/,$_);
	$patient_bio{$p_info[0]}[0]=$p_info[1]; #Gender
	$patient_bio{$p_info[0]}[1]=$p_info[2]; #Age
	$patient_bio{$p_info[0]}[2]=$p_info[3]; #AD/Control
	$patient_bio{$p_info[0]}[3]=$p_info[4]; #APOE
}
close(FILE1);	
#print Dumper \%patient_bio;
open FILE2, "$ARGV[2]" or die $!;	#BioR or Gene Specific BioR Output (after, for eg: $ grep "APOE" Batch1_1000_samples_AD.single_allele.xls > APOE_variants.txt)
open FILE3, ">$ARGV[4]" or die $!;	#Flat file for output
open FILE4, ">$ARGV[5]" or die $!;	#Allele counts per variant output

my $annotation_header=<FILE2>;
$annotation_header=~ s/s_//g;
my @anno_info=split(/\t/,$annotation_header);
print FILE3 "Subject_ID\tSeries\tSex\tAge\tCase\/control\tAPOE\tGene\tSNP_ID\tCHR\tPosition\tRef\tAlt\tCodon\tAA\tCall\tGenotype\n";
print FILE4 "Gene\tSNP_ID\tChr\tPosition\tRef\tAlt\tCodon\tAA\tControl(0/0)\tControl(0/1)\tControl(1/1)\tControl(NA)\tAD(0/0)\tAD(0/1)\tAD(1/1)\tAD(NA)\n";
while(<FILE2>)
{
	my %call=('0/0',11,'0/1',12,'1/1',22);
	#print Dumper \%call;
	$_=~s/\n|\r//g;	
	my @var=split(/\t/,$_);
	if($var[$sample_count+60] eq "MISSENSE")
	{
		my %c_genotype=('0/0',0,'0/1',0,'1/1',0,'NA',0);	#Control	(0)
		my %t_genotype=('0/0',0,'0/1',0,'1/1',0,'NA',0);	#Case	(1,2)
		for(my $i=9;$i<(9+$sample_count);$i++)
		{
			$var[$i]=~ s/\"//g;
			my @gt=split(':',$var[$i]);
#			print $gt[0]."\t";	
			if($patient_bio{$anno_info[$i]}[2]==0)
			{
				if($gt[0]=~/\./)
				{
					$c_genotype{'NA'}++;
				}
				else
				{
					$c_genotype{$gt[0]}++;
				}
			}
			else
			{
				if($gt[0]=~/\./)
				{
					$t_genotype{'NA'}++;
				}
				else
				{
					$t_genotype{$gt[0]}++;
				}
			}
		}
	#	print "\n";
	#	print Dumper \%c_genotype;
	#	print Dumper \%t_genotype;
	#	print Dumper \%call;
		if(($c_genotype{'0/0'})&&($c_genotype{'1/1'}))
		{
			if($c_genotype{'0/0'} < ($c_genotype{'1/1'}*0.3))
			{
#				print Dumper \%call;
#				print Dumper \%c_genotype;
				$call{'0/0'}=22;
				$call{'1/1'}=11;
#				print Dumper \%call;
	
			}
		}
		if(($t_genotype{'0/0'})&&($t_genotype{'1/1'}))
		{
			if($t_genotype{'0/0'} < ($t_genotype{'1/1'}*0.3))
			{
#				print Dumper \%call;
#				print Dumper \%t_genotype;
				$call{'0/0'}=22;
				$call{'1/1'}=11;
#				print Dumper \%call;
			}
		}
	#	print Dumper \%call;
		for(my $i=9;$i<(9+$sample_count);$i++)
		{
			my @gt=split(':',$var[$i]);
			if(!$call{$gt[0]})
			{
				$call{$gt[0]}='NA';	
			}
			print FILE3 "$anno_info[$i]\t$series\t$patient_bio{$anno_info[$i]}[0]\t$patient_bio{$anno_info[$i]}[1]\t$patient_bio{$anno_info[$i]}[2]\t$patient_bio{$anno_info[$i]}[3]\t$var[$sample_count+63]\t$var[$#var]\t$var[0]\t$var[1]\t$var[3]\t$var[4]\t$var[$sample_count+61]\t$var[$sample_count+62]\t$call{$gt[0]}\t$gt[0]\n";
		}
		my @freq=();
		print FILE4 "$var[$sample_count+63]\t$var[$#var]\t$var[0]\t$var[1]\t$var[3]\t$var[4]\t$var[$sample_count+61]\t$var[$sample_count+62]\t";
		foreach my $key (sort keys %c_genotype)
		{
			print FILE4 "$c_genotype{$key}\t";
		}
		my $c_typed=join(",",@freq);
		@freq=();	
		foreach my $key (sort keys %t_genotype)
		{
			print FILE4 "$t_genotype{$key}\t";
		}
		my $t_typed=join(",",@freq);
		print FILE4 "\n";
	}
}
close(FILE2);
close(FILE3);
close(FILE4);
