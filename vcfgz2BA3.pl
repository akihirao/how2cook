#!/bin/perl
# vcfgz2BA3_input.pl
# Converter from vcf.gz to a input file for BayesAss3 (BA3)
# 
# how to use:
# I. Output With population-indexing 
# i) perl vcf2BA3_input.pl hogehoge.vcf.gz ID_Pop.map > hogehoge.BA3
#
# A map file provides links between individual IDs and population IDs. For example,
#
# ID001	Pop1
# ID002	Pop1
# ID003	Pop2
# ...
#
# II. Output Without population-indexing
# i) perl vcf2BA3_input.pl hogehoge.vcf.gz > hogehoge.BA3
# After the above process, users can execute population-indexing with using ugnix 
# ugnix: https://github.com/brannala/ugnix/wiki/VCF-to-BA3-File-Conversion


use Compress::Zlib;


$no_argv=@ARGV;
#------------------------------------------
# ID-Population

if($no_argv== 2){

	$FILE_MAP_name = $ARGV[1];
	print $FILE_MAP_name, "\n";
	open (FILE_MAP, $FILE_MAP_name) or die "Failed to open a map file\n";

	%id_pop_list = ();

	while ($line = <FILE_MAP>){
		chomp $line;
		($ID, $Population)= split /\s+/, $line;
		$id_pop_list{$ID} = $Population;
	}

	close(FILE_MAP);
}
#------------------------------------------


#------------------------------------------
$temp_pop_ID = 'proxypop';
$FILE_VCF_name = $ARGV[0];

#open (FILE_VCF, $FILE_VCF_name) or die "Failed to open a vcf file\n";
$vcfgz = gzopen($FILE_VCF_name, "rb") or die "File not found.";

#while ($line = <FILE_VCF>){

while($vcfgz -> gzreadline($line) > 0) {
	chomp $line;

	if($line =~ m/^#CHROM/){
		($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @INDV)= split /\s+/, $line;
		$no_indiv = @INDV;
	}


	if($line !~ m/^#/){
		$ID_count = 0;

		($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @G_INFO)= split /\s+/, $line;
		$underbar_lab="_";
		$loc_name = $CHROM.$underbar_lab.$POS;

		foreach $genotype_info (@G_INFO){
			($genotype, @other_info) = split /:/, $genotype_info;
			$allele_A = substr($genotype, 0, 1);
			$allele_B = substr($genotype, 2, 1);

			if($allele_A eq "0"){
				$allele_A_out = $REF;
			}elsif($allele_A eq "1"){
				$allele_A_out  = $ALT;
			}elsif($allele_A eq "." ){
				$allele_A_out  = 0;
			}else{
				$allele_A_out  = 0;
			}


			if($allele_B eq "0"){
				$allele_B_out = $REF;
			}elsif($allele_B eq "1"){
				$allele_B_out  = $ALT;
			}elsif($allele_B eq "." ){
				$allele_B_out  = 0;
			}else{
				$allele_B_out  = 0;
			}

			$target_indiv_ID = $INDV[$ID_count];

			if($no_argv== 2){
				$pop_ID_output = $id_pop_list{$target_indiv_ID};
			}else{
				$pop_ID_output = $temp_pop_ID;
			}

			print $target_indiv_ID, "\t", $pop_ID_output, "\t", $loc_name, "\t", $allele_A_out, "\t",$allele_B_out, "\n";
			$ID_count = $ID_count + 1;

		}
	}

}

$vcfgz -> gzclose();
