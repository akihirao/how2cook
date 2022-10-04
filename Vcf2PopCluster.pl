#!/usr/bin/perl
# Vcf2PopCluster.pl
# conversion script from vcf to PopCluster (Wang 2022)
# by HIRAO Akira
# how to use: perl Vcf2PopCluster.pl <- hogehoge.vcf > hogehoge.dat
# The data format for PopCluster is individual genotypes in 2 rrows.


@indiv_genotype = ();

$no_locus = 0;
while ($line = <>) {
	chomp $line;

	if($line =~ m/^#CHROM/){
		($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @INDV)= split /\s+/, $line;
		$no_indiv = @INDV;
	}


	if($line !~ m/^#/){
		$target_loc_genotype_across_indiv = "";

		($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @G_INFO)= split /\s+/, $line;

		foreach $genotype_info (@G_INFO){
			($genotype, @other_info) = split /:/, $genotype_info;
			$allele_A = substr($genotype, 0, 1);
			$allele_B = substr($genotype, 2, 1);

			if($allele_A =~ /[0-1]/ && $allele_B =~ /[0-1]/){
				$genotype_out = $allele_A + $allele_B;
			}else{
				$genotype_out = 3;
			}

			$target_loc_genotype_across_indiv = $target_loc_genotype_across_indiv.$genotype_out;

		}

		$ref_indiv_genotype = \@indiv_genotype;

		for($i = 0; $i < $no_indiv; $i++) {
			$target_genotype = substr($target_loc_genotype_across_indiv, $i, 1);

			$pre_indiv_genotype = $$ref_indiv_genotype[$i];

			$post_indiv_genotype= $pre_indiv_genotype.$target_genotype;

			$$ref_indiv_genotype[$i] = $post_indiv_genotype;

		}

	}
}



#Output PopCluster input form
for($j = 0; $j < $no_indiv; $j++) {
	print $INDV[$j], "\n";
	print $$ref_indiv_genotype[$j], "\n";
}
