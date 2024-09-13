#!/usr/bin/perl
#Vcf2allele_depth_count_hetero.pl
#by HIRAO Akira
# how to use: ./vcf2allele_depth_count_het.pl < hogehoge.vcf > hogehoge.allele_depth_ratio.txt


while ($line = <>) {
	chomp $line;
	@allele_depth_ratio_locus = ();

	if($line =~ m/^#CHROM/){

		($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@sample_id) = split /\s+/, $line; 
		$n_sample_id = @sample_id;
		$last_id = pop(@sample_id);
		foreach $out_sample_id (@sample_id){
			print $out_sample_id, "\t";
		}
		print $last_id, "\n";

	}elsif($line !~ m/^#/){	    

		($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@info) = split /\s+/, $line; 
		$n_sample = @info;
		
		if($n_sample==$n_sample_id){


			for($i = 0; $i < $n_sample; $i++){
				$each_sample_info = $info[$i];
				(@sample_info) = split /:/, $each_sample_info;
				$sample_GT = $sample_info[0];
				$sample_DP = $sample_info[1];
				$sample_AD = $sample_info[2];
				@sample_allele_depth = split /,/, $sample_AD;

				if($sample_GT ne "0/1"){
					$allele_depth_ratio_hetero_GT = "NA";	
				}else{
					$allele_depth_ratio_hetero_GT = $sample_allele_depth[1]/$sample_DP;
				}	
			
				push (@allele_depth_ratio_locus, $allele_depth_ratio_hetero_GT);
			}

			$last_allele_depth_ratio_locus =  pop(@allele_depth_ratio_locus);
			foreach $out_allele_depth_ratio (@allele_depth_ratio_locus){
				if($out_allele_depth_ratio ne "NA"){
					printf ("%.3f\t", $out_allele_depth_ratio);	
				}else{
					print "NA", "\t";
				}
			
			}
			if($last_allele_depth_ratio_locus ne "NA"){
					printf ("%.3f\n", $last_allele_depth_ratio_locus);	
			}else{
				print "NA", "\n";
			}
		}
	}
}
