#!/usr/bin/perl
# make_annotation_file_DDBJ.pl

# how to use: make_annotation_file_DDBJ.pl hogehoge.fasta hogehoge.bed > hogehoge.ann.txt
# hogehoge.fasta: a fasta file of genome assembly
# hogehoge.bed: a bed file of assembly gap
# hogehoge.ann.txt: a annotation text file as output


# input files
$entry_fasta = $ARGV[0]; # a fasta file of genome assembly
$gap_bed_file = $ARGV[1]; # a bed file of assembly gap


# items
@ab_name = ("Hogehoge", "Hogehogehoge");#authors
$contact = "Hogehoge";# contact person
$email = "hogehoge"."\@"."fra.go.jp"; #email for contact person
$phone = "xx-xx-xx-xx";
$institute = "Fisheries Resources Research Institute";
$country = "Japan";
$city = "Yokohama";
$street = "2-12-4, Fukuura, Kanazawa";
$zip = "236-8648";

$ref_title = "Chu-chu tako kaina";
@ref_ab_name = ("Hogehoge", "Hogehogehoge");
$ref_year = "202x";
$ref_status = "Unpublished";

$organism = "Trachurus japonicus";
$mol_type = "genomic DNA";
$collection_date = "2019";
$collection_country = "Japan";

$hold_date = "yyyymmdd"; # e.g. 20250131


# entry_length
@entry_length = ();
open(FASTA, $entry_fasta);
$line = <FASTA>;
chomp $line;
$chr_lab = $line;
$chr_lab =~ s/>//;
$chr_lab =~ s/sca/chr/;
@sequence = ();
@entry_names = ();

while($line = <FASTA>){
	chomp $line;
	if($line !~ /^>/){
		push(@sequence, $line);
	}else{
		push(@entry_names,$chr_lab);
		$sequence_out = join("", @sequence);
		$seq_length = length($sequence_out);
		push(@entry_length,$seq_length);
		$chr_lab = $line;
		$chr_lab =~ s/>//;
		$chr_lab =~ s/sca/chr/;
		@sequence = ();
	}
}
$sequence_out = join("", @sequence);
$seq_length = length($sequence_out);
push(@entry_length,$seq_length);
push(@entry_names,$chr_lab);
close (FASTA);



# hash gap position
@gap_id_vec = ();
$gap_id = 1;
open(BED, $gap_bed_file);
while($line = <BED>){
	push(@gap_id_vec,$gap_id);
	chomp $line;
	($chr_gap, $start_gap, $end_gap,@bed_info) = split /\s/, $line;
	$chr_gap_hash{$gap_id} = $chr_gap;
	$start_gap_hash{$gap_id} = $start_gap;
	$end_gap_hash{$gap_id} = $end_gap;
	$gap_id = $gap_id + 1;
}
close(BED);
$no_gap = $gap_id;
push(@gap_id_vec, $gap_id);





# output
print "Entry", "\t", "Feature", "\t", "Location", "\t", "Qualifier", "\t", "Value", "\n";

# print ab_names
$ab_name_1 = shift (@ab_name);
print "COMMON", "\t", "SUBMITTER", "\t",  "\t", "ab_name", "\t", $ab_name_1, "\n";
foreach $ab_name_out (@ab_name){
	print "\t",  "\t", "\t",  "ab_name", "\t", $ab_name_out, "\n";
}

print "\t",  "\t", "\t",  "contact", "\t", $contact, "\n";
print "\t",  "\t", "\t",  "email", "\t", $email, "\n";
print "\t",  "\t", "\t",  "phone", "\t", $phone, "\n";
print "\t",  "\t", "\t",  "institute", "\t", $institute, "\n";
print "\t",  "\t", "\t",  "country", "\t", $country, "\n";
print "\t",  "\t", "\t",  "city", "\t", $city, "\n";
print "\t",  "\t", "\t",  "street", "\t", $street, "\n";
print "\t",  "\t", "\t",  "zip", "\t", $zip, "\n";
print "\t",  "REFERENCE",  "\t", "\t",  "title", "\t", $ref_title, "\n";

# print ref_ab_names
foreach $ref_ab_name_out (@ref_ab_name){
	print "\t",  "\t", "\t", "ab_name", "\t", $ref_ab_name_out, "\n";
}
print "\t",  "\t", "\t",  "year", "\t", $ref_year, "\n";
print "\t",  "\t", "\t",  "status", "\t", $ref_status, "\n";
print "\t",  "DATE", "\t", "\t", "hold_date", "\t", $hold_date, "\n";

$no_entry = @entry_names;



$entry_count = 0;
foreach $each_entry (@entry_names){
	$each_entry_length = $entry_length[$entry_count];
	print $each_entry, "\t", "source", "\t", "1..", $each_entry_length, "\t", "organism", "\t", $organism, "\n";
	print "\t",  "\t", "\t",  "mol_type", "\t", $mol_type, "\n";
	print "\t",  "\t", "\t",  "collection_date", "\t", $collection_date, "\n";
		print "\t",  "\t", "\t",  "country", "\t", $collection_country, "\n";
	if($each_entry =~ /chr/){
		$chr_no = $each_entry;
		$chr_no =~ s/chr//;
		print "\t",  "\t", "\t",  "chromosome", "\t", $chr_no, "\n";
	}

	foreach $each_gap_id (@gap_id_vec){
		$target_entry = $chr_gap_hash{$each_gap_id};
		if($each_entry eq $target_entry){
			$target_start = $start_gap_hash{$each_gap_id} + 1;
			$target_end = $end_gap_hash{$each_gap_id};
			$gap_length = $target_end - $target_start;
			print "\t", "assembly_gap", "\t", $target_start,"..",$target_end, "\t", "estimated_length","\t", "unknown", "\n";
			print "\t", "\t", "\t", "gap_type","\t", "within scaffold", "\n";
			print "\t", "\t", "\t", "linkage_evidence","\t", "paired-ends", "\n";

		}
	}
	$entry_count = $entry_count + 1;
}


