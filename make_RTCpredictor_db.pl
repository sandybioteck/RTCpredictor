$usage = <<EOF;
	This script makes database files for RTCpredictor. The script is computationally demanding. To give an estimate, using hg38 genome and annotation files from Ensembl v105, and providing 8 cores/processes, it may take up to 16 hrs to finish while consuming ~64 GB RAM.

	The script can be run using following command:
	perl make_RTCpredictor_db.pl -a genes_transcripts_exons.tsv -p paralogs.tsv -g genome.fa -d distance-in-bp [Optional parameters]

	Required parameters
	-a annotation file named "genes_transcripts_exons.tsv". Please refer to README.txt file for the instructions on how to make this file for the genome of interest.
	-p paralog file named "paralogs.tsv". Please refer to README.txt file for the instructions on how to make this file for the genome of interest.
	-g genome fasta file named "genome.fa". Please refer to README.txt file for the instructions on how to download this file for the genome of interest.
	-d distance (in base pairs. for example 70000 means 70KB). This is the maximum distance allowed between 5' and 3' breakpoints of a read-through chimeric RNA which can be predicted by RTCpredictor.

	Optional parameters
	-c number of cores to use. (default: 1; allowed values between 1 and 8). This parameter can be used to make blast search step faster. However some blast versions only accept this value between 1 and 8.
	-o option to add SNP variants in the database. Enter 1 to add SNP variants and 0 to skip it. (default: 0) (If entered 1 then -s and -v parameters are required).
	-s dbSNP bed file. Please refer to README.txt file for the instructions on how to download this file for the genome of interest.
	-v VCF file. Please refer to README.txt file for the instructions on how to download this file for the genome of interest.
	-h  Prints help and provides instructions on running the script

EOF

$whichblastn = `which blastn`; chomp($whichblastn); if(!(-e $whichblastn)) { print "\n\n\tError! blastn not found. Please install blast+ (2.7.1+ or later) before running. README.txt file has instructions on how to install blastn \n"; exit; }
$whichbedtools = `which bedtools`; chomp($whichbedtools); if(!(-e $whichbedtools)) { print "\n\n\tError! bedtools not found. Please install bedtools (2.27.1 or later) before running. README.txt file has instructions on how to install bedtools \n"; exit; }
use Getopt::Std;
getopts('a:p:g:d:c:o:s:v:h');
$pa=$opt_a; $pp=$opt_p; $pg=$opt_g; $pd=$opt_d; $pc=$opt_c; $po=$opt_o; $ps=$opt_s; $pv=$opt_v; $ph=$opt_h;
if($pa eq "" and $pp eq "" and $pg eq "" and $pd eq "" and $pc eq "" and $po eq "" and $ps eq "" and $pv eq "") { print "$usage"; exit; }
if($pa eq "" and $pp eq "" and $pg eq "" and $pd eq "" and $pc eq "" and $po eq "" and $ps eq "" and $pv eq "" and $ph ne "") { print "$usage"; exit; }

if($pa eq "") { print "$usage\n\n\tError! -a parameter is missing.\n"; exit; }
if($pp eq "") { print "$usage\n\n\tError! -p parameter is missing.\n"; exit; }
if($pg eq "") { print "$usage\n\n\tError! -g parameter is missing.\n"; exit; }
if($pd eq "") { print "$usage\n\n\tError! -d parameter is missing.\n"; exit; }
if($pd !~ /\d+/) { print "$usage\n\n\tError! -d $pd is invalid. It should be a number\n"; exit; }
if($pd < 1) { print "$usage\n\n\tError! -d $pd is invalid. It should be a number >0\n"; exit; }
if($pc eq "") { $pc = 1; }
if($pc !~ /\d+/) { print "$usage\n\n\tError! -c $pc is invalid. It should be a number between 1 and 8\n"; exit; }
if($pc < 1 or $pc > 8) { print "$usage\n\n\tError! -c $pc is invalid. It should be a number between 1 and 8\n"; exit; }
if($po eq "") { $po = 0; }
if($po !~ /\d+/) { print "$usage\n\n\tError! -o $po is invalid. It should be a number 1 or 0\n"; exit; }
if($po == 1)
{
	if($ps eq "") { print "$usage\n\n\tError! -s parameter is missing. With -o 1, -s and -v parameters are required\n"; exit; }
	if($pv eq "") { print "$usage\n\n\tError! -v parameter is missing. With -o 1, -s and -v parameters are required\n"; exit; }
}
#file checks
if(!(-e $pa)) { print "-a $pa file does not exist. Exiting...\n"; exit; }
if(!(-e $pp)) { print "-p $pp file does not exist. Exiting...\n"; exit; }
if(!(-e $pg)) { print "-g $pg file does not exist. Exiting...\n"; exit; }
if($po == 1)
{
	if(!(-e $ps)) { print "-s $ps file does not exist. Exiting...\n"; exit; }
	if(!(-e $pv)) { print "-v $pv file does not exist. Exiting...\n"; exit; }
}

$genes_transcripts_exons = $pa;
$para = $pp;
$genomefa = $pg;
$distance_cut_off = $pd;
$cpu = $pc;
$optsnp = $po;
$dbsnpfile = $ps;
$allvcffile = $pv;
$bedtools = "bedtools";
$makeblastdb = "makeblastdb";
$blastn = "blastn";
##################################################
$cutcols = "1,2,3,10,11,12";

@tmparr = (0..100000);
$flag = 0;
for($r=1; $r<=1000; $r++)
{
	$ran = int(rand(@tmparr));
	$tmpdir = "tmp$ran";
	if(-d $tmpdir) { next; }
	else { $flag = 1; last; }
}
if($flag==0)
{
	print "\tSomething wrong while making temporary directory!! Exiting!\n";
	exit;
}
print "Making temporary directory ./$tmpdir\n";
system "mkdir $tmpdir";

system "cut -f1,11,12 $genes_transcripts_exons |grep -v 'Chromosome' > $tmpdir/123; cut -f3 $genes_transcripts_exons |grep -v 'Gene' > $tmpdir/4; cut -f2 $genes_transcripts_exons |grep -v 'Strand' | perl -p -e 's/-1/-/g' | perl -p -e 's/1/+/g' > $tmpdir/6; cut -f4 $genes_transcripts_exons |grep -v 'Gene' > $tmpdir/7";
$wcl = `cat $tmpdir/123 | wc -l `; chomp($wcl);
system "echo \$(seq -s \",\" 1 $wcl)|perl -p -e 's/,/\\n/g' |perl -p -e 's/^.+/./g' > $tmpdir/5";
system "paste $tmpdir/123 $tmpdir/4 $tmpdir/5 $tmpdir/6 $tmpdir/7 | sort -u | sort -k1,1 -k2,2n > $tmpdir/ifp";
system "rm $tmpdir/123 $tmpdir/4 $tmpdir/5 $tmpdir/6 $tmpdir/7";

print "Processing gene pairs capable of generating read-throughs with breakpoints within <=$distance_cut_off bp ...\n";
$mainfile = "$tmpdir/ifp";
$outfile = "$tmpdir/out.pairs";
`sort -k1,1 -k2,2n $mainfile > $tmpdir/tmp.bed`;
system "rm $tmpdir/ifp";

open(FF, "<$tmpdir/tmp.bed");@fileabc = <FF>; close(FF);
open(OOP, ">$tmpdir/pos.genes.bed");
open(OON, ">$tmpdir/neg.genes.bed");
for($i=0; $i<@fileabc; $i++)
{
	$lineabc = $fileabc[$i]; chomp($lineabc);
	@a = split(/\t/, $lineabc);
	if($a[5] eq "+") { print OOP "$lineabc\n"; }
	if($a[5] eq "-") { print OON "$lineabc\n"; }
}
close OOP; close OON;
`rm $tmpdir/tmp.bed`;

%hash = ();
%hash_track = ();
open(FF, "<$tmpdir/pos.genes.bed");
@file = <FF>; close FF;
`rm $tmpdir/pos.genes.bed`;
$wcl = @file;
for($i=0; $i<@file; $i++)
{
	$ii = $i + 1;
	if($ii == $wcl) { last; }
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$chr = $a[0]; $start = $a[1]; $end = $a[2]; $ensid = $a[3]; $col5 = $a[4]; $strand = $a[5];
	for($j=$ii; $j<@file; $j++)
	{
		$line2 = $file[$j]; chomp($line2);
		@a2 = split(/\t/, $line2);
		$chr2 = $a2[0]; $start2 = $a2[1]; $end2 = $a2[2]; $ensid2 = $a2[3]; $col5_2 = $a2[4]; $strand2 = $a2[5];
		if($ensid eq $ensid2) { next; }
		if($chr ne $chr2) { next; }
		if(exists($hash_track{$ensid}{$ensid2})) { next; }
		$diff = -1000; $diff = $start2 - $end;
		if($diff > $distance_cut_off) { last; }
		if($diff < 0) { next; }
		$hash_track{$ensid}{$ensid2}++;
		$hash{$ensid}{$ensid2}{$diff}++;
	}
}
open(FF, "<$tmpdir/neg.genes.bed");
@file2 = <FF>; close FF;
`rm $tmpdir/neg.genes.bed`;
$wcl = @file2;
for($i=0; $i<@file2; $i++)
{
	$ii = $i + 1;
	if($ii == $wcl) { last; }
	$line = $file2[$i]; chomp($line);
	@a = split(/\t/, $line);
	$chr = $a[0]; $start = $a[1]; $end = $a[2]; $ensid = $a[3]; $col5 = $a[4]; $strand = $a[5];
	for($j=$ii; $j<@file2; $j++)
	{
		$line2 = $file2[$j]; chomp($line2);
		@a2 = split(/\t/, $line2);
		$chr2 = $a2[0]; $start2 = $a2[1]; $end2 = $a2[2]; $ensid2 = $a2[3]; $col52 = $a2[4]; $strand2 = $a2[5];
		if($ensid eq $ensid2) { next; }
		if($chr ne $chr2) { next; }
		if(exists($hash_track{$ensid}{$ensid2})) { next; }
		$diff = -1000; $diff = $start2 - $end;
		if($diff > $distance_cut_off) { last; }
		if($diff < 0) { next; }
		$hash_track{$ensid}{$ensid2}++;
		$hash{$ensid2}{$ensid}{$diff}++;
	}
}
open(OO, ">$outfile");
foreach $k (sort keys(%hash))
{
	foreach $k2 (sort keys(%{$hash{$k}}))
	{
		$k1k2 = "$k\t$k2";
		foreach $d (sort {$a <=> $b} keys(%{$hash{$k}{$k2}}))
		{
			print OO "$k1k2\n";
		}
	}
}
close OO;

print "Implementing Paralog filter ...\n";
$pairs = "$tmpdir/out.pairs";
open(FF, "<$para");
@file = <FF>; close(FF);
%hashp = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	if($line =~ /Gene/) { next; }
	@a = split(/\t/, $line);
	$g1 = $a[0]; $g2 = $a[1];
	if($g1 eq "" or $g2 eq "") { next; }
	$hashp{$g1}{$g2}++;
	$hashp{$g2}{$g1}++;
}

open(FF, "<$pairs");
@file = <FF>; close(FF);
open(OO, ">$tmpdir/out.pairs.nopara");
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$g1 = $a[0]; $g2 = $a[1];
	if(exists($hashp{$g1}{$g2})) { next; }
	if(exists($hashp{$g2}{$g1})) { next; }
	print OO "$line\n";
}
close OO;

print "Generating initial database.tsv file ...\n";
$eachside = 20;
$pairs = "$tmpdir/out.pairs.nopara";
####################################################
$in3 = "gn.ftmp";
system "cut -f3-4 $genes_transcripts_exons | sort -u | grep -v Gene  > $tmpdir/$in3";
%hash_gid_to_gname = ();
open(FF3, "<$tmpdir/$in3");
@filethree = <FF3>; close(FF3);
system "rm $tmpdir/$in3";
for($i=0; $i<@filethree; $i++)
{
	$line3 = $filethree[$i]; chomp($line3);
	($gid3, $gname3) = split(/\t/, $line3);
	$hash_gid_to_gname{$gid3} = $gname3;
}

$in2 = "ftmp";
if(!(-e "$in2")) {
system "cut -f1-3,10- $genes_transcripts_exons | sort -u | sort -nk3 | grep -v 'Exon' > $tmpdir/$in2";
}

open(FF, "<$tmpdir/$in2"); @file = <FF>; close(FF);
system "rm $tmpdir/$in2";
%hash = (); %hash_estart = (); %hash_eend = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$chr = $a[0]; $strand = $a[1]; $gid = $a[2]; $eid = $a[3]; $estart = $a[4]; $eend = $a[5]; $strand =~ s/-1/-/g; $strand =~ s/1/+/g;
	if(!(exists($hash{$gid}))) { $ii = 1; }
	$hash{$gid}{$ii} = "$chr,$strand,$ii,$eid,$estart,$eend";
	$hash_estart{$gid}{$estart} = $eend; $hash_eend{$gid}{$eend} = $estart;
	$ii++;
}

open(FF, "<$pairs");
@file = <FF>; close(FF);
%main_hash1 = ();
%main_hash2 = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$g1 = $a[0]; $g2 = $a[1];
	if( !(exists($hash{$g1})) or !(exists($hash{$g2})) ) { print "NOT\t$line\n"; next; }
	%hash1 = ();
	foreach $e1 (sort {$a <=> $b} keys(%{$hash{$g1}}))
	{
		($chr,$strand,$num,$eid,$estart,$eend) = split(/,/, $hash{$g1}{$e1});
		$len_exon = $eend-$estart+1;
		if($strand eq "+") {
			$upstart = $eend-$eachside+1; $upend = $eend;
		}
		if($strand eq "-") {
			$upstart = $estart; $upend = $upstart+$eachside-1;
		}
		$hkey1 = "$chr#$upstart#$upend#$strand#$g1";
		$hash1{$hkey1}{$len_exon} = $eid;
	}
	%hash11 = ();
	$totnum = 0;
	foreach $kk1 (keys(%hash1))
	{
		($chr, $upstart, $upend, $strand, $g1) = split(/#/, $kk1);
		foreach $len1 (sort {$b <=> $a} keys(%{$hash1{$kk1}}))
		{
			$eid = $hash1{$kk1}{$len1};
			$len11 = $len1;
			last;
		}
		$hash11{$upstart} = "$upend#$chr#$strand#$g1#$eid#$len11";
		$totnum++;
	}
	$count_exon = 0;
	foreach $kkk1 (sort {$a <=> $b} keys(%hash11))
	{
		$count_exon++;
		$upstart = $kkk1;
		($upend, $chr, $strand, $g1, $eid, $len_exon) = split(/#/, $hash11{$kkk1});
		$unique_key = "$g1#e$count_exon"; $main_hash1{$g1}{$unique_key} = "$chr\t$upstart\t$upend\t$eid#len$len_exon#$totnum\t.\t$strand";
	}

	%hash2 = ();
	foreach $e2 (sort {$a <=> $b} keys(%{$hash{$g2}})) { $tnum = $e2; }
	foreach $e2 (sort {$a <=> $b} keys(%{$hash{$g2}}))
	{
		($chr,$strand,$num,$eid,$estart,$eend) = split(/,/, $hash{$g2}{$e2});
		$len_exon = $eend-$estart+1;
		if($strand eq "+") {
			$upstart = $estart; $upend = $upstart+$eachside-1;
		}
		if($strand eq "-") {
			$upstart = $eend-$eachside+1; $upend = $eend;
		}
		$hkey2 = "$chr#$upstart#$upend#$strand#$g2";
		$hash2{$hkey2}{$len_exon} = $eid;
	}
	%hash22 = ();
	$totnum = 0;
	foreach $kk2 (keys(%hash2))
	{
		($chr, $upstart, $upend, $strand, $g2) = split(/#/, $kk2);
		foreach $len2 (sort {$b <=> $a} keys(%{$hash2{$kk2}}))
		{
			$eid = $hash2{$kk2}{$len2};
			$len22 = $len2;
			last;
		}
		$hash22{$upstart} = "$upend#$chr#$strand#$g2#$eid#$len22";
		$totnum++;
	}
	$count_exon = 0;
	foreach $kkk2 (sort {$a <=> $b} keys(%hash22))
	{
		$count_exon++;
		$upstart = $kkk2;
		($upend, $chr, $strand, $g2, $eid, $len_exon) = split(/#/, $hash22{$kkk2});
		$unique_key = "$g2#e$count_exon"; $main_hash2{$g2}{$unique_key} = "$chr\t$upstart\t$upend\t$eid#len$len_exon#$totnum\t.\t$strand";
	}
}

$just_count = 0;
open(FOUT, ">$tmpdir/tmp.bed");
open(FF, "<$pairs");
@file = <FF>; close(FF);
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$g1 = $a[0]; $g2 = $a[1];
	foreach $k1 (sort keys(%{$main_hash1{$g1}}))
	{
		$v1 = $main_hash1{$g1}{$k1};
		($chr1,$start1,$end1,$name1,$col5_1,$strand1) = split(/\t/, $v1);
		$start1_cp = $start1;
		$start1--; # make it 0-based
		$v1 = "$chr1\t$start1\t$end1\t$name1\t$col5_1\t$strand1";
		foreach $k2 (sort keys(%{$main_hash2{$g2}}))
		{
			$v2 = $main_hash2{$g2}{$k2};
			($chr2,$start2,$end2,$name2,$col5_2,$strand2) = split(/\t/, $v2);
			$start2_cp = $start2;
			$start2--; # make it 0-based
			$v2 = "$chr2\t$start2\t$end2\t$name2\t$col5_2\t$strand2";
			$distance_between_pairs = -1000;
			if($strand1 eq '+')
			{
				$distance_between_pairs = $start2_cp - $end1 -1; # because of 1-based
			}
			else
			{
				$distance_between_pairs = $start1_cp - $end2 -1; # because of 1-based
			}
			if($distance_between_pairs > $distance_cut_off) { next; }
			if($distance_between_pairs < 0) { next; }
			print FOUT "$k1#$k2\t$v1\t$v2\n";
			$just_count++;
		}
	}
}
close FOUT;

print "Implementing Same isoform filter ...\n";
$inf = "$tmpdir/tmp.bed"; open(FF, "<$inf"); @file = <FF>; close(FF);
open(FLAGTMP, ">$tmpdir/tmpfl.bed");
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	($id1, $en1, $id2, $en2) = split(/#/, $a[0]);
	($eid1, $col2, $col3) = split(/#/, $a[4]);
	($eid2, $col4, $col5) = split(/#/, $a[10]);
	$start1 = $a[2]; $end1 = $a[3];
	$start2 = $a[8]; $end2 = $a[9];
	$strand = $a[6];
	# breakpoint coordinates are 1-based; final.bed coordinates are 0-based;
	if($strand eq "+") { $break1 = $end1; $break2 = $start2+1; }
	if($strand eq "-") { $break1 = $start1+1; $break2 = $end2; }

	$flag_to_discard = 0; 	# 0 means don't discard
	if($strand eq "+")
	{
		if(exists($hash_eend{$id1}{$break1}) and exists($hash_estart{$id1}{$break2})) { $flag_to_discard = 1; }
		if(exists($hash_eend{$id2}{$break1}) and exists($hash_estart{$id2}{$break2})) { $flag_to_discard = 1; }
	}
	if($strand eq "-")
	{
		if(exists($hash_eend{$id1}{$break2}) and exists($hash_estart{$id1}{$break1})) { $flag_to_discard = 1; }
		if(exists($hash_eend{$id2}{$break2}) and exists($hash_estart{$id2}{$break1})) { $flag_to_discard = 1; }
			
	}
	if($flag_to_discard == 0)
	{
		print FLAGTMP "$line\n";
	}
}
close FLAGTMP;
system "mv $tmpdir/tmpfl.bed $tmpdir/tmp.bed";

$inf = "$tmpdir/tmp.bed";
open(FF, "<$inf");
@file = <FF>; close(FF);
%hash = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$k = "$a[1]#$a[2]#$a[3]#$a[6]#$a[7]#$a[8]#$a[9]#$a[12]";
	$id = $a[0];
	push(@{$hash{$k}}, $id);
}

%id_to_remove_from_final = ();
open(FEF, ">$tmpdir/efinal.tsv");
foreach $hk (sort keys %hash)
{
	$num = @{$hash{$hk}};
	if($num < 2) { next; }
	$get_ids = "";
	for($ia=0; $ia<@{$hash{$hk}}; $ia++)
	{
		$get_ids .= "${$hash{$hk}}[$ia],";
		if($ia>0) { $id_to_remove_from_final{${$hash{$hk}}[$ia]} = 1; }
	}
	$get_ids =~ s/,$//g;
	print FEF "$hk\t$get_ids\t$num\n";
}
close FEF;

open(NFINAL, ">$tmpdir/ntmp.bed");
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$id = $a[0];
	if(exists($id_to_remove_from_final{$id})) { next; }
	print NFINAL "$line\n";
}
close NFINAL;

system "mv $tmpdir/ntmp.bed $tmpdir/tmp.bed";
system "rm $tmpdir/efinal.tsv";

system "cut -f2-7 $tmpdir/tmp.bed > $tmpdir/tmp1.bed";
system "cut -f8-13 $tmpdir/tmp.bed > $tmpdir/tmp2.bed";
system "$bedtools getfasta -name -fi $genomefa -bed $tmpdir/tmp1.bed -s > $tmpdir/tmp1.fasta";
system "$bedtools getfasta -name -fi $genomefa -bed $tmpdir/tmp2.bed -s > $tmpdir/tmp2.fasta";
system "grep -v '^>' $tmpdir/tmp1.fasta > $tmpdir/s1";
system "grep -v '^>' $tmpdir/tmp2.fasta > $tmpdir/s2";
system "paste -d '' $tmpdir/s1 $tmpdir/s2 > $tmpdir/s1s2";
system "paste $tmpdir/tmp.bed $tmpdir/s1s2 > $tmpdir/final.bed";
system "rm $tmpdir/tmp.bed $tmpdir/tmp1.bed $tmpdir/tmp2.bed $tmpdir/tmp1.fasta $tmpdir/tmp2.fasta $tmpdir/s1 $tmpdir/s2 $tmpdir/s1s2";

print "Implementing Junction sequence length filter ...\n";
open(FF, "<$tmpdir/final.bed");@finalfile = <FF>; close(FF);
open(OFFO, ">$tmpdir/fEfinal.bed");
for($i=0; $i<@finalfile; $i++)
{
	$line = $finalfile[$i]; chomp($line);
	@a = split(/\t/, $line);
	$lid = $a[0]; $ljseq = $a[13];
	$name1 = $a[4]; $name2 = $a[10];
	($lexid1, $lenexon1) = split(/#/, $name1); $lenexon1 =~ s/^len//g;
	($lexid2, $lenexon2) = split(/#/, $name2); $lenexon2 =~ s/^len//g;
	if($lenexon1 < $eachside or $lenexon2 < $eachside) { next; }

	# adding gene names
	($llgid1, $lle1, $llgid2, $lle2) = split(/#/, $lid);
	$llgn1 = $hash_gid_to_gname{$llgid1}; $llgn2 = $hash_gid_to_gname{$llgid2};
	$gnamepair = $llgn1."#".$llgn2;
	print OFFO "$line\t$gnamepair\n";
}
close OFFO;
system "mv $tmpdir/fEfinal.bed $tmpdir/final.bed";

print "Generating exons.tsv file ...\n";
system "cut -f1,11,12 $genes_transcripts_exons |grep -v 'Chromosome' > $tmpdir/123; cut -f10 $genes_transcripts_exons |grep -v 'Exon' > $tmpdir/4; cut -f2 $genes_transcripts_exons |grep -v 'Strand' | perl -p -e 's/-1/-/g' | perl -p -e 's/1/+/g' > $tmpdir/6;";
$wcl = `cat $tmpdir/123 | wc -l `; chomp($wcl);
system "echo \$(seq -s \",\" 1 $wcl)|perl -p -e 's/,/\\n/g' |perl -p -e 's/^.+/./g' > $tmpdir/5";
system "paste $tmpdir/123 $tmpdir/4 $tmpdir/5 $tmpdir/6 | sort -u | sort -k1,1 -k2,2n > $tmpdir/ife";
system "rm $tmpdir/123 $tmpdir/4 $tmpdir/5 $tmpdir/6";
open(FF, "<$tmpdir/ife");
@file = <FF>; close(FF);
system "rm $tmpdir/ife";
open(OO, ">$tmpdir/ife.bed");
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$start = $a[1];
	$start--;
	print OO "$a[0]\t$start\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t\n";
}
close OO;

system "$bedtools getfasta -name -fi $genomefa -bed $tmpdir/ife.bed -s > $tmpdir/exons.fa";
system "rm $tmpdir/ife.bed";
system "perl -pi -e 's/\\(.\\)//g' $tmpdir/exons.fa";
$f2t = fasta2tsv("$tmpdir/exons.fa");
system "mv $tmpdir/exons.fa.tsv $tmpdir/exons.tsv";

print "Implementing Exon filter (This step uses blast against all exons and therefore is a time consuming step).  ...\n";
$in = "$tmpdir/final.bed";
open(FF, "<$in");
@file = <FF>; close(FF);
%hash = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$id = $a[0]; $jseq = $a[13];
	$rjseq = $jseq; $rjseq =~ tr/ATGC/TACG/; $rjseq = reverse $rjseq;
	$hash{$jseq}{$id} = 1;
	$hash{$rjseq}{$id} = 1;
}

$counter = 0;
open(OO, ">$tmpdir/jseqs.fasta");
foreach $kseq (sort keys(%hash))
{
	$counter++;
	$allids = "";
	foreach $kid (sort keys(%{$hash{$kseq}}))
	{
		$allids .= "$kid;";
	}
	$allids =~ s/;$//g;
	print OO ">jseq$counter:$allids\n$kseq\n";
}
close OO;

$in = "$tmpdir/exons.tsv";
open(FF, "<$in");
@file = <FF>; close(FF);
open(OO, ">$tmpdir/ids.exon.length");
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$id = $a[0]; $seq = $a[1];
	$len = length($seq);
	print OO "$id\t$len\n";
}
close OO;
system "cp $tmpdir/ids.exon.length $tmpdir/1.ids.exon.length";

system "$makeblastdb -dbtype nucl -in $tmpdir/exons.fa";
system "$blastn -query $tmpdir/jseqs.fasta -num_threads $cpu -db $tmpdir/exons.fa -outfmt 7 -strand plus -task blastn > $tmpdir/out";

$eachside = 20;
$inb = "$tmpdir/out";
$inl = "$tmpdir/ids.exon.length";
$inf = "$tmpdir/jseqs.fasta";
$cut_off_pid = 90;
$cut_off_aper = 90;
########################################
open(FF, "<$inl"); @file = <FF>; close(FF);
%hashl = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$id1 = $a[0]; $leng = $a[1];
	$hashl{$id1} = $leng;
}

# blast file
open(FF, "<$inb");
%hashb = ();
while(<FF>)
{
	$line = $_; chomp($line);
	if($line =~ /^#/) { next; }
	@a = split(/\t/, $line);
	$query = $a[0]; $target = $a[1]; $pid = $a[2]; $alen = $a[3];
	$qs = $a[6]; $qe = $a[7]; $ts = $a[8]; $te = $a[9];  # q=query, t=target
	$len_q = $eachside*2;
	$len_t = $hashl{$target};
	$minlen = $len_q; if($len_t < $len_q) { $minlen = $len_t; }
	if($minlen < $len_q) { next; }
	$aper = ($alen*100)/$minlen;
	$hashb{$query}{$target} = "$pid,$aper";
}
close FF;

# jseqs.fasta file
open(FF, "<$inf"); @file = <FF>; close(FF);
open(OTM, ">$tmpdir/tmp9.9");
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	if($line !~ /^>/) { next; }
	$id = $line; $id =~ s/^>//g;
	$nid = $id; $nid =~ s/^jseq.+://g;
	$flag = 0;	# if 1 it means to discard
	foreach $ktarget (sort keys (%{$hashb{$id}}))
	{
		($pidd, $aperr) = split(/,/, $hashb{$id}{$ktarget});
		if($pidd >= $cut_off_pid and $aperr >= $cut_off_aper)
		{
			$flag = 1;
			last;
		}
	}
	if($flag == 1)
	{
		@anid = split(/;/, $nid);
		for($ianid=0; $ianid<@anid; $ianid++)
		{
			print OTM "$anid[$ianid]\n";
		}
	}
}
close OTM;

$allfg = "$tmpdir/tmp9.9,1,$tmpdir/final.bed,1,$tmpdir/filt_final.bed,1";
$runfg = fastgrep($allfg);

print "Implementing Similarity filter. (This step uses blast and therefore is a time consuming step). ...\n";
system "cut -f1 $tmpdir/out.pairs.nopara > $tmpdir/up.tmp";
system "cut -f2 $tmpdir/out.pairs.nopara > $tmpdir/down.tmp";
system "cut -f$cutcols $genes_transcripts_exons | grep -v '^Chromosome' |grep -v 'Strand' | sort -u > $tmpdir/ex.tmp";
open(FF, "<$tmpdir/ex.tmp"); @exfile = <FF>; close(FF);
open(OFF, ">$tmpdir/ex.tmp.bed");
for($i=0; $i<@exfile; $i++)
{
	$line = $exfile[$i]; chomp($line);
	@a = split(/\t/, $line);
	$exchr = $a[0]; $exstrand = $a[1]; $exstart = $a[4]; $exend = $a[5];
	$exgid = $a[2]; $exeid = $a[3];
	$bid = $exgid."#$exeid";
	if($exstrand eq "1") { $exstrand = "+"; }
	if($exstrand eq "-1") { $exstrand = "-"; }
	$exstart--;	# to make compatible with bed file
	print OFF "$exgid\t$exchr\t$exstart\t$exend\t$bid\t.\t$exstrand\n";
}
close OFF;
$fg1i = "$tmpdir/up.tmp,0,$tmpdir/ex.tmp.bed,1,$tmpdir/fup"; $fg1 = fastgrep($fg1i);
system "cat $tmpdir/fup |cut -f2-7 |sort -u  |sort -k1,1 -k2,2n > $tmpdir/up.bed";
$fg2i = "$tmpdir/down.tmp,0,$tmpdir/ex.tmp.bed,1,$tmpdir/fdown"; $fg2 = fastgrep($fg2i);
system "cat $tmpdir/fdown |cut -f2-7 |sort -u  |sort -k1,1 -k2,2n > $tmpdir/down.bed";
system "$bedtools getfasta -name -fi $genomefa -bed $tmpdir/up.bed -s > $tmpdir/eu.fasta";
system "$bedtools getfasta -name -fi $genomefa -bed $tmpdir/down.bed -s > $tmpdir/ed.fasta";
open(FF, "<$tmpdir/ex.tmp.bed");
@file = <FF>; close(FF);
open(OOL, ">$tmpdir/ids.exon.length");
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$len = $a[3] - $a[2];
	$id = $a[4];
	print OOL "$id\t$len\n";
}
close OOL;
system "cp $tmpdir/ids.exon.length $tmpdir/2.ids.exon.length";
system "$makeblastdb -dbtype nucl -in $tmpdir/ed.fasta";
system "$blastn -query $tmpdir/eu.fasta -num_threads $cpu -db $tmpdir/ed.fasta -outfmt 7 -strand plus -task blastn > $tmpdir/out";

$inb = "$tmpdir/out";
$inl = "$tmpdir/ids.exon.length";
$inf = "$tmpdir/filt_final.bed";
$gte = $genes_transcripts_exons;
$cut_off_pid = 70;
$cut_off_aper = 70;
########################################
open(FF, "<$gte");
@file = <FF>; close(FF);
%hashmul = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$chr = $a[0]; $strand = $a[1]; $ensid = $a[2]; $gname = $a[3];
	if($strand eq "1") { $strand = "+"; }
	if($strand eq "-1") { $strand = "-"; }
	$hashmul{$gname}{$chr}{$strand}{$ensid}++;
}

# length file
open(FF, "<$inl"); @file = <FF>; close(FF);
%hashl = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$id1 = $a[0]; $leng = $a[1];
	@aa = split(/#/, $id1);
	$gid = $aa[0]; $eid = $aa[1];
	$hashl{$gid}{$eid} = $leng;
}

# blast file
open(FF, "<$inb");
%hashb = ();
while(<FF>)
{
	$line = $_; chomp($line);
	if($line =~ /^#/) { next; }
	@a = split(/\t/, $line);
	$query = $a[0]; $target = $a[1]; $pid = $a[2]; $alen = $a[3];
	$qs = $a[6]; $qe = $a[7]; $ts = $a[8]; $te = $a[9];  # q=query, t=target
	$query =~ s/\(.\)//g; $target =~ s/\(.\)//g;
	($qgid, $qeid)  = split(/#/, $query);
	($tgid, $teid)  = split(/#/, $target);
	$len_q = $hashl{$qgid}{$qeid};
	$len_t = $hashl{$tgid}{$teid};
	$minlen = $len_q; if($len_t < $len_q) { $minlen = $len_t; }
	$aper = ($alen*100)/$minlen;
	$mid1 = $qgid."#$qeid";
	$mid2 = $tgid."#$teid";
	$hashb{$mid1}{$mid2} = "$pid,$aper";
	$hashb{$mid2}{$mid1} = "$pid,$aper";
}
close FF;

# final.bed file
open(OOLT, ">$tmpdir/tmp7.7");
open(FF, "<$inf"); @file = <FF>; close(FF);
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$id = $a[0]; $gname12 = $a[14];
	($gname1, $gname2) = split(/#/, $gname12);
	$chrom = $a[1]; $strand  = $a[6];
	($eid1, $le1) = split(/#/, $a[4]);
	($eid2, $le2) = split(/#/, $a[10]);
	($gid1, $var1, $gid2, $var2) = split(/#/, $id);
	$mid1 = $gid1."#$eid1"; $mid2 = $gid2."#$eid2";
	$filter_out = 0; # 0 means don't filter out this entry: 1 means filter out
 	foreach $mul2 (sort keys (%{$hashmul{$gname2}{$chrom}{$strand}})) 
	{
		$each_gid2 = $mul2; $gid2 = $each_gid2;
		foreach $k1 (sort keys (%{$hashl{$gid2}})) 
		{
			$make_mid2 = "$gid2#$k1";
			($pid, $aper) = split(/,/, $hashb{$mid1}{$make_mid2});
			if($pid >= $cut_off_pid and $aper >= $cut_off_aper) { $filter_out = 1; last; }
		}
	}
 	foreach $mul1 (sort keys (%{$hashmul{$gname1}{$chrom}{$strand}})) 
	{
		$each_gid1 = $mul1; $gid1 = $each_gid1;
		foreach $k1 (sort keys (%{$hashl{$gid1}})) 
		{
			$make_mid1 = "$gid1#$k1";
			($pid, $aper) = split(/,/, $hashb{$mid2}{$make_mid1});
			if($pid >= $cut_off_pid and $aper >= $cut_off_aper) { $filter_out = 1; last; } 
		}
	}
	if($filter_out == 1) # 0 means retain it, 1 means filter it out
	{
		#print OOLT "$line\t$pid\t$aper\n";
		print OOLT "$line\n";
	}
}
close OOLT;
print "Generating final database.tsv file\n";
$fg3i = "$tmpdir/tmp7.7,1,$tmpdir/filt_final.bed,1,$tmpdir/filtered_final.bed,1"; $fg3 = fastgrep($fg3i);
system "mv $tmpdir/filtered_final.bed $tmpdir/database.tsv";

if($optsnp == 1)
{
	print "Mapping SNPs on junction sequences\n";
	system "cut -f2-7 $tmpdir/database.tsv | sort -u | sort -k1,1 -k2,2n > $tmpdir/a.bed";
	system "cut -f8-13 $tmpdir/database.tsv | sort -u | sort -k1,1 -k2,2n > $tmpdir/b.bed";
	system "$bedtools intersect -wa -a $dbsnpfile -b $tmpdir/a.bed -s > $tmpdir/dbsnp.a.bed";
	system "$bedtools intersect -wa -a $dbsnpfile -b $tmpdir/b.bed -s > $tmpdir/dbsnp.b.bed";
	system "sort -u $tmpdir/dbsnp.a.bed | sort -k1,1 -k2,2n > $tmpdir/u.bed; mv $tmpdir/u.bed $tmpdir/dbsnp.a.bed";
	system "sort -u $tmpdir/dbsnp.b.bed | sort -k1,1 -k2,2n > $tmpdir/u.bed; mv $tmpdir/u.bed $tmpdir/dbsnp.b.bed";
	system "cat $tmpdir/dbsnp.a.bed $tmpdir/dbsnp.b.bed |sort -u |sort -k1,1 -k2,2n > $tmpdir/final.dbsnp.bed";

	$in = "$tmpdir/final.dbsnp.bed";
	$in2 = $allvcffile;
	open(FF, "<$in");
	@file = <FF>; close(FF);
	%hash = ();
	for($i=0; $i<@file; $i++)
	{
		$line = $file[$i]; chomp($line);
		@a = split(/\t/, $line);
		$id = $a[3];
		$hash{$id} = 1;
	}

	open(OO, ">$tmpdir/final.dbsnp.vcf");
	open(FF, "<$in2");
	while(<FF>)
	{
		$line = $_; chomp($line);
		if($line =~ /^#/) { next; }
		@a = split(/\t/, $line);
		$pos = $a[1]; $id = $a[2];
		$ref = $a[3]; $alt = $a[4];
		if(exists($hash{$id})) {
			print OO "$line\n";
		}
	}
	close FF;
	close OO;
	system "grep ';CAF=' $tmpdir/final.dbsnp.vcf > $tmpdir/freq.final.dbsnp.vcf";

	$in = "$tmpdir/freq.final.dbsnp.vcf";
	$caf_cut_off = 0.1;
	##################################################
	open(FF, "<$in");
	@file = <FF>; close(FF);
	open(VOO, ">$tmpdir/caf$caf_cut_off.final.vcf");
	for($i=0; $i<@file; $i++)
	{
		$line = $file[$i]; chomp($line);
		if($line  !~ /CAF=/) { next; }
		@a = split(/\t/, $line);
		$id = $a[2]; $info = $a[7];
		$rv = 0; if($info =~ /;RV;/) { $rv = 1; }
		$info =~ s/^.+;CAF=//g;
		$info =~ s/;.+$//g;
		@ai = split(/,/, $info);
		$flag = 1; # 0 means don't discard and 1 means discard
		for($j=0; $j<@ai; $j++)
		{
			if($j == 0) { next; } # ref
			if($ai[$j] >= $caf_cut_off) { $flag = 0; last; } # if any alternate allele is >= 0.1 then don't discard
		}
		if($flag == 0)
		{
			print VOO "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$info\tRV$rv\n";
		}
	}
	close VOO;
	system "cut -f3 $tmpdir/caf$caf_cut_off.final.vcf > $tmpdir/tmp.caf";
	$fgvi = "$tmpdir/tmp.caf,0,$tmpdir/final.dbsnp.bed,4,$tmpdir/ftmp.bed";
	$fgv = fastgrep($fgvi);
	system "sort -u $tmpdir/ftmp.bed | sort -k1,1 -k2,2n  > $tmpdir/final2.dbsnp.bed;"; 

	$final_file = "$tmpdir/database.tsv";
	$caf_file = "$tmpdir/caf$caf_cut_off.final.vcf";
	$bed_file = "$tmpdir/final2.dbsnp.bed";
	print "Adding junction sequence variants in the database\n";
	##################################################
	# save caf file in hash
	open(FF, "<$caf_file");
	@file = <FF>; close(FF);
	%hashcaf = ();
	%hashcafrv = ();
	for($i=0; $i<@file; $i++)
	{
		$line = $file[$i]; chomp($line);
		@a = split(/\t/, $line);
		$id = $a[2]; $rv = $a[6];
		$hashcaf{$id} = $line;
		if($rv eq "RV1") { $hashcafrv{$id}++; }
	}
	
	# save bed file in hash
	open(FF, "<$bed_file");
	@file = <FF>; close(FF);
	%hashbed = ();
	for($i=0; $i<@file; $i++)
	{
		$line = $file[$i]; chomp($line);
		@a = split(/\t/, $line);
		$chr = $a[0]; $start = $a[1]; $end = $a[2]; $id = $a[3]; $strand = $a[5];
		$hashbed{$chr}{$strand}{$start}{$end} = $id;
	}
	
	# now final.bed file
	open(FF, "<$final_file");
	open(VOOO, ">$tmpdir/vdatabase.tsv");
	while(<FF>)
	{
		$line = $_; chomp($line);
		@a = split(/\t/, $line);
		$id = $a[0]; $c1 = $a[1]; $s1 = $a[2]; $e1 = $a[3]; $str1 = $a[6]; $c2 = $a[7]; $s2 = $a[8]; $e2 = $a[9]; $str2 = $str1; $jseq = $a[13];
		$each_side = $e1-$s1;
		@ajseq = split(//, $jseq);
		$allnjseq1 = "";
		foreach $ks (sort keys(%{$hashbed{$c1}{$str1}}))
		{
			foreach $ke (sort keys(%{$hashbed{$c1}{$str1}{$ks}}))
			{
				if($ks >= $s1 and $ke <= $e1)
				{
					$diff = $ks-$s1;
					$diff1 = ($e1-$s1)-$diff-1; # if RV1 i.e. -ve
					$thisid1 = $hashbed{$c1}{$str1}{$ks}{$ke}; # $str1 is same as $str2
					if(!(exists($hashcaf{$thisid1}))) { next; }
					@acaf = split(/\t/, $hashcaf{$thisid1});
					$cafpos = $acaf[1]; # this is 1-based
					$cafref = $acaf[3]; $cafalt = $acaf[4]; $caffreq = $acaf[5];
					@aalt = split(/,/, $cafalt);
					@afreq = split(/,/, $caffreq); shift(@afreq); # remove reference frequency
					$reflen = length($cafref); if($reflen != 1) { next; }
					for($alt=0; $alt<@aalt; $alt++)
					{
						$alternate_allele = $aalt[$alt]; # reference allele is $cafref
						if($alternate_allele eq ".") { next; }
						$len_alternate_allele = length($alternate_allele);
						if($len_alternate_allele != 1) { next; }
						if($afreq[$alt] < $caf_cut_off) { next; }
						@cpajseq = @ajseq;
						if(exists($hashcafrv{$thisid1})) { $nt_from_jseq = $cpajseq[$diff1]; $nt_from_jseq =~ tr/ATGC/TACG/; $alternate_allele =~ tr/ATGC/TACG/; $cpajseq[$diff1] = $alternate_allele; $newjseq = join("", @cpajseq); }
						else { $nt_from_jseq = $cpajseq[$diff]; $cpajseq[$diff] = $alternate_allele; $newjseq = join("", @cpajseq); }
						$allnjseq1 .= "$newjseq#";
						#print "a$diff\t$id\t$thisid1\t$nt_from_jseq\t$cafref\n";
						if($nt_from_jseq ne $cafref) { 
							#print "a$diff\t$id\t$thisid1\t$nt_from_jseq\t$cafref\n";
							#print "\nExiting!\n"; #exit;
						}
					}
	
				}
			}
		}
		$allnjseq1 =~ s/#$//g;
		if($allnjseq1 eq "") { $allnjseq1 = "-NA-"; }
	
		$get_left_half = "";
		for($sf=1; $sf<=$each_side; $sf++) { $hseq = shift(@ajseq); $get_left_half .= $hseq; }
		$allnjseq2 = "";
		foreach $ks (sort keys(%{$hashbed{$c2}{$str1}}))
		{
			foreach $ke (sort keys(%{$hashbed{$c2}{$str1}{$ks}}))
			{
				if($ks >= $s2 and $ke <= $e2)
				{
					$diff = $ks-$s2;
					$diff1 = ($e2-$s2)-$diff-1; # if RV1 i.e. -ve
					$thisid1 = $hashbed{$c2}{$str1}{$ks}{$ke}; # $str1 is same as $str2
					if(!(exists($hashcaf{$thisid1}))) { next; }
					@acaf = split(/\t/, $hashcaf{$thisid1});
					$cafpos = $acaf[1]; # this is 1-based
					$cafref = $acaf[3]; $cafalt = $acaf[4]; $caffreq = $acaf[5];
					@aalt = split(/,/, $cafalt);
					@afreq = split(/,/, $caffreq); shift(@afreq); # remove reference frequency
					$reflen = length($cafref); if($reflen != 1) { next; }
					for($alt=0; $alt<@aalt; $alt++)
					{
						$alternate_allele = $aalt[$alt]; # reference allele is $cafref
						if($alternate_allele eq ".") { next; }
						$len_alternate_allele = length($alternate_allele);
						if($len_alternate_allele != 1) { next; }
						if($afreq[$alt] < $caf_cut_off) { next; }
						@cpajseq = @ajseq;
						if(exists($hashcafrv{$thisid1})) { $nt_from_jseq = $cpajseq[$diff1]; $nt_from_jseq =~ tr/ATGC/TACG/; $alternate_allele =~ tr/ATGC/TACG/; $cpajseq[$diff1] = $alternate_allele; $newjseq = join("", @cpajseq); }
						else { $nt_from_jseq = $cpajseq[$diff]; $cpajseq[$diff] = $alternate_allele; $newjseq = join("", @cpajseq); }
						$allnjseq2 .= $get_left_half."$newjseq#";
						#print "b$diff\t$id\t$thisid1\t$nt_from_jseq\t$cafref\n";
						if($nt_from_jseq ne $cafref) { 
							#print "b$diff\t$id\t$thisid1\t$nt_from_jseq\t$cafref\n";
							#print "\nExiting!\n"; #exit;
						}
					}
				}
			}
		}
		$allnjseq2 =~ s/#$//g;
		if($allnjseq2 eq "") { $allnjseq2 = "-NA-"; }
	
		print VOOO "$line\t$allnjseq1\t$allnjseq2\n";
	}
	close FF;
	close VOOO;
	system "mv $tmpdir/vdatabase.tsv $tmpdir/database.tsv";
} # end of if statement --> if($optsnp == 1)

opendir(DIR,"$tmpdir");
while($tfile=readdir(DIR))
{
	chomp($tfile);
	if($tfile eq "." or $tfile eq "..") { next; }
	if($tfile eq "database.tsv" or $tfile eq "exons.tsv") { next; }
	system "rm $tmpdir/$tfile";
}
close DIR;
$pwd = `pwd`; chomp($pwd);
print "\tRTCpredictor database files \"database.tsv\" and \"exons.tsv\" are present in $pwd/$tmpdir/\n";
print "====== Done ======\n";
#################### subroutines ####################
sub fasta2tsv
{
my $inp = shift;
my $out = $inp.".tsv";
###############################
open(FF, "<$inp");
@file = <FF>; close FF;
open(OO1, ">$out");
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	if($i==0 and $line !~ /^>/) { print "\tNot a valid fasta file!\n"; exit; }
	if($line =~ /^>/)
	{
		$cpline = $line;
		$hdr = ""; $hdr = $cpline; $hdr =~ s/^>//g;
		$catseq = "";
		for($j=$i+1; $j<@file; $j++)
		{
			$jline = $file[$j]; chomp($jline);
			if($jline =~ /^>/) { $i = $j - 1; last; }
			$catseq .= $jline;
		}
		print OO1 "$hdr\t$catseq\n";
	}
}
close OO1;
return 1;
}

sub fastgrep
{
my $data = shift;
($fl1, $fl1c, $fl2, $fl2c, $fl3, $ooptv) = split(/,/, $data);
$inhash = $fl2;
$inhash_n = $fl2c - 1;
$infile = $fl1;
$infile_n = $fl1c - 1;
$outfile = $fl3;
$minusv = $ooptv;
#$minusi = $ARGV[6];
if($inhash_n == -1) { $save_col_or_wholefile = 0; } else { $save_col_or_wholefile = 1; }
if($infile_n == -1) { $search_col_or_wholefile = 0; } else { $search_col_or_wholefile = 1; }
$deli1 = "\t";
$deli2 = "\t";
######################################
open(T6, "<$inhash");
@ft6 = <T6>; close T6;
open(T11, "<$infile");
@ft11 = <T11>; close T11;
open(TMP12, ">$outfile");
#print TMP12 "started\n";

%ht6 = ();
for($t6=0; $t6<@ft6; $t6++)
{
        $lt6 = $ft6[$t6]; chomp($lt6);
        @alt6 = split(/$deli2/, $lt6);
        $fusn = $alt6[$inhash_n];
	if($save_col_or_wholefile == 0) { $fusn = $lt6; }
	if($minusi == 1) { $fusn = lc($fusn); }
        push(@{$ht6{$fusn}},$lt6);
}

%present_grepwf = ();
for($t11=0; $t11<@ft11; $t11++)
{
        $lt11 = $ft11[$t11]; chomp($lt11);
        @alt11 = split(/$deli1/, $lt11);
        $colsearch = $alt11[$infile_n];
	if($search_col_or_wholefile == 0) { $colsearch = $lt11; }
	if($minusi == 1) { $colsearch = lc($colsearch); }
        if(exists($ht6{$colsearch}))
        {
                for($phfe=0; $phfe<@{$ht6{$colsearch}}; $phfe++)
                {
                        if($minusv == 0) { print TMP12 "${$ht6{$colsearch}}[$phfe]\n"; }
			$present_grepwf{$colsearch}++;
                }
        }
}

for($t6=0; $t6<@ft6; $t6++)
{
	$lt6 = $ft6[$t6]; chomp($lt6);
	@alt6 = split(/$deli2/, $lt6);
	$fusn = $alt6[$inhash_n];
	if($save_col_or_wholefile == 0) { $fusn = $lt6; }
	if(!(exists($present_grepwf{$fusn})))
	{
		if($minusv == 1) { print TMP12 "$lt6\n"; }
	}
}
close TMP12;
return 1;
}
#################### end of subroutines ####################
