BEGIN { $program_name="RTpred"; $home=$ENV{'HOME'}; $hcatfile="$home/.$program_name/.$program_name.rc"; open(FF, "$hcatfile"); $varlib = <FF>; chomp($varlib); close FF; if($varlib eq "") { $varlib=$home; } } use lib $varlib;
$whichrg = `which rg`; chomp($whichrg); if(!(-e $whichrg)) { print "\n\n\tError! ripgrep not found. Please install ripgrep before running. README.txt file has instructions on how to install ripgrep \n"; exit; }
$usage = <<EOF;
	$program_name: A software tool for the prediction of read-through chimeric RNAs from RNA-Seq data.

	Command to run on paired-end RNAseq fastq files
	perl $program_name.pl -d database.tsv -e exons.tsv -t paired -i read_1.fastq.gz -j read_2.fastq.gz [Optional parameters]

	Command to run on single-end RNAseq fastq files
	perl $program_name.pl -d database.tsv -e exons.tsv -t single -k read.fastq.gz [Optional parameters]

	Required parameters
	-d database.tsv
	-e exons.tsv
	-t single/paired end fastq files.
	-i read_1.fastq.gz (or read_1.fastq). First read pair of the paired-end RNAseq fastq files. (Use with option -t paired)
	-j read_2.fastq.gz (or read_2.fastq). Second read pair of the paired-end RNAseq fastq files. (Use with option -t paired)
	-k read.fastq.gz (or read.fastq). Single-end RNAseq fastq file. (Use with option -t single)

	Optional parameters
	-c	Minimum number of split reads required (default: 1)
	-p	number of cores to use. (default: 1)
	-r	Approximate RAM memory (in GB) to be allocated per core. Value should be >=5. The more the RAM, the fast the software will run. (default: 5).
	-n	Prefix name of the output file. Extension .csv will be added to the output file. (default: result_$program_name.csv)
	-h	Prints help and provides instructions on running the software
EOF

use Getopt::Std;
getopts('d:e:c:r:t:i:j:k:n:p:h');

$pd=$opt_d; $pe=$opt_e; $pc=$opt_c; $pr=$opt_r; $pt=$opt_t; $pi=$opt_i; $pj=$opt_j; $pk=$opt_k; $pn=$opt_n; $ph=$opt_h; $pp=$opt_p;
$pa=$pd;

if($pp eq "" and $pa eq "" and $pe eq "" and $pc eq "" and $pr eq "" and $pt eq "" and $pi eq "" and $pj eq "" and $pk eq "" and $pn eq "") {  print "$usage"; exit; }
if($pp eq "" and $pa eq "" and $pe eq "" and $pc eq "" and $pr eq "" and $pt eq "" and $pi eq "" and $pj eq "" and $pk eq "" and $pn eq "" and $ph == 1) {  print "$usage"; exit; }
if($pa eq "-p" or $pa eq "-e" or $pa eq "-c" or $pa eq "-r" or $pa eq "-t" or $pa eq "-i" or $pa eq "-j" or $pa eq "-k" or $pa eq "-n") { print "$usage\n\n\tError! -d parameter is missing.\n"; exit; }
if($pe eq "-p" or $pe eq "-d" or $pe eq "-c" or $pe eq "-r" or $pe eq "-t" or $pe eq "-i" or $pe eq "-j" or $pe eq "-k" or $pe eq "-n") { print "$usage\n\n\tError! -e parameter is missing.\n"; exit; }
if($pi eq "-p" or $pi eq "-d" or $pi eq "-e" or $pi eq "-c" or $pi eq "-r" or $pi eq "-t" or $pi eq "-j" or $pi eq "-k" or $pi eq "-n") { print "$usage\n\n\tError! -i parameter is missing.\n"; exit; }
if($pj eq "-p" or $pj eq "-d" or $pj eq "-e" or $pj eq "-c" or $pj eq "-r" or $pj eq "-t" or $pj eq "-i" or $pj eq "-k" or $pj eq "-n") { print "$usage\n\n\tError! -j parameter is missing.\n"; exit; }
if($pk eq "-p" or $pk eq "-d" or $pk eq "-e" or $pk eq "-c" or $pk eq "-r" or $pk eq "-t" or $pk eq "-i" or $pk eq "-j" or $pk eq "-n") { print "$usage\n\n\tError! -k parameter is missing.\n"; exit; }
if($pn eq "-p" or $pn eq "-d" or $pn eq "-e" or $pn eq "-c" or $pn eq "-r" or $pn eq "-t" or $pn eq "-i" or $pn eq "-j" or $pn eq "-k") { print "$usage\n\n\tError! -n parameter is missing.\n"; exit; }
if($pr eq "-p" or $pr eq "-d" or $pr eq "-e" or $pr eq "-c" or $pr eq "-t" or $pr eq "-i" or $pr eq "-j" or $pr eq "-k" or $pr eq "-n") { print "$usage\n\n\tError! -r parameter is missing.\n"; exit; }
if($pp eq "-r" or $pp eq "-d" or $pp eq "-e" or $pp eq "-c" or $pp eq "-t" or $pp eq "-i" or $pp eq "-j" or $pp eq "-k" or $pp eq "-n") { print "$usage\n\n\tError! -p parameter is missing.\n"; exit; }
if($pn eq "") { $pn = "result_".$program_name; }
if($pp eq "") { $pp = 1; }
if($pp !~ /\d*/) { print "$usage\n\n\tError! -p $pp is invalid. It should be a number >= 1\n"; exit; }
if($pp < 1) { print "$usage\n\n\tError! -p $pp is invalid. It should be a number >= 1\n"; exit; }
if($pc eq "") { $pc = 1; }
if($pc !~ /\d*/) { print "$usage\n\n\tError! -c $pc is invalid. It should be a number >= 1\n"; exit; }
if($pc < 1) { print "$usage\n\n\tError! -c $pc is invalid. It should be a number >= 1\n"; exit; }
if($pa eq "") { print "$usage\n\n\tError! -d parameter is empty. Please provide the database.tsv file\n"; }
if($pe eq "") { print "$usage\n\n\tError! -e parameter is empty. Please provide the exons.tsv file\n"; }
if($pt eq "") { print "$usage\n\n\tError! -t parameter is missing.\n"; exit; }
if($pt ne "" and $pt ne "single" and $pt ne "paired") { print "$usage\n\n\tError! -t $pt is invalid. It should be either single or paired\n"; exit; }
if($pt eq "single" and $pk eq "") { print "$usage\n\n\tError! -k parameter is missing.\n"; exit; }
if($pt eq "single" and $pi ne "") { print "$usage\n\n\tError! -i parameter should not be used with -t single.\n"; exit; }
if($pt eq "single" and $pj ne "") { print "$usage\n\n\tError! -j parameter should not be used with -t single.\n"; exit; }
if($pt eq "paired" and $pi eq "") { print "$usage\n\n\tError! -i parameter is missing.\n"; exit; }
if($pt eq "paired" and $pj eq "") { print "$usage\n\n\tError! -j parameter is missing.\n"; exit; }
if($pt eq "paired" and $pk ne "") { print "$usage\n\n\tError! -k parameter should not be used with -t paired.\n"; exit; }
if($pr eq "") { $pr = 5; }
if($pr !~ /\d*/) { print "$usage\n\n\tError! -r $pr is invalid. It should be a number >= 5\n"; exit; }
if($pr < 5) { print "$usage\n\n\tError! -r $pr is invalid. It should be a number >= 5\n"; exit; }

# file checks
if(!(-e $pa)) { print "-d $pa file does not exist. Exiting...\n"; exit; }
if(!(-e $pe)) { print "-e $pe file does not exist. Exiting...\n"; exit; }
if($pt eq "paired") { 
if(!(-e $pi)) { print "-i $pi file does not exist. Exiting...\n"; exit; }
if(!(-e $pj)) { print "-j $pj file does not exist. Exiting...\n"; exit; }
}
if($pt eq "single") {
if(!(-e $pk)) { print "-k $pk file does not exist. Exiting...\n"; exit; }
}

$prpp = $pr*$pp;
$nsub = 180/$prpp;
if($nsub == int($nsub)) { $nsub_files = $nsub; }
else { $nsub_files = int($nsub)+1; }
#print "$pr/$pp/$nsub_files\n"; exit;
$ncores = $pp; $nc = $ncores;
$in = $pa;
$exon_file = $pe;
if($pt eq "single") { $read1 = $pk; $read2 = ""; }
if($pt eq "paired") { $read1 = $pi; $read2 = $pj; }
$cut_off = 1;
$cut_col = 14;
$ripgrep = "rg";
$prefix_output = "tmp.RTpred";
$pwd = `pwd`; chomp($pwd);
$jeach_side = 50;
$final_output_file = $pn.".csv";

##################################################
$compressed_option = "";
if($read1 =~ /\.gz/) { $compressed_option = "-z "; }
$cola1 = $cut_col + 2;
$cola2 = $cut_col + 3;
##################################################

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

print "Processing RTpred database...\n";
system "cut -f$cut_col $in |perl -p -e 's/#//g' > $tmpdir/$prefix_output.seqs";
system "cut -f$cola1 $in | grep -v '\\-NA\\-' | perl -p -e 's/#/\\n/g' >> $tmpdir/$prefix_output.seqs";
system "cut -f$cola2 $in | grep -v '\\-NA\\-' | perl -p -e 's/#/\\n/g' >> $tmpdir/$prefix_output.seqs";
system "sed -i -e '/^\$/d' $tmpdir/$prefix_output.seqs";
$get_each_side = `head -1 $tmpdir/$prefix_output.seqs`; chomp($get_each_side);
@ges = split(//, $get_each_side);
$bothside = @ges;
$eachside = int($bothside/2);
open(FF, "<$tmpdir/$prefix_output.seqs");
open(OO, ">$tmpdir/$prefix_output.seqs.revcomp");
while(<FF>)
{
	chomp($_);
	$_ =~ tr/ATGC/TACG/;
	$_ = reverse $_;
	print OO "$_\n";
}
close FF; close OO;
system "cat $tmpdir/$prefix_output.seqs $tmpdir/$prefix_output.seqs.revcomp | sort -u > $tmpdir/$prefix_output.fseqs";
system "rm $tmpdir/$prefix_output.seqs.revcomp $tmpdir/$prefix_output.seqs";

$number_of_files_after_splitting = $nc;
$inseq = "$tmpdir/$prefix_output.fseqs";
$wcl = `wc -l $inseq | perl -p -e 's/^\\s+//g' | cut -f1  -d ' '`; chomp($wcl);
$check = $wcl/$number_of_files_after_splitting;
$checkint = int($check);
if($check eq $checkint) { $lines_each_part = $check; }
else { $lines_each_part = int($check)+1; }
if($lines_each_part == 1) { $lines_each_part = 2; }
open(FF, "<$inseq");
@file = <FF>; close(FF);
$counter = 0;
for($i=0; $i<@file; $i++)
{
	$ii = $i+1;
	$line = $file[$i]; chomp($line);
	if($ii % $lines_each_part == 1)
	{
		$counter++;
		open(OO, ">$tmpdir/tmp.file.$counter");
		print OO "$line\n";
	}
	if($ii % $lines_each_part != 1 and $ii % $lines_each_part != 0)
	{
		print OO "$line\n";
	}
	if($ii % $lines_each_part != 1 and $ii % $lines_each_part == 0)
	{
		print OO "$line\n";
		close OO;
	}
}
system "rm $tmpdir/$prefix_output.fseqs";
$number_of_files_after_splitting = $counter;

print "Searching for potential read-throughs in input fastq file(s)...\n";
$etcat12 = ""; $etcattmp = ""; for($pqi=1; $pqi<=$number_of_files_after_splitting; $pqi++) { $etcat12 .= "$tmpdir/$prefix_output.12.$pqi "; $etcattmp .= "$tmpdir/tmp.file.$pqi "; }
use Parallel::ForkManager;
 
my @tempfiles=(1..$number_of_files_after_splitting);
my $pm = Parallel::ForkManager->new($number_of_files_after_splitting);
 
MULTI:
foreach my $etfile (@tempfiles)
{
	$pm->start and next MULTI;
	$cetfile = $etfile;
	$etfile = "tmp.file.$etfile";
	$ewcl = `wc -l $tmpdir/$etfile | perl -p -e 's/^\\s+//g' | cut -f1  -d ' '`; chomp($ewcl);
	$lewcl = $ewcl/$nsub_files; $ilewcl = int($lewcl);
	if($lewcl == $ilewcl) { $lines_all = $lewcl; }
	else { $lines_all = $ilewcl+1; }
	if($lines_all == 1) { $lines_all = 2; }

	open(EFF, "<$tmpdir/$etfile");
	@fetf = <EFF>; close(EFF);
	$etcounter = 0;
	for($eti=0; $eti<@fetf; $eti++)
	{
		$etii = $eti+1;
		$etline = $fetf[$eti]; chomp($etline);
		if($etii % $lines_all == 1)
		{
			$etcounter++;
			open(EOO, ">$tmpdir/$etfile.sub.$etcounter");
			print EOO "$etline\n";
		}
		if($etii % $lines_all != 1 and $etii % $lines_all != 0)
		{
			print EOO "$etline\n";
		}
		if($etii % $lines_all != 1 and $etii % $lines_all == 0)
		{
			print EOO "$etline\n";
			close EOO;
		}
	}
	$nsub_files = $etcounter;
	system "rm $tmpdir/$etfile";

	$etcatcom = "cat ";
	$etrmocom = "rm ";
	for($eti=1; $eti<=$nsub_files; $eti++)
	{
		$etinseq = "$tmpdir/$etfile.sub.$eti";
		system "$ripgrep -f $etinseq $compressed_option $read1 > $tmpdir/tmpo.$cetfile.1";
		if($read2 ne "") { system "$ripgrep -f $etinseq $compressed_option $read2 > $tmpdir/tmpo.$cetfile.2";
		system "cat $tmpdir/tmpo.$cetfile.1 $tmpdir/tmpo.$cetfile.2 > $tmpdir/o.$cetfile.$eti";
		}
		else {
			system "cat $tmpdir/tmpo.$cetfile.1 > $tmpdir/o.$cetfile.$eti";
		}
		open(RTMP, "<$tmpdir/o.$cetfile.$eti");
		open(WTMP, ">$tmpdir/wo.$cetfile.$eti");
		$flag_quality = 0; # 0 means sequence; 1 means rest
		while($rtmpline = <RTMP>)
		{
			chomp($rtmpline);
			@rtmpa = split(//, $rtmpline);
			for($rtmpi=0; $rtmpi<@rtmpa; $rtmpi++)
			{
				$rnt = uc($rtmpa[$rtmpi]);
				if($rnt ne "T" and $rnt ne "C" and $rnt ne "G" and $rnt ne "A" and $rnt ne "N") { $flag_quality = 1; last; }
			}
			if($flag_quality == 1) { next; }
			print WTMP "$rtmpline\n";
		}
		close RTMP;
		close WTMP;
		system "mv $tmpdir/wo.$cetfile.$eti $tmpdir/o.$cetfile.$eti";
		$etcatcom .= "$tmpdir/o.$cetfile.$eti ";
		system "rm $tmpdir/$etfile.sub.$eti $tmpdir/tmpo.$cetfile.1";
		if($read2 ne "") { system "rm $tmpdir/tmpo.$cetfile.2"; }
		$etrmocom .= "$tmpdir/o.$cetfile.$eti ";
	}
	$etcatcom .= "| sort -u > $tmpdir/$prefix_output.12.$cetfile";
	system "$etcatcom";
	system "$etrmocom";

	$pm->finish;
}
$pm->wait_all_children;

system "cat $etcat12 |sort -u > $tmpdir/$prefix_output.12";
system "rm $etcat12";
#-###########################################################
$inputforscript = "$tmpdir/$prefix_output.12";
$filetogrep = $inputforscript;
open(FF, "<$filetogrep");
@file = <FF>; close(FF);
%hash_patterns = ();
for($i=0; $i<@file; $i++)
{
	$ii = $i+1;
	$line = $file[$i]; chomp($line);
	@a = split(//, $line);
	for($j=0; $j<@a; $j++)
	{
		$jj = $j+$bothside-1;
		$pattern = join("", @a[$j..$jj]);
		$len_pattern = length($pattern);
		if($len_pattern != $bothside) { last; }
		$hash_patterns{$pattern}{$ii}++;
	}
}

open(FF, "<$in");
@annotation_file = <FF>; close(FF);
open(OOO, ">$tmpdir/$prefix_output.counts");
for($i=0; $i<@annotation_file; $i++)
{
	$line = $annotation_file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$mainid = $a[0];
	#@aa = split(/#/, $a[0]);
	#$g1 = $aa[0]; $g2 = $aa[2];
	$jseq = $a[$cut_col-1]; $jseq =~ s/#//g;
	#print "$a[$cut_col-1]\t$a[$cut_col+1]\t$a[$cut_col+2]\n";
	if($a[$cut_col+1] eq "-NA-" and $a[$cut_col+2] eq "-NA-") { $allseqs = $jseq; }
	if($a[$cut_col+1] eq "-NA-" and $a[$cut_col+2] ne "-NA-") { $allseqs = $jseq."#".$a[$cut_col+2]; }
	if($a[$cut_col+1] ne "-NA-" and $a[$cut_col+2] eq "-NA-") { $allseqs = $jseq."#".$a[$cut_col+1]; }
	if($a[$cut_col+1] ne "-NA-" and $a[$cut_col+2] ne "-NA-") { $allseqs = $jseq."#".$a[$cut_col+1]."#".$a[$cut_col+2]; }
	@ar_seqs = split(/#/, $allseqs);
	$count_seq = ""; $get_line_numbers1 = "";
	$count_seq_revcomp = ""; $get_line_numbers2 = "";
	for($se=0; $se<@ar_seqs; $se++)
	{
		$each_sequence = $ar_seqs[$se]; $reverse_each_sequence = $each_sequence;
		$reverse_each_sequence =~ tr/ATGC/TACG/; $reverse_each_sequence = reverse $reverse_each_sequence;
		$seq = $each_sequence;
		if(exists($hash_patterns{$seq}))
		{
			$its_counts = 0;
			foreach $each_line (sort {$a <=> $b} keys(%{$hash_patterns{$seq}}))
			{
				$get_line_numbers1 .= "$each_line;";
				$its_counts++;
			}
			$get_line_numbers1  =~  s/;$//g;
			$count_seq += $its_counts;
		}
	
		$seq = $reverse_each_sequence;
		if(exists($hash_patterns{$seq}))
		{
			$its_counts = 0;
			foreach $each_line (sort {$a <=> $b} keys(%{$hash_patterns{$seq}}))
			{
				$get_line_numbers2 .= "$each_line;";
				$its_counts++;
			}
			$get_line_numbers2  =~ s/;$//g;
			$count_seq_revcomp += $its_counts;
		}
	}
	
	if($count_seq eq "" and $count_seq_revcomp eq "") { next; }
	if($count_seq eq "") { $count_seq = 0; }
	if($count_seq_revcomp eq "") { $count_seq_revcomp = 0; }
	$count_total = $count_seq + $count_seq_revcomp;
	$get_line_numbers = $get_line_numbers1.";".$get_line_numbers2;
	$get_line_numbers  =~  s/;$//g;
	$get_line_numbers  =~  s/^;//g;
	print OOO "$mainid\t$count_total\t$get_line_numbers\n";
}
close OOO;
#-###########################################################
$input_file_to_modify = "$tmpdir/$prefix_output.12";
$input_file_with_line_numbers = "$tmpdir/$prefix_output.counts";
open(FF, "<$input_file_with_line_numbers");
@file = <FF>; close(FF);
%hash = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	($mainid, $counts, $line_numbers) = split(/\t/, $line);
	@a = split(/;/, $line_numbers);
	for($j=0; $j<@a; $j++)
	{
		push(@{$hash{$a[$j]}}, $mainid);
	}
}

open(FF, "<$input_file_to_modify");
@file = <FF>; close(FF);
open(OO, ">$input_file_to_modify.revised");
$sequence_number = 0;
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	$line_number = $i+1;
	$isoforms_it_supports = @{$hash{$line_number}};
	if($isoforms_it_supports > $cut_off)
	{
		next;
	}
	$sequence_number++;
	$header = ">S$sequence_number";
	print OO "$line\n";
}
close OO;

system "mv $input_file_to_modify $input_file_to_modify.prev; cp $input_file_to_modify.revised $input_file_to_modify";
system "mv $input_file_with_line_numbers $input_file_with_line_numbers.prev";
#system "rm $input_file_to_modify.revised $input_file_with_line_numbers.prev";
#-###########################################################
$inputforscript = "$tmpdir/$prefix_output.12";
$filetogrep = $inputforscript;
open(FF, "<$filetogrep");
@file = <FF>; close(FF);
%hash_patterns = ();
for($i=0; $i<@file; $i++)
{
	$ii = $i+1;
	$line = $file[$i]; chomp($line);
	@a = split(//, $line);
	for($j=0; $j<@a; $j++)
	{
		$jj = $j+$bothside-1;
		$pattern = join("", @a[$j..$jj]);
		$len_pattern = length($pattern);
		if($len_pattern != $bothside) { last; }
		$hash_patterns{$pattern}{$ii}++;
	}
}

open(OOO, ">$tmpdir/$prefix_output.counts");
for($i=0; $i<@annotation_file; $i++)
{
	$line = $annotation_file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$mainid = $a[0];
	#@aa = split(/#/, $a[0]);
	#$g1 = $aa[0]; $g2 = $aa[2];
	$jseq = $a[$cut_col-1]; $jseq =~ s/#//g;
	#print "$a[$cut_col-1]\t$a[$cut_col+1]\t$a[$cut_col+2]\n";
	if($a[$cut_col+1] eq "-NA-" and $a[$cut_col+2] eq "-NA-") { $allseqs = $jseq; }
	if($a[$cut_col+1] eq "-NA-" and $a[$cut_col+2] ne "-NA-") { $allseqs = $jseq."#".$a[$cut_col+2]; }
	if($a[$cut_col+1] ne "-NA-" and $a[$cut_col+2] eq "-NA-") { $allseqs = $jseq."#".$a[$cut_col+1]; }
	if($a[$cut_col+1] ne "-NA-" and $a[$cut_col+2] ne "-NA-") { $allseqs = $jseq."#".$a[$cut_col+1]."#".$a[$cut_col+2]; }
	@ar_seqs = split(/#/, $allseqs);
	$count_seq = ""; $get_line_numbers1 = "";
	$count_seq_revcomp = ""; $get_line_numbers2 = "";
	for($se=0; $se<@ar_seqs; $se++)
	{
		$each_sequence = $ar_seqs[$se]; $reverse_each_sequence = $each_sequence;
		$reverse_each_sequence =~ tr/ATGC/TACG/; $reverse_each_sequence = reverse $reverse_each_sequence;
		$seq = $each_sequence;
		if(exists($hash_patterns{$seq}))
		{
			$its_counts = 0;
			foreach $each_line (sort {$a <=> $b} keys(%{$hash_patterns{$seq}}))
			{
				$get_line_numbers1 .= "$each_line;";
				$its_counts++;
			}
			$get_line_numbers1  =~  s/;$//g;
			$count_seq += $its_counts;
		}
	
		$seq = $reverse_each_sequence;
		if(exists($hash_patterns{$seq}))
		{
			$its_counts = 0;
			foreach $each_line (sort {$a <=> $b} keys(%{$hash_patterns{$seq}}))
			{
				$get_line_numbers2 .= "$each_line;";
				$its_counts++;
			}
			$get_line_numbers2  =~ s/;$//g;
			$count_seq_revcomp += $its_counts;
		}
	}
	
	if($count_seq eq "" and $count_seq_revcomp eq "") { next; }
	if($count_seq eq "") { $count_seq = 0; }
	if($count_seq_revcomp eq "") { $count_seq_revcomp = 0; }
	$count_total = $count_seq + $count_seq_revcomp;
	$get_line_numbers = $get_line_numbers1.";".$get_line_numbers2;
	$get_line_numbers  =~  s/;$//g;
	$get_line_numbers  =~  s/^;//g;
	print OOO "$mainid\t$count_total\t$get_line_numbers\n";
}
close OOO;
#-###########################################################
print "Generating the result file...\n";
$inc = "$tmpdir/$prefix_output.counts";
%hash_final = ();
for($i=0; $i<@annotation_file; $i++)
{
	$line = $annotation_file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$mid = $a[0];
	$hash_final{$mid} = $line;
}

open(FF, "<$exon_file"); @file = <FF>; close(FF);
%hash_eseq = ();
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$eid = $a[0]; $eseq = $a[1];
	$hash_eseq{$eid} = $eseq;
}


open(FF, "<$inc");
@file = <FF>; close(FF);
open(OOFI, ">$final_output_file");
$col_split_reads = 11;
print OOFI "GeneName1\tGeneName2\tChr1\tBreakpoint1\tStrand1\tChr2\tBreakpoint2\tStrand2\tEnsemblGeneID1\tEnsemblGeneID2\tsplit_reads\tJunctionSequence\tfull_exon_seqs\n";
for($i=0; $i<@file; $i++)
{
	$line = $file[$i]; chomp($line);
	@a = split(/\t/, $line);
	$mainid = $a[0]; $count = $a[1];
	($gid1, $e1, $gid2,$e2) = split(/#/, $mainid);
	$get_final_line = $hash_final{$mainid};
	@ga = split(/\t/, $get_final_line);
	$name1 = $ga[4]; $name2 = $ga[10];
	($exid1, $len1, $num1) = split(/#/, $name1);
	($exid2, $len2, $num2) = split(/#/, $name2);
	$get_eseq1 = $hash_eseq{$exid1}; $get_eseq1 = lc($get_eseq1); $fullj = $get_eseq1; $get_eseq1 =~ s/^.*(.{$jeach_side})$/$1/g;
	$get_eseq2 = $hash_eseq{$exid2}; $fullj .= $get_eseq2; $get_eseq2 =~ s/(^.{$jeach_side}).*$/$1/g;
	$jseq = $get_eseq1.$get_eseq2;
 	$gnames = $ga[14];
	($gn1, $gn2)  = split(/#/, $gnames);
	$chr1 = $ga[1]; $chr2 = $ga[7];
	$strand1 = $ga[6]; $strand2 = $ga[12];
	$start1 = $ga[2]; $end1 = $ga[3];
	$start2 = $ga[8]; $end2 = $ga[9];
	if($strand1 eq "+") { $b1 = $end1; $b2 = $start2 + 1; }
	if($strand1 eq "-") { $b1 = $start1 + 1; $b2 = $end2; }
	if($count < $pc) { next; }
	if($gn1 eq "") { $gn1 = "N/A"; }
	if($gn2 eq "") { $gn2 = "N/A"; }
	print OOFI "$gn1\t$gn2\t$chr1\t$b1\t$strand1\t$chr2\t$b2\t$strand2\t$gid1\t$gid2\t$count\t$jseq\t$fullj\n";
}
close OOFI;
$lfof = `wc -l $final_output_file | perl -p -e 's/^\\s+//g' | cut -f1  -d ' '`; chomp($lfof); if($lfof<2) { `echo 'Chimeric RNA not found' >> $final_output_file`; }
#-###########################################################
#system "sort -nrk$col_split_reads $final_output_file > $final_output_file.s; mv $final_output_file.s $final_output_file";
system "perl -pi -e 's/\t/,/g' $final_output_file";
print "Deleting temporary directory ./$tmpdir\n";
use File::Path; rmtree $tmpdir;
print "\nResult file is $final_output_file\n";
