#!/usr/bin/perl -w

use threads; 

my $usage=<<EOF;
--------------------------------------
Self-introduction:
	I do what pal2nal.pl do, but in a different way. pal2nal.pl fails and stops when frameshit exist (ERROR:inconsistency), while I'll inform you that, but do not fail out.
	I use genewise to retrieve cds for each protein, this means, you can even feed me with intron-containning DNA sequences.

Request:
	Pleaes let me have genewise program in the PATH
	Please let the protein and its corresponding cds in the same name
	Please let the alignment and sequence file in fasta format

Usage: perl $0 pep.align cds.fa (-n num)
	-n set the threads number, default 15
	NO REDIRECTION ">" is needed, the result file would be *.pal4nal
	while warning message would be revealed on the screen, or you can re-direct it

												Du Kang 2019-3-10

v2.
	Genewise is much more time-consuming than pal2nal.pl, so here I took adventage of pal2nal.pl.
	For those with aa-cds consistency, they will be well looked after by the code from pal2nal.pl.
	For those with aa-cds inconsistency, let's feed them to Genewise

													2019-3-19
--------------------------------------
EOF

$ARGV[0] or die $usage;
$ARGV[1] or die $usage;

$basename= $ARGV[0]=~/.*\/(.*)/? $1 : $ARGV[0];
$outfile="$basename.pal4nal";
unlink $outfile and print "Prompt: file $outfile exist, and is removed\n" if -e $outfile;

$max_process = 15;
foreach $i (0..@ARGV-1){
        $max_process= $ARGV[$i+1] if $ARGV[$i] eq "-n";
}

print ">>>>>>>>>> start for $ARGV[0]\n";

##################### readin sequences

open PEP, $ARGV[0] or die $!;
while (<PEP>) {
	chomp;
	if (/>/) {
		s/>//;
		$name=$_;
	}else{
		$pseq{$name} .= uc $_;
		$pseq{$name} =~ s/-//g;
		$pseqal{$name} .= uc $_;

	}
}
close PEP;

open DNA, $ARGV[1] or die $!;
while (<DNA>) {
	chomp;
	if (/>/) {
		s/>//;
		$name=$_;
	}else{
		$dseq{$name} .= $_;
	}
}
close DNA;

foreach $key (keys %pseq){
        $pseq{$key} or die "ERROR: no protein sequence for $key!!!\n";
        $dseq{$key} or die "ERROR: no cds sequence for $key!!!\n";
}

foreach $key (keys %pseq){	# for each protein sequence in this alignment file

	############################## to see if the aa-cds consistent or not
	%codonout = &pn2codon($pseq{$key}, $dseq{$key}, 1);

		#$$$$$$$$$$$$$$$$$$$ if consistent, pal2nal.pl
	if ($codonout{'result'} == 1 || $codonout{'result'} == 2) {
		$codonseq = lc $codonout{'codonseq'};

		$cseqal = $pseqal{$key};
		$cseqal =~ s/-/---/g;
		$l=length $pseq{$key};
		for $i (0..$l-1){
			$site=substr ($pseq{$key},$i,1);
			if ($site =~ /[A-Z]/ || $site eq '*') {
				$csite = substr($codonseq,0,3);
				substr($codonseq,0,3) = "";
				$cseqal =~ s/$site/$csite/;
			}
		}

		open OUT, ">>$outfile" or die $!;
		print OUT ">$key\n$cseqal\n";
		close OUT;

		print "pal2nal.pl says:\n\t\t@{$codonout{'message'}}\n" if $codonout{'message'};

	}else{
		#$$$$$$$$$$$$$$$$$$$$ if inconsistent, launch Genewise in parallel
		$wisepep="$key.190319.pep";
		$wisecds="$key.190319.cds";
		`echo ">&P_$key\n$pseq{$key}" >$wisepep`;
		`echo ">&C_$key\n$dseq{$key}" >$wisecds`;
		sleep rand 2 while threads->list(threads::running) >= $max_process;
		$thread{$key} = threads->new(sub{`genewise $wisepep $wisecds -quiet`});
		push @wise, $key;
	}
}

#################### recycle Genewise if launched, and parse the result
if (@wise) {
	foreach $key (@wise){
		$wise{$key} = $thread{$key}->join();
	}

	foreach $key (keys %wise){
		$inconsistency=0;
		$wisepep="$key.190319.pep";
                $wisecds="$key.190319.cds";	
		unlink ($wisepep, $wisecds);

		#$$$$$$$$$$$$$ parse genewise outcomes
		$wise{$key} =~ s/(.*?)\/\/.*/$1/s;
		@wise=split /\n/, $wise{$key};
		$queseq="";
		$tarpep="";
		$tarcds1="";
		$tarcds2="";
		$tarcds3="";
		foreach $i (0..@wise-1){
			$_=$wise[$i];
			@_=split;
			if (/^&P_\S*?\s+(\d+?)\s.*/){
				$start=$1 if ! defined $start;
				$queseq .= substr ($_, 21, 49);
        		}elsif (/^&C_/){
				$tarcds1 .= substr ($_, 21, 49);
                		$tarpep .= substr ($wise[$i-1], 21, 49);
                		$tarcds2 .= substr ($wise[$i+1], 21, 49);
                		$tarcds3 .= substr ($wise[$i+2], 21, 49);
			}
		}
		#	print "$queseq\n$tarpep\n$tarcds1\n$tarcds2\n$tarcds3\n";

		#$$$$$$$$$$$$$$$ output cds alignment

		$cseqal=$pseqal{$key};	# now we start to transform $cseqal

		$count=0;
		for $site ($cseqal=~/(.)/g){
			if ($site=~/[A-Z]/){
				$count++;
				$cseqal =~ s/$site/-/ if $count < $start;	# no cds corresponds to the head part, change them to "-"s
			}
		}

		$cseqal =~ s/-/---/g;

		$l=length $queseq;

		for $i (0..$l-1){
			$site=substr ($queseq,$i,1);
			$csite=substr($tarcds1,$i,1).substr($tarcds2,$i,1).substr($tarcds3,$i,1);
			$csitep=substr ($tarpep,$i,1);
			if ($site=~/[A-Z]/) {
				if ($csitep =~ /X/i){
					$inconsistency=1;
					$item= $csitep eq "X" ? "a premature stop condon" : "an unknown aa";
					print "Genewise Warning: $key in $ARGV[0] contains $item at site $i\n";
					$cseqal =~ s/$site/---/;
				} elsif ($csitep eq "!") {
					$inconsistency=1;
					print "Genewise Warning: $key in $ARGV[0] presend a frameshift with its cds at site $i\n";
					$cseqal =~ s/$site/---/;
				} elsif ($csitep eq "-") {
					$inconsistency=1;
					print "Genewise Warning: $key in $ARGV[0] have no cds corresponding at site $i\n";
					$cseqal =~ s/$site/---/;
				} elsif ($csitep ne $site) {
					$inconsistency=1;
					print "Genewisewise Warning: $key in $ARGV[0] is inconsistent with its cds sequence at site $i: $site <=> $csite\n";
					$cseqal =~ s/$site/$csite/;
				} else {
					$cseqal =~ s/$site/$csite/;
				}
			}
		}
		print "\n".$wise{$key} if $inconsistency==1;
		$cseqal =~ s/[A-Z]/---/g;	# no cds corresponds to the tail part, change them to "---"s

		open OUT, ">>$outfile" or die $!;
		print OUT ">$key\n$cseqal\n";
		close OUT;
	}
}


print "<<<<<<<<<< end for $ARGV[0] \n\n";


############################ sub ##############################
sub pn2codon {
    #    pn2codon v6
    #
    #    input:   $pep    protein sequence
    #                         termination -> "_" or "*";
    #                         frameshift  -> digit
    #                         "-" or "."  -> gap
    #             $nuc    DNA or RNA sequence (lower/upper case letters)
    #
    #             $codontable  (corresponds to codon tables used in GenBank
    #                1  Universal code
    #                2  Vertebrate mitochondrial code
    #                3  Yeast mitochondrial code
    #                4  Mold, Protozoan, and Coelenterate Mitochondrial code
    #                   and Mycoplasma/Spiroplasma code
    #                5  Invertebrate mitochondrial
    #                6  Ciliate, Dasycladacean and Hexamita nuclear code
    #                9  Echinoderm and Flatworm mitochondrial code
    #               10  Euplotid nuclear code
    #               11  Bacterial, archaeal and plant plastid code
    #               12  Alternative yeast nuclear code
    #               13  Ascidian mitochondrial code
    #               14  Alternative flatworm mitochondrial code
    #               15  Blepharisma nuclear code
    #               16  Chlorophycean mitochondrial code
    #               21  Trematode mitochondrial code
    #               22  Scenedesmus obliquus mitochondrial code
    #               23  Thraustochytrium mitochondrial code
    #
    #    return:  hash
    #                $retval{'codonseq'}: codon seq (w/o gap)
    #                $retval{'message'}:  error/warning messame (array)
    #                $retval{'result'}:   1: OK, 2: mismatch, -1: no match found
    #
    #
    #                                      05/05/2002    Mikita Suyama
    #                                      12/07/2004
    #
    #    v2                                22/04/2005
    #    - reverse translation table (%p2c) was replaced
    #    - if there is no exact match, try fragment anchors
    #    - return value is changed to hash
    #
    #    v3                                26/03/2006
    #    - codon table for vertebrate mitochondria (vmitochondria) is added
    #
    #    v4
    #    - from: 10 residue window (incl '-')
    #      to:   10 residue window (excl '-')
    #                                      19/06/2006
    #    v5
    #    - bug fix:
    #        just before 'does not correspond'
    #        /$p2c{$tmpaa}/  ->   /$p2c{$tmpaa}/i
    #                                           ^
    #    - the first Met can be (A|G|C|R)TG
    #                                      01/08/2006
    #    v5.1
    #    - but fix:
    #        if ($tmpcodon !~ /((A|C|G|R)TG)/i) {
    #                                        ^
    #                                      2009/06/22
    #    v6
    #    - start Met -> aa-code "B" in %p2c
    #    - NCBI codon tables
    #                                      2011/11/13
    #      


    local($pep, $nuc, $codontable) = @_;


    local(%p2c);
    local($peplen, $qcodon, $codon);
    local($tmpaa, $message, @qcodon, @fncodon, $wholecodon);
    local($i, $j, $anclen, @anchor, $peppos, $tmpcodon, $codonpos);

    local(%retval);


    if ($codontable == 1) {
        #-----------#
        # Universal Code (transl_table=1)
        #-----------#
        %p2c = (
            "B" => "((U|T|C|Y|A)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 2) {
        #--------------------------#
        # Vertebrate Mitochondrial Code (transl_table=2)
        #--------------------------#
        %p2c = (
            "B" => "((A(U|T).)|G(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "(((U|T)A(A|G|R))|(AG(A|G|R)))",
            "*" => "(((U|T)A(A|G|R))|(AG(A|G|R)))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 3) {
        #--------------------------#
        # Yeast Mitochondrial Code (transl_table=3)
        #--------------------------#
        %p2c = (
            "B" => "(A(U|T)(A|G|R))",
            "L" => "((U|T)(U|T)(A|G|R))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "((AC.)|(C(U|T).))",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 4) {
        #-----------#
        # Mold, Protozoan, and Coelenterate Mitochondrial Code
        # and Mycoplasma/Spiroplasma Code (transl_table=4)
        #-----------#
        %p2c = (
            "B" => "((A(U|T).)|((U|T)(U|T)(A|G|R))|(C(U|T)G)|(G(U|T)G))",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 5) {
        #-----------#
        # Invertebrate Mitochondrial Code (transl_table=5)
        #-----------#
        %p2c = (
            "B" => "((A(U|T).)|((U|T|A|G|R)(U|T)G))",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 6) {
        #-----------#
        # Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((T|U)GA)",
            "*" => "((T|U)GA)",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "((CA(A|G|R))|((U|T)A(A|G|R)))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 9) {
        #-----------#
        # Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
        #-----------#
        %p2c = (
            "B" => "((A|G|R)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AAG)",
            "N" => "(AA(U|T|C|Y|A))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 10) {
        #-----------#
        # Euplotid Nuclear Code (transl_table=10)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 11) {
        #-----------#
        # Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
        #-----------#
        %p2c = (
            "B" => "((A(U|T)G)|(.(U|T)G))",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 12) {
        #-----------#
        # Alternative Yeast Nuclear Code (transl_table=12)
        #-----------#
        %p2c = (
            "B" => "((A|C)(U|T)G)",
            "L" => "((C(U|T)(U|T|C|Y|A))|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y))|(C(U|T)G))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 13) {
        #-----------#
        # Ascidian Mitochondrial Code (transl_table=13)
        #-----------#
        %p2c = (
            "B" => "((T|U|A|G|R)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "((GG.)|(AG(A|G|R)))",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 14) {
        #-----------#
        # Alternative Flatworm Mitochondrial Code (transl_table=14)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)AG)",
            "*" => "((U|T)AG)",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AAG)",
            "N" => "(AA(U|T|C|Y|A))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y|A))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)G(G|A|R))",
            "X" => "...",
        );
    } elsif ($codontable == 15) {
        #-----------#
        # Blepharisma Nuclear Code (transl_table=15)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)AA)|((T|U)GA))",
            "*" => "(((U|T)AA)|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "((CA(A|G|R))|((U|T)AG))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 16) {
        #-----------#
        # Chlorophycean MitochondrialCode (transl_table=16)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R))|((U|T)AG))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)AA)|((T|U)GA))",
            "*" => "(((U|T)AA)|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 21) {
        #-----------#
        # Trematode Mitochondrial Code (transl_table=21)
        #-----------#
        %p2c = (
            "B" => "((A|G|R)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AAG)",
            "N" => "(AA(U|T|C|Y|A))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 22) {
        #-----------#
        # Scenedesmus obliquus mitochondrial Code (transl_table=22)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R))|((T|U)AG))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C(U|T|C|Y|G))|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)(C|A|G|R)A)",
            "*" => "((U|T)(C|A|G|R)A)",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 23) {
        #-----------#
        # Thraustochytrium Mitochondrial Code (transl_table=23)
        #-----------#
        %p2c = (
            "B" => "(((A|G|R)(U|T)G)|(A(U|T)(U|T)))",
            "L" => "((C(U|T).)|((U|T)(U|T)G))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA)|((T|U)(T|U)A))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA)|((T|U)(T|U)A))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    }


    #---------------------------------------------------------------#
    # make codon sequence, $qcodon,  with all possible combinations
    #---------------------------------------------------------------#

    $peplen = length($pep);
    $qcodon="";
    foreach $i (0..$peplen - 1) {
        $peppos = $i + 1;
        $tmpaa = substr($pep, $i, 1);
        if ($tmpaa =~ /[ACDEFGHIKLMNPQRSTVWY_\*XU]/) {
            if ($qcodon !~ /\w/ && substr($pep, $i, 1) eq "M") {

                #---------------#
                # the first Met
                #---------------#

                $qcodon .= $p2c{"B"};
            } else {
                $qcodon .= $p2c{substr($pep, $i, 1)};
            }
        } elsif ($tmpaa =~ /\d/) {
            $qcodon .= "." x $tmpaa;
        } elsif ($tmpaa =~ /[-\.]/) {
            # nothing to do
        } else {
            $message = "pepAlnPos $peppos: $tmpaa unknown AA type. Taken as 'X'";
            push(@{$retval{'message'}}, $message);
            $qcodon .= $p2c{'X'};
        }
    }
    # print "$qcodon\n";


    #-----------------------------#
    # does $nuc contain $qcodon ?
    #-----------------------------#

    if ($nuc =~ /$qcodon/i) {
        $codon = $&;

        $retval{'codonseq'} = $codon;
        $retval{'result'} = 1;

    } else {
        #-------------------#
        # make 10 aa anchor
        #-------------------#

#        undef(@{$retval{'message'}});

#        $modpep = $pep;
#        1 while $modpep =~ s/(.{10})(.{10,})/$1\n$2/;
#        @anchor = split(/\n/, $modpep);

        undef(@preanchor);
        undef($tmpanc);
        $nanc = 0;
        foreach $i (0..$peplen - 1) {
            $tmpaa = substr($pep, $i, 1);
            $tmpanc .= $tmpaa;
            ++$nanc if ($tmpaa !~ /-/);
            if ($nanc == 10 || $i == $peplen - 1) {
                push(@preanchor, $tmpanc);
                undef($tmpanc);
                $nanc = 0;
            }
        }
        undef(@anchor);
        $lastanchorlen = length($preanchor[-1]);
        if ($lastanchorlen < 10) {
            foreach $i (0..$#preanchor - 1) {
                if ($i < $#preanchor - 1) {
                    push(@anchor, $preanchor[$i]);
                } else {
                    push(@anchor, "$preanchor[$i]$preanchor[$i + 1]");
                }
            }
        } else {
            @anchor = @preanchor;
        }

        undef($wholecodon);
        foreach $i (0..$#anchor) {
            # print "    $anchor[$i]\n";
            $anclen = length($anchor[$i]);
            $qcodon[$i]="";
            undef(@fncodon);
            foreach $j (0..$anclen - 1) {
                $peppos = $i * 10 + $j + 1;
                $tmpaa = substr($anchor[$i], $j, 1);
                if ($tmpaa =~ /[ACDEFGHIKLMNPQRSTVWY_\*XU]/) {
                    if ($i == 0 && $qcodon[$i] !~ /\w/ && $tmpaa eq "M") {

                        #----------------------------------#
                        # the first Met can be AGT|GTG|CTG
                        #----------------------------------#

                        $qcodon[$i] .= "((A|C|G|R)TG)";
                    } else {
                        $qcodon[$i] .= $p2c{$tmpaa};
                    }
                    $fncodon[$i] .= $p2c{'X'};
                } elsif ($tmpaa =~ /\d/) {
                    $qcodon[$i] .= "." x $tmpaa;
                    $fncodon[$i] .= "." x $tmpaa;
                } elsif ($tmpaa =~ /[-\.]/) {
                    # nothing to do
                } else {
                    #del $message = "pepAlnPos $peppos: $tmpaa unknown AA type. Replaced by 'X'";
                    #del push(@{$retval{'message'}}, $message);
                    $qcodon[$i] .= $p2c{'X'};
                    $fncodon[$i] .= $p2c{'X'};
                }
            }
            if ($nuc =~ /$qcodon[$i]/i) {
                $wholecodon .= $qcodon[$i];
            } else {
                $wholecodon .= $fncodon[$i];
            }
        }

        if ($nuc =~ /$wholecodon/i) {
            $codon = $&;
            $codonpos = 0;
            $tmpnaa = 0;
            foreach $i (0..$peplen - 1) {
                $peppos = $i + 1;
                $tmpaa = substr($pep, $i, 1);
                undef($tmpcodon);
                if ($tmpaa !~ /\d/ && $tmpaa !~ /-/) {
                    ++$tmpnaa;
                    $tmpcodon = substr($codon, $codonpos, 3);
                    $codonpos += 3;
                    if ($tmpnaa == 1 && $tmpaa eq "M") {
                        if ($tmpcodon !~ /((A|C|G|R)TG)/i) {
                            $message = "pepAlnPos $peppos: $tmpaa does not correspond to $tmpcodon";
                            push(@{$retval{'message'}}, $message);
                        }
                    } elsif ($tmpcodon !~ /$p2c{$tmpaa}/i) {
                        $message = "pepAlnPos $peppos: $tmpaa does not correspond to $tmpcodon";
                        push(@{$retval{'message'}}, $message);
                    }
                } elsif ($tmpaa =~ /\d/i) {
                    $tmpcodon = substr($codon, $codonpos, $tmpaa);
                    $codonpos += $tmpaa;
                }
                # print "$tmpaa    $tmpcodon\n";
            }
#        print    "$codon\n";

            $retval{'codonseq'} = $codon;
            $retval{'result'} = 2;

        } else {

            $retval{'result'} = -1;

        }

    }

    return(%retval);
}
