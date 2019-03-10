#!/usr/bin/perl -w

use threads; 

my $usage=<<EOF;
--------------------------------------
Self-introduction:
I do what pal2nal.pl do, but in a different way. pal2nal.pl fails and stops when frameshit exist (ERROR:inconsistency), while I'll inform you that, but do not fail out.
I use genewise to retrieve cds for each protein, this means, you can even feed me with intron-containning DNA sequences.

Request:
Pleaes let me have genewise program in the PATH;
Please let the protein and its corresponding cds in the same name.

Usage: perl $0 pep.align cds.fa (-n num)
-n set the threads number

the result file would be pep.align.pal4nal
while warning message would be revealed on the screen, or you can re-direct it

						           Du Kang 2019-3-10
--------------------------------------
EOF

$ARGV[0] or die $usage;
$ARGV[1] or die $usage;

$basename= $ARGV[0]=~/.*\/(.*)/? $1 : $ARGV[0];
$outfile="$basename.pal4nal";
`rm $outfile` and print "Prompt: exist $outfile removed\n" if -e $outfile;

$max_process = 15;
foreach $i (0..@ARGV-1){
        $max_process= $ARGV[$i+1] if $ARGV[$i] eq "-n";
}

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

############################## genewise/estwise in parallel

foreach $key (keys %pseq){
	`echo ">$key\n$pseq{$key}" >$key.pep`;
	`echo ">$key\n$dseq{$key}" >$key.cds`;
	$thread{$key} = threads->new(sub{`genewise $key.pep $key.cds -quiet`});
}

foreach $key (keys %pseq){
	$wise{$key} = $thread{$key}->join();
}

foreach $key (keys %pseq){
	$inconsistency=0;	
	`rm $key.pep $key.cds`;

	############################# parse genewise outcomes
	$wise{$key} =~ s/(.*?)\/\//$1/;
	@wise=split /\n/, $wise{$key};
	$queseq="";
	$tarpep="";
	$tarcds1="";
	$tarcds2="";
	$tarcds3="";
	foreach $i (0..@wise-1){
		$_=$wise[$i];
		@_=split;
		if (/^Query protein/){
			$queflag= length($_[-1])>10? substr($_[-1],0,10) : $_[-1];
		}
		if (/^Target Sequence/){
			$tarflag= length($_[-1])>10? substr($_[-1],0,10) : $_[-1];
		}
		if (defined $queflag and defined $tarflag){
			if (/^$queflag.*\s(\d+)\s.*[A-Z]+/ and !/^$queflag.*\s\d+\s.*[a-z]+/){
				$start=$1 if ! defined $start;
                		$queseq .= substr ($_, 21, 49);
        		}elsif (/^$tarflag.*\d\s.*[a-z]+/){
				$tarcds1 .= substr ($_, 21, 49);
                        	$tarpep .= substr ($wise[$i-1], 21, 49);
                        	$tarcds2 .= substr ($wise[$i+1], 21, 49);
                        	$tarcds3 .= substr ($wise[$i+2], 21, 49);
			}
		}
	}
#	print "$queseq\n$tarpep\n$tarcds1\n$tarcds2\n$tarcds3\n";

	########################### output cds alignment

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
		if ($site=~/[A-Z]/) {
			if (substr ($tarpep,$i,1) =~ /X/i){
				$inconsistency=1;
				print ">>>Warning: $key in $ARGV[0] contains a premature stop condon at site $i\n";
				$cseqal =~ s/$site/---/;
			} elsif (substr ($tarpep,$i,1) eq "!") {
				$inconsistency=1;
				print ">>>Warning: $key in $ARGV[0] presend a frameshift with its cds at site $i\n";
				$cseqal =~ s/$site/---/;
			} elsif (substr ($tarpep,$i,1) eq "-") {
				$inconsistency=1;
				print ">>>Warning: $key in $ARGV[0] have no cds corresponding at site $i\n";
				$cseqal =~ s/$site/---/;
			} elsif (substr ($tarpep,$i,1) ne $site) {
				$inconsistency=1;
				print ">>>Warning: $key in $ARGV[0] is inconsistent with its cds sequence at site $i: $site <=> $csite\n";
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
}

print "^^^^^^^^^^^^^^^^^^^^^ end ^^^^^^^^^^^^^^^^^^^^^^\n";
