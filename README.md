# pal4nal.pl
similar as what pal2nal.pl do, i.e. change protein alignment into CDS alignment, except some difference.

Story start like this: 
  CDS alignment for gene is very important for my research, you'll need that for dN/dS calculation and so on.
  Currently there are two choices to fix this: 
  1)pal2nal.pl (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1538804/)
  2)MACSE (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022594)

The pain goes like this: 
  MACSE take unbearable time when dealing with large alignments.
  pal2nal.pl is fast, however, it have zero tolerance to cds-pep sequence inconsistency. What make it worse? this inconsistency exists even with the data downloaded from Ensembl/NCBI(genome).

The last straw crashed me: 
  Recently I found that my CDS alignments, gap-removed using Gblocks -t=c, contain lots of premature stop codons when being    translated into pretein (maybe I use Gblocks wrong).

How does pal4nal.pl do it?
  1)readin protein alignment and corresponding CDS sequence
  2)align each protein alignment with its corresponding CDS sequence using genewise
  3)return the CDS sequence according to the genewise alignemnt

At last:
  Bugs may exist, I would be very appreciated if you could let me know.

Many thanks and best wishes.

                                                     Kang @ Würzburg
                                                     2019-3-11
