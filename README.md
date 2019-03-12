# pal4nal.pl
similar as what pal2nal.pl do, i.e. change protein alignment into CDS alignment, except some difference.

Story start like this:    
  CDS alignment for gene is very important for my research, you'll need that for dN/dS calculation and so on.   
  Currently there are two choices to do this:    
  1)pal2nal.pl (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1538804/)      
  2)MACSE (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022594)   

The pain goes like this:    
  MACSE take unbearable time when dealing with large alignments.   
  pal2nal.pl is fast, however, it have zero tolerance to cds-pep sequence inconsistency, thus on this situation give an empty output. What make it worse? cds-pep inconsistency exists even with the data downloaded from Ensembl/NCBI(genome).   

The last straw crashed me:    
  Recently I found that my CDS alignments, gap-removed using Gblocks -t=c, contain lots of premature stop codons when being    translated into proteins (maybe I use Gblocks wrong). Can I removing gaps from protein alignment first then translate them to CDS alignment? (I have tested it afterward, better not do this!!!! genewise failed to retrieve the correct condon, however, the good news is, when using pal4nal.pl with "Gblock -t=c" afterwards, premature stop codons never pop out again #@2019-3-12)

How does pal4nal.pl do it?   
  1)readin protein alignment and corresponding CDS sequence   
  2)align each protein with its corresponding CDS sequence using genewise    
  3)retrieve the CDS sequence according to the genewise alignemnt for each protein in the alignment   

At last:    
  Bugs may exist, I would be very appreciated if you could let me know.   

Many thanks and best wishes.   

                                                     Kang (dukang1117@outlook.com)
                                                     2019-3-11
