# Barkemeyer_Python_Portfolio_2
This is the portfoltio of python code that I learned in my Advanced Python course.
```python
# Sequence Objects
# In terminal: pip install biopython

from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" %(index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence:
print(len(my_seq))
```

    5



```python
# we can access the different positions in seq objects:
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
# We can do a .count (does not count overlapping)

Seq("AAAA").count("AA")
```




    2




```python
# another example
my_seq = "GATACGATTGCATGCAGCATAAACGTA"
```


```python
len(my_seq)
```




    27




```python
# To find the count of Gs in the sequence
my_seq.count("G")
```




    6




```python
# To find the percentage of Gs and Cs:

100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    40.74074074074074




```python
# these calculations are already builin in biopython:
# Seq package
# Seq() makes something a seq object which allows us to perform biological applications

from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATACGATTGCATGCAGCATAAACGTA")
```


```python
gc_fraction(my_seq)
```




    0.4074074074074074




```python
# we can slice sequences into multiple parts:
# To cut out every 3rd nucleotide and choose our start position:

my_seq[0::3]
```




    Seq('GAAGTAAAG')




```python
my_seq[1::3]
```




    Seq('ACTCGGTAT')




```python
my_seq[2::3]
```




    Seq('TGTACCACA')




```python
# To print the string backwards:

my_seq[::-1]

```




    Seq('ATGCAAATACGACGTACGTTAGCATAG')




```python
# To change the object seq back into a string:

str(my_seq)
```




    'GATACGATTGCATGCAGCATAAACGTA'




```python
fasta_format_string = ">Name\n%s\n" % my_seq

# Here we use >Name as a placeholder to label our seq
# > is fasta format
```


```python
print(fast_format_string)
```

    >Name
    GATACGATTGCATGCAGCATAAACGTA
    



```python
# Concatenating (adding two scripts together)

seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
seq2 + seq1 
```




    Seq('AACCGGACGT')




```python
# manipulating strings

contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTTGCA")]
```


```python
# N is where sequencer was not confident in nucleotide id
# N can be ACT or G
spacer = Seq("N" * 10)
```


```python
# We take the spacer object and join it with contigs
# It puts 10 Ns to join each one

spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTTGCA')




```python
# case sensitivity issue:

dna_seq = Seq("acgtAcGT")
```


```python
# to fix:

dna_seq
```




    Seq('acgtAcGT')




```python
# to make our dna_seq all uppercase:

dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# to make it lowercase:

dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq = dna_seq.upper()
```


```python
# if we want to find a certain seq

"gtac" in dna_seq

# we get false because we asked for lowercase but our object is in uppercase
```




    False




```python
"GTAC" in dna_seq
```




    True




```python
# in biopython

my_seq = Seq("GATCCATTTGGCCATGCATGACCCGATCAAATTGC")
```


```python
# To get the complement to my_seq:

my_seq.complement()
```




    Seq('CTAGGTAAACCGGTACGTACTGGGCTAGTTTAACG')




```python
# reverse complement

my_seq.reverse_complement()
```




    Seq('GCAATTTGATCGGGTCATGCATGGCCAAATGGATC')




```python
# a protein example:

protein_seq = Seq("EVRNAK")
protein_seq.complement()

# there are ambiguity codes for nucleotides
```




    Seq('EBYNTM')




```python
# To create an object for coding dna and its template

coding_dna = Seq("ACGCAGTTCCAAGGTTAGCAGATGACGTAAATTGCAT")
```


```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('ATGCAATTTACGTCATCTGCTAACCTTGGAACTGCGT')




```python
coding_dna
```




    Seq('ACGCAGTTCCAAGGTTAGCAGATGACGTAAATTGCAT')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna

# we see that python took Ts and changed them to Us with transcribe
```




    Seq('ACGCAGUUCCAAGGUUAGCAGAUGACGUAAAUUGCAU')




```python
template_dna.reverse_complement().transcribe()

# gives same result bc template is reverse of coding_dna
```




    Seq('ACGCAGUUCCAAGGUUAGCAGAUGACGUAAAUUGCAU')




```python
# reverse transcription example

messenger_rna.back_transcribe()

# we took coding changed it to messenger rna and then back to coding with back_transcribe
```




    Seq('ACGCAGTTCCAAGGTTAGCAGATGACGTAAATTGCAT')




```python
messenger_rna
```




    Seq('ACGCAGUUCCAAGGUUAGCAGAUGACGUAAAUUGCAU')




```python
# to translate messenger_rna into amino acids

messenger_rna.translate()

# we get * for stop codons
# premature stop codon
```




    Seq('TQFQG*QMT*IA')




```python
# Seq objects pt 3

# continuing with translation and proteins

# we can specify the codon tables to determine if nuclear or mintochondrial genome codons
# to translate in mito genome

coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('TQFQG*QMT*IA')




```python
# to translate up until the first in frame stop codon, and then stop

coding_dna.translate(to_stop = True)
```




    Seq('TQFQG')




```python
# we get the same thing with this

coding_dna.translate(table=2, to_stop=True)
```




    Seq('TQFQG')




```python
# To change how stop codons are displayed:

coding_dna.translate(table=2, stop_symbol = "!")
```




    Seq('TQFQG!QMT!IA')




```python
# to handle a coding sequence with a non-standard start codon (like in bacteria) with a complete cds
# a long nucleotide seq:

gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTTCTGGTTCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
# Count individual nucleotides
count_A = gene.count("A")
count_T = gene.count("T")
count_G = gene.count("G")
count_C = gene.count("C")

print(f"Count of A: {count_A}")
print(f"Count of T: {count_T}")
print(f"Count of G: {count_G}")
print(f"Count of C: {count_C}")

total_N = count_A + count_T + count_G + count_C

print(total_N)

gene.translate(table = "Bacterial")
```

    Count of A: 77
    Count of T: 61
    Count of G: 74
    Count of C: 82
    294





    Seq('VKKMQSIVLALSLVSGSPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
# To show translation up until the stop codon:

gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVSGSPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
# To tell biopython that this is a complete cds and there will be a different start codon 
# than methionine

gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVSGSPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# To find the stop and start codons in the mito_table:

mito_table.stop_codons

```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
# Seq comparisons are difficult to determine if two sequences are truly equal

seq1 = Seq("ACGT")
```


```python
# ask the question two ways:

"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
# if length of seq is known but not the number of letters
# we can create an arugment with unknown nucleotides but 10 total
# this might be useful to create an object that we can fill later
# or if we need to just compare lengths of genes

unknown_seq = Seq(None, 10)

```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
# Seq Objects pt 4
```


```python
# a sequence from a seq alignment printout of various organisms for comparison
# human chromsome example:
# The first number is the start postion, then the seq, and legnth of chromosome

seq = Seq({117512683:"TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
# to find the seq from 1000-1020

seq[1000:1020]

# this seq is not defined in our sequence
```




    Seq(None, length=20)




```python
# to find the nucleotides at specified positions:

seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
# from a postion to the end:

seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACTG")
```


```python
undefined_seq = Seq(None, length = 10)
```


```python
# To put them together

seq + undefined_seq + seq
```




    Seq({0: 'ACTG', 14: 'ACTG'}, length=18)




```python
# important to not change the sequence data
# we are not able to change our data (my_seq) as seen in above error
# data is immutable in python

my_seq = Seq("AGTCAGGGGAACCCCTTTTTTTTAGCAGAGTCGAAGACTGAAC")
```


```python
# For example, this will not work:

my_seq[5] = "G"
```


    ---------------------------------------------------------------------------

    TypeError                                 Traceback (most recent call last)

    <ipython-input-13-50e848a55f5c> in <module>
    ----> 1 my_seq[5] = "G"
    

    TypeError: 'Seq' object does not support item assignment



```python
# to create a mutable seq

from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq[5] = "C"
```


```python
# we have changed the fifth position

mutable_seq
```




    MutableSeq('AGTCACGGGAACCCCTTTTTTTTAGCAGAGTCGAAGACTGAAC')




```python
# To remove the first T

mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('AGCACGGGAACCCCTTTTTTTAGCAGAGTCGAAGACTGAAC')




```python
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('CAAGTCAGAAGCTGAGACGATTTTTTTCCCCAAGGGCACGA')




```python
# to change gene and then make immutable again

new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('CAAGTCAGAAGCTGAGACGATTTTTTTCCCCAAGGGCACGA')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCAGACAGAGTTTTGAGAGTCAGAGTC"
```


```python
reverse_complement(my_string)
```




    'GACTCTGACTCTCAAAACTCTGTCTGC'




```python
transcribe(my_string)
```




    'GCAGACAGAGUUUUGAGAGUCAGAGUC'




```python
back_transcribe(my_string)
```




    'GCAGACAGAGTTTTGAGAGTCAGAGTC'




```python
translate(my_string)
```




    'ADRVLRVRV'




```
# Advanced Python - Biopython suite
# Sequence Annotation video pt. 1 
# In terminal: pip install biopython

from Bio.SeqRecord import SeqRecord
```


```python
# To find help file:
# help(SeqRecord)
```


```python
from Bio.Seq import Seq
```


```python
# first create a seq object:

simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
# We see that we created it, but it is blank
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# Let's add info:

simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequence for the VDB Comp Bio class"
```


```python
# print just the description:

print(simple_seq_r.description)
```

    Made up sequence for the VDB Comp Bio class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
# we can see that we added a seq, id, and description; we are storing sequences with annotation
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Comp Bio class', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"] = "None. This is just an example"
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Comp Bio class', dbxrefs=[])




```python
# scores for sequencing; letter annotations
# scores for how confident the algorithm is for each call of the sequencing

simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
# Let's use a real fasta file:
# https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna


```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# To pull pieces of information from the file:
# 3 examples below
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# we can pull an empty feature and see there is no annotation provided in file:
# 3 examples below
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
# pt 2 video
# to get genbank file for data:
# https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
# to see how many features there are

len(record.features)
```




    41




```python
# to add "fuzzy" info about genes (not verified yet)

from Bio import SeqFeature
```


```python
# example: if we think that the start of a gene is at 5 we create an object

start_pos = SeqFeature.AfterPosition(5)
```


```python
# If we think the end position is 
end_pos = SeqFeature.BetweenPosition(9, left = 8, right = 9)
```


```python

```


```python
# if we think a postion is one of something (ambiguous)

my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
# to simplify

int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
# to pass the numbers to the verbs

exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
# pt 3 video

from Bio.SeqRecord import SeqRecord
```


```python
# to create a seq record for a protein (AA seq)

record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSVITAVTFRGPSETHLDSMVGQALFGD"
                       "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
                       "NIEKSLKEAFTPLGISDWINSTFWIAHPGGPAILDQVEAKLGLEEKMRATREVLSEYGNM"
                       "SSAC"),
                                id="gi|14150838|gb|AAK54648.1|AF376133_1",
                   description="chalcone synthase [Cucumis sativus]",
                  )
```


```python
# To print as a FASTA file:

print(record.format("fasta"))
```

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSVITAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWINSTFWIAHPGGPAILDQVEAKLGLEEKMRATREVLSEYGNM
    SSAC
    



```python
print(record)
# note that we did not add any annotated features
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSVITAVTFRGPSETHLDSMVG...SAC')



```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# To see how many basepairs:

len(record)
```




    9609




```python
# To see the number of features

len(record.features)
```




    41




```python
print(record.features[20])
# wee can see the 20th feature
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# To subdivide the record (focus on a smaller part of the sequence)

sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
len(sub_record.features)
```




    2




```python
# To see the first feature in the sub_record
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_record.dbxrefs
```




    []




```python
# To add to our annotations:
# we have a linear piece of DNA

sub_record.annotations["topology"] = "linear"
```


```python
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
sub_record.id
```




    'NC_005816.1'




```python
sub_record.name
```




    'NC_005816'




```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# To update the description feature

sub_record.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
```


```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
# To see everything thru 200 and ... as genbank file:

print(sub_record.format("genbank")[:200] +"...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...



```python
# pt 4
# A circular genome:

record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# to shift the genome to change where we start reading the genome:
# because it's circular, where to start isn't specified

shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# length hasn't been changed

len(shifted)
```




    9609




```python
# one less feature

len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
# This feature was not preserved

shifted.dbxrefs
```




    []




```python
# To explicitly 
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# reverse complement
# %s means to print the first value of a string
# %i means to present the next four as integers

print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
# if you run reverse complement on this sequence, are the annotations conserved?

rc = record.reverse_complement(id = "Testing")
```


```python
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
# length and features are the same, but refs and annotations are different (by default to make sure I am not creating duplicate ids with sequences I've manipulated)
```

    Testing 9609 41 0 0



```python

```


```python
# Blast Challenge
# To explore how closely "x gene" in humans is related to other organisms
# The lower the E-value, the more significant and less likely the similarity occured by chance.
# Human Cancer Gene MYCL (4610)
# After installing biopython in terminal, import the following:

from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML

```


```python
NCBIWWW.email = "melissabarkemeyer@gmail.com"
```


```python
# Obtain the FASTA file from NCBI for gene of interest
# https://www.ncbi.nlm.nih.gov/refseq/
# to blast the uploaded FASTA file, use SeqIO to read it
# Because my file contains more than one FASTA sequence, I need to parse it.

record = next(SeqIO.parse("gene_1.fasta", "fasta"))

```


```python
# Let's view the file
print(record)
```

    ID: NC_000001.11:c39901917-39895428
    Name: NC_000001.11:c39901917-39895428
    Description: NC_000001.11:c39901917-39895428 MYCL [organism=Homo sapiens] [GeneID=4610] [chromosome=1]
    Number of features: 0
    Seq('GAGTGCGGGCCGCGCTCTCGGCGGCGCGCATGTGCGTGTGTGCTGGCTGCCGGG...ACA')



```python
# Now we can blast using the uploaded sequence, which we saved as 'record':

# Submit BLAST request to NCBI 
result_handle = NCBIWWW.qblast(
    program="blastn",
    database="nt",
    sequence=str(record.seq),  
    format_type="XML"           
)


```


```python
# Now to read the result handle results (XML text), write it to 
# the created XML file, then close the connection to NCBI search
# Note: The FASTA file is still FASTA, not overwritten
with open("gene_1_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

result_handle.close()

```



```python
# To create an object from the XML file

with open("gene_1_blast.xml") as result_handle:
    blast_record = NCBIXML.read(result_handle)
```


```python

# To find the best hit (alignment) with chimpanzee:

best_hsp = None

for alignment in blast_record.alignments:
   if "Pan troglodytes" in alignment.title or "chimpanzee" in alignment.title or "Pan" in alignment.title:
        for hsp in alignment.hsps:
            if best_hsp is None or hsp.expect < best_hsp.expect:
                best_hsp = hsp
                best_title = alignment.title

if best_hsp:
    print("Best chimpanzee hit:")
    print("Title:", best_title)
    print("E-value:", best_hsp.expect)
    print(best_hsp.query)
    print(best_hsp.match)
    print(best_hsp.sbjct)
else:
    print("No chimpanzee hits found.")


```

    Best chimpanzee hit:
    Title: gi|2697761507|ref|XM_016959814.4| PREDICTED: Pan troglodytes MYCL proto-oncogene, bHLH transcription factor (MYCL), mRNA
    E-value: 0.0
    GAGAATGAAGAAATTGATGTTGTGACAGTAGAGAAGAGGCAGTCTCTGGGTATTCGGAAGCCGGTCACCATCACGGTGCGAGCAGACCCCCTGGATCCCTGCATGAAGCATTTCCACATCTCCATCCATCAGCAACAGCACAACTATGCTGCCCGTTTTCCTCCAGAAAGCTGCTCCCAAGAAGAGGCTTCAGAGAGGGGTCCCCAAGAAGAGGTTCTGGAGAGAGATGCTGCAGGGGAAAAGGAAGATGAGGAGGATGAAGAGATTGTGAGTCCCCCACCTGTAGAAAGTGAGGCTGCCCAGTCCTGCCACCCCAAACCTGTCAGTTCTGATACTGAGGATGTGACCAAGAGGAAGAATCACAACTTCCTGGAGCGCAAGAGGCGGAATGACCTGCGTTCGCGATTCTTGGCGCTGAGGGACCAGGTGCCCACCCTGGCCAGCTGCTCCAAGGCCCCCAAAGTAGTGATCCTAAGCAAGGCCTTGGAATACTTGCAAGCCCTGGTGGGGGCTGAGAAGAGGATGGCTACAGAGAAAAGACAGCTCCGATGCCGGCAGCAGCAGTTGCAGAAAAGAATTGCATACCTCACTGGCTACTAACTGACCAAAAAGCCTGACAGTTCTGTCTTACGAAGACACAAGTTTATTTTTTAACCTCCCTCTCCCCTTTAGTAATTTGCACATTTTGGTTATGGTGGGACAGTCTGGACAGTAGATCCCAGAATGCATTGCAGCCGGTGCACACACAATAAAGGCTTGCATTCTTGGAAACCTTGAAACCCAGCTCTCCCTCTTCCCTGACTCATGGGAGTGCTGTATGTTCTCTGGCGCCTTTGGCTTCCCAGCAGGCAGCTGACTGAGGAGCCTTGGGGTCTGCCTAGCTCACTAGCTCTGAAGAAAAGGCTGACAGATGCTATGCAACAGGTGGTGGATGTTGTCAGGGGCTCCAGCCTGCATGAAATCTCACACTCTGCATGAGCTTTAGGCTAGGAAAGGATGCTCCCAACTGGTGTCTCTGGGGTGATGCAAGGACAGCTGGGCCTGGATGCTCTCCCTGAGGCTCCTTTTTCCAGAAGACACACGAGCTGTCTTGGGTGAAGACAAGCTTGCAGACTTGATCAACATTGACCATTACCTCACTGTCAGACACTTTACAGTAGCCAAGGAGTTGGAAACCTTTATATATTATGATGTTAGCTGACCCCCTTCCTCCCACTCCCAATGCTGCGACCCTGGGAACACTTAAAAAGCTTGGCCTCTAGATTCTTTGTCTCAGAGCCCTCTGGGCTCTCTCCTCTGAGGGAGGGACCTTTCTTTCCTCACAAGGGACTTTTTTGTTCCATTATGCCTTGTTATGCAATGGGCTCTACAGCACCCTTTCCCACAGGTCAGAAATATTTCCCCAAGACACAGGGAAATCGGTCCTAGCCTGGGGCCTGGGGATAGCTTGGAGTCCTGGCCCATGAACTTGATCCCTGCCCAGGTGTTTTCCGAGGGGCACTTGAGGCCCAGTCTTTTCTCAAGGCAGGTGTAAGACACCTCAGAGGGAGAACTGTACTGCTGCCTCTTTCCCACCTGCCTCATCTCAATCCTTGAGCGGCAAGTTTGAAGTTCTTCTGGAACCATGCAAATCTGTCCTCCTCATGCAATTCCAAGGAGCTTGCTGGCTCTGCAGCCACCCTTGGGCCCCTTCCAGCCTGCCATGAATCAGATATCTTTCCCAGAATCTGGGCGTTTCTGAAGTTTTGGGGAGAGCTGTTGGGACTCATCCAGTGCTCCAGAAGGTGGACTTGCTTCTGGTGGGTTTTAAAGGAGCCTCCAGGAGATATGCTTAGCCAACCATGATGGATTTTACCCCAGCTGGACTCGGCAGCTCCAAGTGGAATCCACGTGCAGCTTCTAGTCTGGGAAAGTCACCCAACCTAGCAGTTGTCATGTGGGTAACCTCAGGCACCTCTAAGCCTGTCCTGGAAGAAGGACCAGCAGCCCCTCCAGAACTCTGCCCAGGACAGCAGGTGCCTGCTGGCTCTGGGTTTGGAAGTTGGGGTGGGTAGGGGGTGGTAAGTACTATATATGGCTCTGGAAAACCAGCTGCTACTTCCAAATCTATTGTCCATAATGGTTTCTTTCTGAGGTTGCTTCTTGGCCTCAGAGGACCCCAGGGGATGTTTGGAAATAGCCTCTCTACCCTTCTGGAGCATGGTTTACAAAAGCCAGCTGACTTCTGGAATTGTCTATGGAGGACAGTTTGGGTGTAGGTTACTGATGTCTCAACTGAATAGCTTGTGTTTTATAAGCTGCTGTTGGCTATTATGCT-GGGGGAGTCTTTTTTTTTTATATTGTATTTTTGTATGCCTTTTGCAAAGTGGTGTTAACTGTTTTTGTACAAGGAAAAAAACTCTTGGGGCAATT-TCCTGTTGCAAGGGTCTGATTTATTTTGAAAGGCAAGTTCACCTGAAATTTTGTATTTAGTTGTGATTACTGATTGCCTGATTTTAAAATGTTGCCTTCTGGGACATCTTCTAATAAAAGATTTCTCAAACA
    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    GAGAATGAAGAAATTGATGTTGTGACAGTAGAGAAGAGGCAGTCTCTGGGTATTCGGAAGCCGGTCACCATCACGGTGCGAGCAGACCCCCTGGATCCCTGCATGAAGCATTTCCACATCTCCATCCATCAGCAACAGCACAACTATGCTGCCCGTTTTCCTCCAGAAAGCTGCTCCCAAGAAGAGGCTTCAGAGAGGGGTCCCCAAGAAGAGGTTCTGGAGAGAGATGCTGCAGGGGAAAAGGAAGATGAGGAGGATGAAGAGATTGTGAGTCCCCCACCTGTAGAAAGTGAGGCTGCCCAGTCCTGCCACCCCAAACCTGTCAGTTCTGATACTGAGGATGTGACCAAGAGGAAGAATCACAACTTCCTGGAGCGCAAGAGGCGGAATGACCTGCGTTCGCGATTCTTGGCGCTGAGGGACCAGGTGCCCACCCTGGCCAGCTGCTCCAAGGCCCCCAAAGTAGTGATCCTAAGCAAGGCCTTGGAATACTTGCAAGCCCTGGTGGGGGCTGAGAAGAGGATGGCTACAGAGAAAAGACAGCTCCGATGCCGGCAGCAGCAGTTGCAGAAAAGAATTGCATACCTCAGTGGCTACTAACTGACCAAAAAGCCTGACAGTTCTGTCTTACAAAGACACAAGTTTATTTTTTAACCTCCCTCTCCCCTTTAGTAATTTGCACATTTTGGTTATGGTGGGACAGTCTGGACAGTAGATCCCAGAATGCATTGCAGCCGGTGCACACACAATAAAGGCTTGCATTCTTGGAAACCTTGAAACCCAGCTCTCCCTCTTCCCTGACTCATGGGAGTGCTGTATGTTCTCTGGCGCCTTTGGCTTCCCAGCAGGCAGCTGACTGAGGAGCCTTGGGGTCTGCCTAGCTCACTAGCTCCGAAGAAAAGGCTGACAGATGCTATGCAACAGGTGGTGGATGTTGTCAGGGGCTCCAGCCTGCTTGAAATCTCACACTCTGCATGAGCTTTAGGCTAGGAAAGGATGCTCCCAACTGGTGTCTCTGGGGTGATGCAAGGACAGCTGGGCCTGGATGCTCTCCCTGAGGCTCCTTTTTCCAGAAGACACACGAGCTGTCTTGGGTGAAGACAAGCTTGCAGACTTGATCAACATTGACCATTACCTCACTGTCAGACACTTTACAGTAGCCAAGGAGTTGGAAACCTTTATATATTATGATGTTAGCTGA-CCCCTTCCTCCCACTCCCAATGCTGCGACCCTGGGAACACTTAAAAAGCTTGGCCTCTAGATTCTTTGTCTCAGAGCCCTCTGGGCTCTCTCCTCTGAGGGAGGGACCTTTCTTTCCTCACAAGGGACTTTTTTGTTCCATTATGCCTTGTTATGCAATGGGCTCTACAGCACCCTTTCTCACAGGTCAAAAATATTTCCCCAAGACACAGGGAAATCGGTCCTAGCCTGGGGCCTGGGGATAGCTTGGAGTCCTGGCCCATGAACTTGATCCCTGCCCAGGTGTTTTCCGAGGGGCACTTGAGGCCCAGTCTTTTCTCAAGGCAGGTGTAAGACACCTCAGAGGGAGAACTGTACTGCTGCCTCTTTCCCACCTGCCTCATCTCAATCCTTGAGCGGCAAGTTTGAAGTTCTTCTGGAACCATGCAAATCTGTCCTCCTCATGCAATTCCAAGGAGCTTGCTGGCTCTGCAGCCACCCTTGGGCCCCTTCCAGCCTGCCATGAATCAGATATCTTTCCCAGAATCTGGGCGTTTCTGAAGTTTTGGGGAGAGCTGTTGGGACTCATCCAGTGCTCCAGAAGGTGGACTTGCTTCTGGTGGGTTTTAAAGGAGCCTCCAGGAGATATGCTTAGCCAACCATGATGGATTTTACCCCAGCTGGACTCGGCAGCTCCAAGTGGAATCCACGTGCAGCCTCTAGTCTGGGAAAGTCACCCAACCTAGCAGTTGTCATGTGGGTAACCTCAGGCACCTCTAAGCCTGTCCTGGAAGAAGGACCAGCAGCCCCTCCAGAACTCTGCCCAGGACAGCAGGTGCCTGCTGGCTCTGGGTTTGGAAGTTGGGGTGGGTAGGGGGTGGTAAGTACTATATATGGCTCTGGAAAACCAGCTGCTACTTCCAAATCTATTGTCCATAATGGTTTCTTTCTGAGGTTGCTTCTTGGCCTCAGAGGACCCCAGGGGATGTTTGGAAATAGCCTCTCTACCCTTCTGGAGCATGGTTTACAAAAGCCAGCTGACTTCTGGAATTGTCTATGGAGGACAGTTTGGGTGTAGGTTACTGATGTCTCAACTGAATAGCTTGTGTTTTATAAGCTGCTGTTGGCTATTATGCTGGGGGGAGTCTTTTTTTTTTATATTGTATTTTTGTATGCCTTTTGCAAAGTGGTGTTAACTGTTTTTGTACAAGGAAAAAAACTCTT-GGGCAATTCTCCTGTTGCAAGGGTCTGATTTATTTTGAAAGGCAAGTTCACCTGAAATTTTGTATTTAGTTGTGATTACTGATTGCCTGATTTTAAAATGTTGCCTTCTGGGACATCTTCTAATAAAAGATTTCTCAAACA



```python

```
