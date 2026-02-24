# Barkemeyer_Python_Portfolio_2
This is the portfolio of python code that I learned in my Advanced Python course.
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
# We can access the different positions in seq objects:
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
# We can do a .count (does not count overlapping):

Seq("AAAA").count("AA")
```




    2




```python
# Another example:
my_seq = "GATACGATTGCATGCAGCATAAACGTA"
```


```python
len(my_seq)
```




    27




```python
# To find the count of Gs in the sequence:
my_seq.count("G")
```




    6




```python
# To find the percentage of Gs and Cs:

100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    40.74074074074074




```python
# These calculations are already built in biopython:
# Seq package
# Seq() makes something a seq object which allows us to perform biological applications.

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
# We can slice sequences into multiple parts.
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

# Here we use >Name as a placeholder to label our seq:
# > is fasta format
```


```python
print(fast_format_string)
```

    >Name
    GATACGATTGCATGCAGCATAAACGTA
    



```python
# Concatenating (adding two scripts together):

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
# Manipulating Strings

contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTTGCA")]
```


```python
# N is where sequencer was not confident in nucleotide id.
# N can be ACT or G.
spacer = Seq("N" * 10)
```


```python
# We take the spacer object and join it with contigs.
# It puts 10 Ns to join each one.

spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTTGCA')




```python
# Case sensitivity issue:

dna_seq = Seq("acgtAcGT")
```


```python
# To fix:

dna_seq
```




    Seq('acgtAcGT')




```python
# To make our dna_seq all uppercase:

dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# To make it lowercase:

dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq = dna_seq.upper()
```


```python
# If we want to find a certain seq:

"gtac" in dna_seq

# We get false because we asked for lowercase but our object is in uppercase.
```




    False




```python
"GTAC" in dna_seq
```




    True




```python
# In biopython

my_seq = Seq("GATCCATTTGGCCATGCATGACCCGATCAAATTGC")
```


```python
# To get the complement to my_seq:

my_seq.complement()
```




    Seq('CTAGGTAAACCGGTACGTACTGGGCTAGTTTAACG')




```python
# Reverse Complement

my_seq.reverse_complement()
```




    Seq('GCAATTTGATCGGGTCATGCATGGCCAAATGGATC')




```python
# A protein example:

protein_seq = Seq("EVRNAK")
protein_seq.complement()

```




    Seq('EBYNTM')




```python
# To create an object for coding DNA and its template:

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

# We see that python took Ts and changed them to Us with 'transcribe'.
```




    Seq('ACGCAGUUCCAAGGUUAGCAGAUGACGUAAAUUGCAU')




```python
template_dna.reverse_complement().transcribe()

# Gives the same result bc template is reverse of coding_dna.
```




    Seq('ACGCAGUUCCAAGGUUAGCAGAUGACGUAAAUUGCAU')




```python
# Reverse transcription example

messenger_rna.back_transcribe()

# We took coding changed it to messenger RNA and then back to coding with back_transcribe
```




    Seq('ACGCAGTTCCAAGGTTAGCAGATGACGTAAATTGCAT')




```python
messenger_rna
```




    Seq('ACGCAGUUCCAAGGUUAGCAGAUGACGUAAAUUGCAU')




```python
# To translate messenger_rna into amino acids:

messenger_rna.translate()

# We get * for stop codons
# premature stop codon
```




    Seq('TQFQG*QMT*IA')




```python
# Seq Objects Pt. 3
# Continuing with translation and proteins
# We can specify the codon tables to determine if nuclear or mintochondrial genome codons.
# To translate in mitochondrial genome codons:

coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('TQFQG*QMT*IA')




```python
# To translate up until the first in frame stop codon, and then stop:

coding_dna.translate(to_stop = True)
```




    Seq('TQFQG')




```python
# We get the same thing with this

coding_dna.translate(table=2, to_stop=True)
```




    Seq('TQFQG')




```python
# To change how stop codons are displayed:

coding_dna.translate(table=2, stop_symbol = "!")
```




    Seq('TQFQG!QMT!IA')




```python
# To handle a coding sequence with a non-standard start codon (like in bacteria) with a complete cds:
# A long nucleotide seq:

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
# than methionine:

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
# Seq comparisons are difficult to determine if two sequences are truly equal.

seq1 = Seq("ACGT")
```


```python
# You can ask the question two ways:

"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
# If length of seq is known but not the number of letters,
# we can create an arugment with unknown nucleotides but 10 total.
# This might be useful to create an object that we can fill later
# or if we need to just compare lengths of genes.

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
# Seq Objects Pt. 4
```


```python
# A sequence from a seq alignment printout of various organisms for comparison
# Human chromsome example
# The first number is the start postion, then the seq, and legnth of chromosome.

seq = Seq({117512683:"TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
# To find the seq from 1000-1020:

seq[1000:1020]

# This seq is not defined in our sequence:
```




    Seq(None, length=20)




```python
# To find the nucleotides at specified positions:

seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
# From a postion to the end:

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
# To put them together:

seq + undefined_seq + seq
```




    Seq({0: 'ACTG', 14: 'ACTG'}, length=18)




```python
# It is important to not change the sequence data.
# We are not able to change our data (my_seq) as seen in above error.
# Data is immutable in python.

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
# To create a mutable seq:

from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq[5] = "C"
```


```python
# We have changed the fifth position.

mutable_seq
```




    MutableSeq('AGTCACGGGAACCCCTTTTTTTTAGCAGAGTCGAAGACTGAAC')




```python
# To remove the first T:

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
# To change the gene sequence and then make immutable again:

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
# Biopython Suite
# Sequence Annotation Pt. 1 
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
# First, create a seq object:

simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
# We see that we created it, but it is blank.
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
# We can see that we added a seq, id, and description; we are storing sequences with annotation.
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Comp Bio class', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"] = "None. This is just an example."
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example.



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Comp Bio class', dbxrefs=[])




```python
# Scores for sequencing; letter annotations
# Scores for how confident the algorithm is for each call of the sequencing

simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
# Let's use a real FASTA file:
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
# 3 examples below:
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
# We can pull an empty feature and see there is no annotation provided in file:
# 3 examples below:
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
# Sequence Annotations Pt. 2 
# To get genbank file for data:
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
# To see how many features there are:

len(record.features)
```




    41




```python
# To add "fuzzy" info about genes (not verified yet):

from Bio import SeqFeature
```


```python
# Example: if we think that the start of a gene is at 5, we create an object:

start_pos = SeqFeature.AfterPosition(5)
```


```python
# If we think we know what the end position is: 
end_pos = SeqFeature.BetweenPosition(9, left = 8, right = 9)
```


```python

```


```python
# If we think a postion is one of something (ambiguous):

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
# To simplify:

int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
# To pass the numbers to the verbs:

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
# Sequence Annotations Pt. 3 

from Bio.SeqRecord import SeqRecord
```


```python
# To create a seq record for a protein (AA seq):

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
# Note that we did not add any annotated features.
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
# To see the number of features:

len(record.features)
```




    41




```python
print(record.features[20])
# We can see the 20th feature.
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
# To subdivide the record (focus on a smaller part of the sequence):

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
# To see the first feature in the sub_record:
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
# Example: We have a linear piece of DNA.

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
# To update the description feature:

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
# Sequence Annotations Pt. 4
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
# To shift the genome to change the start location for reading the genome:
# Because it's circular, where to start isn't specified.

shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# The length hasn't been changed.

len(shifted)
```




    9609




```python
# One less feature:

len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
# This feature was not preserved:

shifted.dbxrefs
```




    []




```python
# To explicitly shift:
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
# Reverse Complement
# %s means to print the first value of a string
# %i means to present the next four as integers

print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
# If you run reverse complement on this sequence, are the annotations conserved?

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



```

```python
# Seq Input/Ouput Pt. 1
# go to help(SeqIO) for help with code
# https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk
# https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.fasta
from Bio import SeqIO
```


```python
# We will manipulate a file that contains genes of many organisms
# we will use parse - an iterator that will go through the file record
# load file
# print id, seq, and length

for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
                          
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# to extract things as list
# go thru the whole file and pull out the identifier number 
# useful for looking for SNPs

identifiers = {seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")}
```


```python
identifiers
```




    {'Z78439.1',
     'Z78440.1',
     'Z78441.1',
     'Z78442.1',
     'Z78443.1',
     'Z78444.1',
     'Z78445.1',
     'Z78446.1',
     'Z78447.1',
     'Z78448.1',
     'Z78449.1',
     'Z78450.1',
     'Z78451.1',
     'Z78452.1',
     'Z78453.1',
     'Z78454.1',
     'Z78455.1',
     'Z78456.1',
     'Z78457.1',
     'Z78458.1',
     'Z78459.1',
     'Z78460.1',
     'Z78461.1',
     'Z78462.1',
     'Z78463.1',
     'Z78464.1',
     'Z78465.1',
     'Z78466.1',
     'Z78467.1',
     'Z78468.1',
     'Z78469.1',
     'Z78470.1',
     'Z78471.1',
     'Z78472.1',
     'Z78473.1',
     'Z78474.1',
     'Z78475.1',
     'Z78476.1',
     'Z78477.1',
     'Z78478.1',
     'Z78479.1',
     'Z78480.1',
     'Z78481.1',
     'Z78482.1',
     'Z78483.1',
     'Z78484.1',
     'Z78485.1',
     'Z78486.1',
     'Z78487.1',
     'Z78488.1',
     'Z78489.1',
     'Z78490.1',
     'Z78491.1',
     'Z78492.1',
     'Z78493.1',
     'Z78494.1',
     'Z78495.1',
     'Z78496.1',
     'Z78497.1',
     'Z78498.1',
     'Z78499.1',
     'Z78500.1',
     'Z78501.1',
     'Z78502.1',
     'Z78503.1',
     'Z78504.1',
     'Z78505.1',
     'Z78506.1',
     'Z78507.1',
     'Z78508.1',
     'Z78509.1',
     'Z78510.1',
     'Z78511.1',
     'Z78512.1',
     'Z78513.1',
     'Z78514.1',
     'Z78515.1',
     'Z78516.1',
     'Z78517.1',
     'Z78518.1',
     'Z78519.1',
     'Z78520.1',
     'Z78521.1',
     'Z78522.1',
     'Z78523.1',
     'Z78524.1',
     'Z78525.1',
     'Z78526.1',
     'Z78527.1',
     'Z78529.1',
     'Z78530.1',
     'Z78531.1',
     'Z78532.1',
     'Z78533.1'}




```python
# another function: next
# python remembers and keeps calling the next record in line
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
second_record = next(record_iterator)
```


```python
print(second_record.id)
```

    gi|2765657|emb|Z78532.1|CCZ78532



```python
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
# if you want to jump to a record without using next

records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
print("Found %i records" % len(records))
```

    Found 94 records



```python
print("The last record")
last_record = records[-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("The first record")
first_record = records[0]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740



```python
# Sequence Input/Output pt 2
# extracting data from seq files
from Bio import SeqIO
```


```python
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
all_species = []
```


```python
# for the file, append organism:

for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
# to see all 94 species being compared in the file
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
# another way to do it:

all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank" )
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
# We get genus abbreviated then species b/c fasta are arranged differently

print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
# to alter data

record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
# we can see part of the first line of the record:

first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
# to change the frst_record.id for helpful descriptions:

first_record.id = "new_id"
```


```python
first_record.id
```




    'new_id'




```python
first_record.description = first_record.id + " " + "desired new description"
```


```python
# to view the change:

print(first_record.format("fasta"[:200]))
```

    >new_id desired new description
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    



```python
# Sequence Input/Output pt 3

# first, write a fasta file

from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord(
    Seq("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
),
    id="sghkajghkjshdgkjasndgkjadgn_1",
    description  = "ubiquitin [homo sapien]",
)
```


```python
rec2 = SeqRecord(
    Seq("HSQGTFTSDYSKYLDSRRAQDFVQWLMNT",
),
    id = "alsdkfhdlkfjhakjsdhfkjh_2",
        description = "glucagon [homo sapien]",
        )
```


```python
rec3 = SeqRecord(
    Seq("CGNLSTCMLGTYTQDFNKFHTFPQTAIGVGAP",
),
    id = "alsdkfhdlkfjhakjeieieieh_3",
        description = "calcitonin [homo sapien]",
        )
```


```python
my_records = [rec1, rec2, rec3]
```


```python
# To create a new file using the records created:

SeqIO.write(my_records, "my_example.faa", "fasta")
```




    3




```python
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
# To convert genbank to a fasta file:
count = SeqIO.write(records, "my_example.fasta", "fasta")
# to check to make sure that all 94 records were converted:
print("Converted %i records" % count)
```

    Converted 94 records



```python
# to print 94 reverse complements of each seq (not applied to id)
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement)
```

    Z78533.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')>
    Z78532.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')>
    Z78531.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')>
    Z78530.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')>
    Z78529.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')>
    Z78527.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')>
    Z78526.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')>
    Z78525.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')>
    Z78524.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')>
    Z78523.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')>
    Z78522.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')>
    Z78521.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')>
    Z78520.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')>
    Z78519.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')>
    Z78518.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')>
    Z78517.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')>
    Z78516.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')>
    Z78515.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')>
    Z78514.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')>
    Z78513.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')>
    Z78512.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')>
    Z78511.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')>
    Z78510.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')>
    Z78509.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')>
    Z78508.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')>
    Z78507.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')>
    Z78506.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')>
    Z78505.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')>
    Z78504.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')>
    Z78503.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')>
    Z78502.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')>
    Z78501.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')>
    Z78500.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')>
    Z78499.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')>
    Z78498.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')>
    Z78497.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78496.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')>
    Z78495.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')>
    Z78494.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')>
    Z78493.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')>
    Z78492.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')>
    Z78491.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')>
    Z78490.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78489.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')>
    Z78488.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')>
    Z78487.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')>
    Z78486.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')>
    Z78485.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')>
    Z78484.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')>
    Z78483.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')>
    Z78482.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')>
    Z78481.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78480.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78479.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')>
    Z78478.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')>
    Z78477.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')>
    Z78476.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')>
    Z78475.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')>
    Z78474.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')>
    Z78473.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')>
    Z78472.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78471.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78470.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')>
    Z78469.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')>
    Z78468.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')>
    Z78467.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78466.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')>
    Z78465.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')>
    Z78464.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')>
    Z78463.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')>
    Z78462.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')>
    Z78461.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')>
    Z78460.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')>
    Z78459.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')>
    Z78458.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')>
    Z78457.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')>
    Z78456.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78455.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')>
    Z78454.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')>
    Z78453.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')>
    Z78452.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')>
    Z78451.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78450.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')>
    Z78449.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')>
    Z78448.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')>
    Z78447.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')>
    Z78446.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')>
    Z78445.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')>
    Z78444.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')>
    Z78443.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')>
    Z78442.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')>
    Z78441.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')>
    Z78440.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')>
    Z78439.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')>



```python
# To specify that this is a reverse complement by adding rc at the beginning of each id and reverse complement as the desc:
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
]
```


```python
print(records)
```

    [SeqRecord(seq=Seq('GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGT...ACG'), id='rc_gi|2765658|emb|Z78533.1|CIZ78533', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAAT...ACG'), id='rc_gi|2765657|emb|Z78532.1|CCZ78532', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCC...ACG'), id='rc_gi|2765656|emb|Z78531.1|CFZ78531', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGG...ACG'), id='rc_gi|2765655|emb|Z78530.1|CMZ78530', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAAC...CGT'), id='rc_gi|2765654|emb|Z78529.1|CLZ78529', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGT...ACG'), id='rc_gi|2765652|emb|Z78527.1|CYZ78527', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAAT...ACG'), id='rc_gi|2765651|emb|Z78526.1|CGZ78526', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATC...ACA'), id='rc_gi|2765650|emb|Z78525.1|CAZ78525', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAAC...ACG'), id='rc_gi|2765649|emb|Z78524.1|CFZ78524', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCA...ACG'), id='rc_gi|2765648|emb|Z78523.1|CHZ78523', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCC...ACG'), id='rc_gi|2765647|emb|Z78522.1|CMZ78522', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAG...TAC'), id='rc_gi|2765646|emb|Z78521.1|CCZ78521', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACC...ACG'), id='rc_gi|2765645|emb|Z78520.1|CSZ78520', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAA...TAT'), id='rc_gi|2765644|emb|Z78519.1|CPZ78519', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATC...ACG'), id='rc_gi|2765643|emb|Z78518.1|CRZ78518', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAAT...ACG'), id='rc_gi|2765642|emb|Z78517.1|CFZ78517', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCA...ACG'), id='rc_gi|2765641|emb|Z78516.1|CPZ78516', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAAT...ACG'), id='rc_gi|2765640|emb|Z78515.1|MXZ78515', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCC...ACG'), id='rc_gi|2765639|emb|Z78514.1|PSZ78514', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGG...ACG'), id='rc_gi|2765638|emb|Z78513.1|PBZ78513', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGT...ACG'), id='rc_gi|2765637|emb|Z78512.1|PWZ78512', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGA...ACG'), id='rc_gi|2765636|emb|Z78511.1|PEZ78511', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAAC...TAG'), id='rc_gi|2765635|emb|Z78510.1|PCZ78510', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACC...ACG'), id='rc_gi|2765634|emb|Z78509.1|PPZ78509', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGG...ACG'), id='rc_gi|2765633|emb|Z78508.1|PLZ78508', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGG...ACG'), id='rc_gi|2765632|emb|Z78507.1|PLZ78507', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCC...ACG'), id='rc_gi|2765631|emb|Z78506.1|PLZ78506', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAAC...ACG'), id='rc_gi|2765630|emb|Z78505.1|PSZ78505', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCA...ACG'), id='rc_gi|2765629|emb|Z78504.1|PKZ78504', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAAC...ACG'), id='rc_gi|2765628|emb|Z78503.1|PCZ78503', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGT...ACG'), id='rc_gi|2765627|emb|Z78502.1|PBZ78502', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGT...ACG'), id='rc_gi|2765626|emb|Z78501.1|PCZ78501', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTC...ACG'), id='rc_gi|2765625|emb|Z78500.1|PWZ78500', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGT...ACG'), id='rc_gi|2765624|emb|Z78499.1|PMZ78499', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAAC...ACG'), id='rc_gi|2765623|emb|Z78498.1|PMZ78498', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765622|emb|Z78497.1|PDZ78497', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCA...ACG'), id='rc_gi|2765621|emb|Z78496.1|PAZ78496', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGA...ACG'), id='rc_gi|2765620|emb|Z78495.1|PEZ78495', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATT...ACG'), id='rc_gi|2765619|emb|Z78494.1|PNZ78494', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATG...ACG'), id='rc_gi|2765618|emb|Z78493.1|PGZ78493', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAAT...ACG'), id='rc_gi|2765617|emb|Z78492.1|PBZ78492', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAA...ACG'), id='rc_gi|2765616|emb|Z78491.1|PCZ78491', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCA...ACG'), id='rc_gi|2765615|emb|Z78490.1|PFZ78490', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765614|emb|Z78489.1|PDZ78489', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGG...CAG'), id='rc_gi|2765613|emb|Z78488.1|PTZ78488', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATC...ACG'), id='rc_gi|2765612|emb|Z78487.1|PHZ78487', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGAT...ACG'), id='rc_gi|2765611|emb|Z78486.1|PBZ78486', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACT...CAG'), id='rc_gi|2765610|emb|Z78485.1|PHZ78485', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCA...ACG'), id='rc_gi|2765609|emb|Z78484.1|PCZ78484', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_gi|2765608|emb|Z78483.1|PVZ78483', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...AGA'), id='rc_gi|2765607|emb|Z78482.1|PEZ78482', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGAT...ACG'), id='rc_gi|2765606|emb|Z78481.1|PIZ78481', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGAT...ACG'), id='rc_gi|2765605|emb|Z78480.1|PGZ78480', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTG...ACG'), id='rc_gi|2765604|emb|Z78479.1|PPZ78479', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAA...ACG'), id='rc_gi|2765603|emb|Z78478.1|PVZ78478', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTC...ACG'), id='rc_gi|2765602|emb|Z78477.1|PVZ78477', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGG...ACG'), id='rc_gi|2765601|emb|Z78476.1|PGZ78476', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCC...ACG'), id='rc_gi|2765600|emb|Z78475.1|PSZ78475', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAAC...ACG'), id='rc_gi|2765599|emb|Z78474.1|PKZ78474', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAAT...ACG'), id='rc_gi|2765598|emb|Z78473.1|PSZ78473', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_gi|2765597|emb|Z78472.1|PLZ78472', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765596|emb|Z78471.1|PDZ78471', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTAT...ACG'), id='rc_gi|2765595|emb|Z78470.1|PPZ78470', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTA...ACG'), id='rc_gi|2765594|emb|Z78469.1|PHZ78469', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAA...ACG'), id='rc_gi|2765593|emb|Z78468.1|PAZ78468', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATC...ACG'), id='rc_gi|2765592|emb|Z78467.1|PSZ78467', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGG...ACG'), id='rc_gi|2765591|emb|Z78466.1|PPZ78466', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCA...ACG'), id='rc_gi|2765590|emb|Z78465.1|PRZ78465', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCA...ACG'), id='rc_gi|2765589|emb|Z78464.1|PGZ78464', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTC...ACG'), id='rc_gi|2765588|emb|Z78463.1|PGZ78463', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAAT...ACG'), id='rc_gi|2765587|emb|Z78462.1|PSZ78462', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCA...ACG'), id='rc_gi|2765586|emb|Z78461.1|PWZ78461', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCA...ACG'), id='rc_gi|2765585|emb|Z78460.1|PCZ78460', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAAC...ACG'), id='rc_gi|2765584|emb|Z78459.1|PDZ78459', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAAC...ACG'), id='rc_gi|2765583|emb|Z78458.1|PHZ78458', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGA...ACG'), id='rc_gi|2765582|emb|Z78457.1|PCZ78457', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAA...ACG'), id='rc_gi|2765581|emb|Z78456.1|PTZ78456', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCA...ACG'), id='rc_gi|2765580|emb|Z78455.1|PJZ78455', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGG...ACG'), id='rc_gi|2765579|emb|Z78454.1|PFZ78454', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765578|emb|Z78453.1|PSZ78453', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATC...ACG'), id='rc_gi|2765577|emb|Z78452.1|PBZ78452', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCA...ACG'), id='rc_gi|2765576|emb|Z78451.1|PHZ78451', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCA...TCC'), id='rc_gi|2765575|emb|Z78450.1|PPZ78450', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765574|emb|Z78449.1|PMZ78449', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTG...ACG'), id='rc_gi|2765573|emb|Z78448.1|PAZ78448', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_gi|2765572|emb|Z78447.1|PVZ78447', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGG...ACG'), id='rc_gi|2765571|emb|Z78446.1|PAZ78446', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCA...ACG'), id='rc_gi|2765570|emb|Z78445.1|PUZ78445', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCT...ACG'), id='rc_gi|2765569|emb|Z78444.1|PAZ78444', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGG...ACG'), id='rc_gi|2765568|emb|Z78443.1|PLZ78443', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATC...TAC'), id='rc_gi|2765567|emb|Z78442.1|PBZ78442', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCA...TCC'), id='rc_gi|2765566|emb|Z78441.1|PSZ78441', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_gi|2765565|emb|Z78440.1|PPZ78440', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATA...ATG'), id='rc_gi|2765564|emb|Z78439.1|PBZ78439', name='<unknown name>', description='reverse complement', dbxrefs=[])]



```python
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]
```


```python
# To see which records have less than 700 length in nucleotides
len(records)
```




    18




```python
# a complete example:

records = (
rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
if len(rec) < 700
)

SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    18




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
