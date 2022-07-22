import numpy
import random
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

## function to create random DNA sequences
_nuc = numpy.array(["A", "T", "C", "G"])
def randSeq(length):
    seqArray = _nuc[[int(random.random()*4) for i in range(int(length))]]
    return(MutableSeq("".join(seqArray), generic_dna))

## make a reference sequence
ref = randSeq(2*10e5)
## write it to ref.fa
recs = [SeqRecord(ref, id='chr1', description='')]
SeqIO.write(recs, "ref.fa", "fasta")

## prepare some SVs
out_vcf = open('svs.vcf', 'w')
out_vcf.write('##fileformat=VCFv4.1\n')
out_vcf.write('##contig=<ID=chr1,length={}>\n'.format(len(ref)))
out_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
out_vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">\n')
out_vcf.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">\n')
out_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP\n')
## one 300bp het insertion at ~20% of the chromosome
ins_pos = int(len(ref)*.2)
ins_seq = randSeq(300)
ref_s = ref[ins_pos]
alt_s = ref_s + str(ins_seq)
out_vcf.write('chr1\t{pos}\t.\t{ref}\t{alt}\t100\t.\tSVTYPE=INS;SVLEN={svl}\tGT\t0/1\n'.format(pos=ins_pos+1, ref=ref_s, alt=alt_s, svl=len(ins_seq)))
## one 500bp hom deletion at ~80% of the chromosome
del_start = int(len(ref)*.8)
del_end = del_start + 501
ref_s = ref[del_start:del_end]
alt_s = ref_s[0]
out_vcf.write('chr1\t{pos}\t.\t{ref}\t{alt}\t100\t.\tSVTYPE=DEL\tGT\t1/1\n'.format(pos=del_start+1, ref=ref_s, alt=alt_s))
out_vcf.close()

## make two assembled haplotypes and let's say the last 90-95% of the reference is not covered
hap1 = ref[:int(len(ref)*.90)]
hap2 = ref[:int(len(ref)*.90)]
## both have the deletion
hap1 = hap1[:del_start] + hap1[del_end:]
hap2 = hap2[:del_start] + hap2[del_end:]
## hap 1 has the insertion fragmented in 3 pieces a few bp apart
ins_seq_1 = ins_seq[:100]
ins_seq_2 = ins_seq[100:200]
ins_seq_3 = ins_seq[200:]
hap1 = hap1[:ins_pos] + ins_seq_1 + hap1[ins_pos:(ins_pos+5)] + ins_seq_2 + hap1[(ins_pos+5):(ins_pos+15)] + \
    ins_seq_3 + hap1[(ins_pos+15):]
## write fastas with the extra 5% end
recs = [SeqRecord(hap1, id='hap1_1', description='')]
recs.append(SeqRecord(ref[int(len(ref)*.95):], id='hap1_2', description=''))
SeqIO.write(recs, "hap1.fa", "fasta")
recs = [SeqRecord(hap2, id='hap2', description='')]
recs.append(SeqRecord(ref[int(len(ref)*.95):], id='hap2_2', description=''))
SeqIO.write(recs, "hap2.fa", "fasta")

## dummy centromere and tandem repeat files
with open('centromere_ref.txt', 'w') as outf :
    outf.write('chr1\t{}\t{}\t.\tacen\n'.format(int(len(ref)*.45), int(len(ref)*.45)+100))
    outf.write('chr1\t{}\t{}\t.\tacen\n'.format(int(len(ref)*.55), int(len(ref)*.55)+100))
with open('trf_ref.bed', 'w') as outf:
    outf.write('chr1\t{}\t{}\n'.format(int(len(ref)*.4), int(len(ref)*.4)+30))
    outf.write('chr1\t{}\t{}\n'.format(int(len(ref)*.7), int(len(ref)*.7)+100))
