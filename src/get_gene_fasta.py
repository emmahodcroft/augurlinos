import numpy as np
from util import load_features, read_in_vcf, safe_translate
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#example: run src/get_gene_fasta.py --path tb_walker21-OG --vcf "tb_walker21-OG/results/tree_nuc_aln.vcf" --vcfRef "tb_walker21-OG/results/ref.fasta" --reference "tb_walker21-OG/data/Mtb_H37Rv_NCBI_Annot.gff" --genes lldD2 esxM esxN

def get_feature_translation(sequences, ref, feature):

    def str_reverse_comp(str_seq):
        #gets the reverse compliment of a string and returns it as
        #a string, saving to-and-fro lines elsewhere
        seq_str = Seq(str_seq)
        return str(seq_str.reverse_complement())

    refNuc = str(feature.extract( SeqRecord(seq=Seq(ref)) ).seq)
    ref_aa_seq = safe_translate(refNuc)

    start = int(feature.location.start)
    end = int(feature.location.end)

    seqNames = sequences.keys()
    tmp_reduced_aln = []
    for pi in range(start,end):

        pattern = [ sequences[k][pi] if pi in sequences[k].keys() else ref[pi] for k,v in sequences.iteritems() ]
        str_pat = "".join(pattern)
        tmp_reduced_aln.append(pattern)

    reduced_alignment = np.array(tmp_reduced_aln).T
    seq_reduce_align = {seqNames[i]:reduced_alignment[i] for i in xrange(len(seqNames)) }

    geneWrite = []
    for i in xrange(len(seqNames)):
        geneWrite.append(">"+seqNames[i])
        if feature.strand == -1:
            geneWrite.append( str_reverse_comp("".join(reduced_alignment[i])) )
        else:
            geneWrite.append( "".join(reduced_alignment[i]) )

    return geneWrite

    #with open("tb_walker21-OG/results/lldD2.fasta", 'w') as the_file:
    #    the_file.write("\n".join(refWrite))


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Creates a fasta alignment of the requested genes from the alignment'
                ' (in correct orientation) which can be translated to amino acids.')
    parser.add_argument('--vcf', required=True, type=str, help="VCF file for the alignment")
    parser.add_argument('--vcfRef', required=True, type=str, help="fasta containing the reference alignment for the VCF")
    parser.add_argument('--reference', required=True, type=str, help="GFF or GenBank file containing the gene information")
    parser.add_argument('--genes', nargs='+', help="genes to be translated")
    parser.add_argument('--path', required=True, type=str, help="path to write out gene alignment(s)")
    args = parser.parse_args()

    #--vcf
    #--vcfRef
    #--reference
    #--genes
    #--path

    path = args.path

    try:
        vcf_dict = read_in_vcf(args.vcf, args.vcfRef, compressed=False )
    except:
        print "Loading input alignment failed!: {}".format(args.vcf)
    ref = vcf_dict['reference']
    sequences = vcf_dict['sequences']

    selected_features = load_features(args.reference)

    for g in args.genes:
        #find the requested gene in selected_features
        geneWrite = get_feature_translation(sequences, ref, selected_features[g])

        with open(path+'/'+g+'.fasta', 'w') as the_file:
            the_file.write("\n".join(geneWrite))


