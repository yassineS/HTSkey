# -*- coding: utf-8 -*-

import os
import pysam

from cosmos import Execution, Cosmos, rel, Recipe, Input
import tools



def _getHeaderInfo(input_bam):
    if   input_bam[-3:] == 'bam':
        header = pysam.Samfile(input_bam,'rb', check_sq = False).header
    elif input_bam[-3:] == 'sam':
        header = pysam.Samfile(input_bam,'r' , check_sq = False).header
    else:
        raise TypeError, 'input file is not a bam or sam'

    return {'rg': [ [tags['ID'], tags['SM'], tags.get('LB','noLBinfo'), tags.get('PL','noPLinfo') ] for tags in header['RG']],
            'sq': [ [tags['SN']                                                                   ] for tags in header['SQ']]
           }

def _getSeqName(header):
    """
    Return sequence names (@SQ SN in header)
    """
    seqNameList = []
    unMapped=''
    for sn in header['sq']:
        if (sn[0].startswith('GL')) or (sn[0].startswith('chrUn')):
            unMapped += " %s" % sn[0]
        else:
            seqNameList.append(sn[0])  # first column is seqName

    if unMapped != '':
        seqNameList.append(unMapped)

    return seqNameList

<<<<<<< HEAD
def genomkey(execution, bams, test_bam=False, chromosome_only_split=False):
=======
def genomkey(execution, test_bam=False, chromosome_only_split=False):
>>>>>>> 8f6e213ebc59785ee1c90b864387ec78c46dc5ea
    recipe = Recipe()
#
# def genomekey(bams, test_bam=False, chromosome_only_split=False):

    # split_ tuples
    chrom  = ('chrom', range(1,23) + ['X', 'Y', 'MT'])
    chrom  = ('chrom', range(1,23))

    glm = ('glm', ['SNP', 'INDEL'])

    dbnames = ('dbname', ['dbSNP135','CytoBand','Target_Scan','mirBase','Self_Chain','Repeat_Masker','TFBS','Segmental_Duplications','SIFT','COSMIC',
                          'PolyPhen2','Mutation_Taster','GERP','PhyloP','LRT','Mce46way','Complete_Genomics_69','The_1000g_Febuary_all','The_1000g_April_all',
                          'NHLBI_Exome_Project_euro','NHLBI_Exome_Project_aa','NHLBI_Exome_Project_all','ENCODE_DNaseI_Hypersensitivity','ENCODE_Transcription_Factor',
                          'UCSC_Gene','Refseq_Gene','Ensembl_Gene','CCDS_Gene','HGMD_INDEL','HGMD_SNP','GWAS_Catalog'])
    bam_seq = None

    for b in bams:
        header = _getHeaderInfo(b)
        sn     = _getSeqName(header)

        rgid = [ h[0] for h in header['rg']]

        # restrict output for testing
        #if test_bam:
            sn    = ['chr1']
            chrom = ('chrom',[1])
            glm   = ('glm',['SNP'])
            skip_VQSR = ('skip_VQSR', [True])

        #else:
        #    skip_VQSR = ('skip_VQSR', [False])

        # if seqName is empty, then let's assume that the input is unaligned bam
        # use everything before extension as part of tag

        sample_name = os.path.splitext(os.path.basename(b))[0]

        #if chromosome_only_split:
            # Stop splitting by rgId
        #    bam_bwa_split = [ ('prevSn', sn), ('chromosome_only_split', [True]) ]
        #    indelrealign_reduce =  ['bam']
        #else:

        bam_bwa_split = [ ('rgId', rgid), ('prevSn', sn), ('chromosome_only_split', [False]) ]

        indelrealign_reduce =  ['bam','rgId']

        #s = sequence_( add_([INPUT(b, tags={'bam':sample_name})], stage_name="Load BAMs"),
            #              split_(bam_bwa_split, pipes.Bam_To_BWA))

        #if bam_seq is None:   bam_seq = s
        #else:                 bam_seq = sequence_(bam_seq, s, combine=True)


    """# Previous pipeline
    #pr_pipeline = sequence_(
        bam_seq,
        reduce_split_(indelrealign_reduce, [chrom], pipes.IndelRealigner),
        map_(                                  pipes.MarkDuplicates),
        reduce_(['bam','chrom'],               pipes.BaseQualityScoreRecalibration),
        map_(                                  pipes.ReduceReads),
        reduce_split_(['chrom'], [glm],        pipes.UnifiedGenotyper),
        reduce_(['glm'],                       pipes.VariantQualityScoreRecalibration, tag={'vcf':'main'}),
    ###   reduce_(['vcf'],                       pipes.CombineVariants, Merge VCF),
        map_(                                  pipes.Vcf2Anno_in),
        split_([dbnames],                      pipes.Annotate, tag={'build':'hg19'}),
        reduce_(['vcf'],                       pipes.MergeAnnotations)
    )

    # HaplotypeCaller Pipeline: official for GATK 3.0
    hc_pipeline = sequence_(
        bam_seq,
        reduce_split_(indelrealign_reduce, [chrom], pipes.IndelRealigner),
        map_(                                  pipes.MarkDuplicates),
        reduce_(['bam','chrom'],               pipes.BaseQualityScoreRecalibration),
        map_(                                  pipes.HaplotypeCaller),
        reduce_(['chrom'],                     pipes.GenotypeGVCFs),
        split_([glm, skip_VQSR],               pipes.VariantQualityScoreRecalibration, tag={'vcf':'main'})
    )

    # 2.0 hc pipeline
    test_pipeline =

    return test_pipeline
    """

### Pipeline

    add                 = recipe.add_source([Input(b, tags={'bam':sample_name})])
    split               = recipe.add_stage(bam_bwa_split, parent=add)
    align               = recipe.add_stage(tools.BamToBWA, parent=split, rel=rel.One2one)
    indel_religner      = recipe.add_stage(tools.IndelRealigner, parent=align, rel=rel.Many2many(['bam', 'rgId'], ['chrom']))
    mark_duplicates     = add_stage(tools.MarkDuplicates, parent=indel_realigner, rel=rel.Many2many(['bam', 'rgId'], ['chrom']))
<<<<<<< HEAD
=======
    align               = recipe.add_stage(tools.Bam_To_BWA, parent=split, rel=rel.One2one)
    indel_religner      = recipe.add_stage(tools.Bam_To_BWA, parent=align, rel=rel.Many2many(['bam', 'rgId'], ['chrom']))
    mark_duplicates     = add_stage(tools.IndelRealigner, parent=indel_realigner, rel=rel.Many2many(['bam', 'rgId'], ['chrom']))
>>>>>>> 8f6e213ebc59785ee1c90b864387ec78c46dc5ea
    bqsr                = add_stage(tools.BQSR, parent=mark_duplicates, rel=rel.Many2one(['bam', 'chrom']))
    haplotype_caller    = add_stage(tools.HaplotypeCaller, parent=bqsr, rel=rel.One2one)
    genotype_gvcfs      = add_stage(tools.GenotypeGVCFs, parent=haplotype_caller, rel=Many2one(['chrom']))
    vqsr                = add_stage(tools.VQSR, parent=genotype_gvcfs, rel=One2many([glm, skip_VQSR]), tag={'vcf':'main'})



    execution.run(recipe)


if __name__ == '__main__':
    cosmos_app = Cosmos('sqlite:///sqlite.db', default_drm='local')
    cosmos_app.initdb()

    ex = Execution.start(cosmos_app, 'Simple', 'out/simple', max_attempts=3, restart=True, skip_confirm=True)
    genomkey(ex)


