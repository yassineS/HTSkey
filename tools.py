from cosmos import Tool, abstract_input_taskfile as inp, abstract_output_taskfile as out
from main import settings as s


cmd_init = r"""
            set -e -o pipefail && tmpDir=$(mktemp -d --tmpdir={s[scratch]}) && export TMPDIR=$tmpDir;
            printf "%s %s at %s\n" "{s[date]}" "$(hostname)" "$tmpDir" | tee /dev/stderr && cd $tmpDir;
            """
cmd_out = r"""

            echo "{s[date]} Moving files to Storage";
            [[ -a $tmpDir/out.bam ]] && mv -f $tmpDir/out.bam $OUT.bam;
            [[ -a $tmpDir/out.bai ]] && mv -f $tmpDir/out.bai $OUT.bai;
            echo "{s[date]} Moving done";
            cd && /bin/rm -rf $tmpDir;
            """
cmd_out_vcf = r"""

            echo "{s[date]} Moving files to Storage";
            [[ -a $tmpDir/out.vcf     ]] && mv -f $tmpDir/out.vcf     $OUT.vcf;
            [[ -a $tmpDir/out.vcf.idx ]] && mv -f $tmpDir/out.vcf.idx $OUT.vcf.idx;
            echo "{s[date]} Moving done";
            cd && /bin/rm -rf $tmpDir;
            """

#######################################
## Non-GATK tools ##

# List2input

def _list2input(l, opt):
    return opt + ("\n" + opt).join(map(lambda x: str(x), l))


# BWA

class Bam_To_BWA(Tool):
    name = "BWA Alignement"
    cpu_req = 4
    mem_req = 12 * 1024
    time_req = 2 * 60
    inputs = [inp(format='bam', n='>0')]
    outputs = [out(format='bam'), out(format='bai')]

    def cmd(self, i, s, p):
        # removed -m MEM option in samtools sort

        if p['chromosome_only_split']:
            # Using first readgroup id

            cmd_rg = r"""
                rg=$({s[samtools]} view -H {i[bam][0]} | grep -w "@RG" | head -n 1 | sed 's/\t/\\t/g') && echo "RG= $rg";
                {s[samtools]} view -f 2 -u              {i[bam][0]} {p[prevSn]} > tmpIn.ubam;
                """

        else:
            cmd_rg = r"""
                rg=$({s[samtools]} view -H {i[bam][0]} | grep -w {p[rgId]} | uniq | sed 's/\t/\\t/g') && echo "RG= $rg";
                {s[samtools]} view -f 2 -u -r {p[rgId]} {i[bam][0]} {p[prevSn]} > tmpIn.ubam;
                """

        cmd_main = r"""
	            cp {s[empty_sam]} empty.sam && echo -e $rg >> empty.sam;
                {s[samtools]} view -u {i[bam][0]} "empty_region" > empty.ubam 2> /dev/null;

                sizeEmpty="$(du -b empty.ubam | cut -f 1)"; printf "Empty bam file size = %-6d\n" "$sizeEmpty";
                sizeTmpIn="$(du -b tmpIn.ubam | cut -f 1)"; printf "Input bam file size = %-6d\n" "$sizeTmpIn";

                [[ "$sizeTmpIn" -gt "$sizeEmpty" ]] &&
                {s[samtools]} sort -n -o -l 0 -@ {self.cpu_req} tmpIn.ubam _shuf |
                {s[bamUtil]} bam2FastQ --in -.ubam --readname --noeof --firstOut /dev/stdout --merge --unpairedout un.fq 2> /dev/null |
                {s[bwa]} mem -p -M -t {self.cpu_req} -R "$rg" {s[reference_fasta]} - |
                {s[samtools]} view -Shu - |
                {s[samtools]} sort    -o -l 0 -@ {self.cpu_req} - _sort > out.bam;

                # If there's no out.bam available, put an empty bam as output;
                [[ ! -a out.bam ]] && ({s[samtools]} view -Sb empty.sam > out.bam 2> /dev/null) || true;

	            {s[samtools]} index out.bam out.bai;

                """

        return (cmd_init + cmd_main + cmd_out), \
               dict(
                   inputs=_list2input(i['bam'], "-I "),
                   s=s,
                   **locals()
               )

# VCF to Annovar

class Vcf_To_Anno_in(Tool):
    name = "Convert VCF to Annovar"
    forward_input = True
    time_req = 12 * 60
    inputs = [inp(format='vcf', n='>0')]
    outputs = [out(format='anno_in')]

    def cmd(self, i, s, p):
        return cmd_init + r"""
             {s[annovarext]} vcf2anno '{i[vcf][0]}' > $OUT.anno_in
            """
# Annovar

class Annotate(Tool):
    name = "Annotate"
    forward_input = True
    time_req = 12 * 60
    mem_req = 12 * 1024
    inputs = [inp(format='anno_in', n='>0')]
    outputs = [out(format='dir')]

    def cmd(self, i, s, p):
        return cmd_init + r"""

             {s[annovarext]} anno {p[build]} {p[dbname]} {i[anno_in][0]} $OUT.dir

              """

# Merge Annotation

class Merge_Annotations(Tool):
    name = "Merge Annotations"
    inputs = ['anno_in', 'dir']
    outputs = ['dir']
    mem_req = 40 * 1024
    time_req = 12 * 60
    forward_input = True
    inputs = [inp(format='anno_in', n='>0'), inp(format='dir')]
    outputs = [out(format='dir')]

    def cmd(self, i, s, p):
        return cmd_init + r"""

              {s[annovarext]} merge {i[anno_in][0]} $OUT.dir {annotated_dir_output}

              """, {'annotated_dir_output': ' '.join(map(str, i['dir']))}


#######################################
## GATK ##

# Indel Realigner

class Indel_Realigner(Tool):
    name = "Indel Realigner"
    cpu_req = 4
    mem_req = 12 * 1024
    time_req = 4 * 60
    inputs = [inp(format='bam', n='>0')]
    outputs = [out(format='bam'), out(format='bai')]

    def cmd(self, i, s, p):
        cmd_main = r"""
            {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx{self.mem_req}M
            -jar {s[gatk]}
            -T RealignerTargetCreator
            -R {s[reference_fasta]}
            -o $tmpDir/{p[chrom]}.intervals
            --known {s[1kindel_vcf]}
            --known {s[mills_vcf]}
            --num_threads {self.cpu_req}
            -L {p[chrom]} {s[gatk_realigntarget]}
            {inputs};

            printf "\n%s RealignerTargetCreator ended.\n" "{s[date]}" | tee -a /dev/stderr;

            {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx{self.mem_req}M
            -jar {s[gatk]}
            -T IndelRealigner
            -R {s[reference_fasta]}
            -o $tmpDir/out.bam
            -targetIntervals $tmpDir/{p[chrom]}.intervals
            -known {s[1kindel_vcf]}
            -known {s[mills_vcf]}
            -model USE_READS
            -compress 0
            -L {p[chrom]} {s[gatk_indelrealign]}
            {inputs};
        """
        return (cmd_init + cmd_main + cmd_out), \
               dict(
                   inputs=_list2input(i['bam'], "-I "),
                   s=s,
                   **locals()
               )


class Mark_Duplicates(Tool):
    name = "MarkDuplicates"
    cpu_req = 2
    mem_req = 4 * 1024
    time_req = 2 * 60
    inputs = [inp(format='bam', n='>0')]
    outputs = [out(format='bam'), out(format='bai'), out(format='metrics')]

    def cmd(self, i, s, p):
        cmd_main = r"""

            {s[java]} -Xmx{self.mem_req}M -jar {s[picard_dir]}/MarkDuplicates.jar
            TMP_DIR=$tmpDir
            OUTPUT=$tmpDir/out.bam
            METRICS_FILE=$tmpDir/out.metrics
            ASSUME_SORTED=True
            CREATE_INDEX=True
            COMPRESSION_LEVEL=0
            MAX_RECORDS_IN_RAM=1000000
            VALIDATION_STRINGENCY=SILENT
            VERBOSITY=INFO
            {inputs};

            mv -f $tmpDir/out.metrics $OUT.metrics;
        """
        return (cmd_init + cmd_main + cmd_out), \
               dict(
                   inputs=_list2input(i['bam'], "INPUT="),
                   s=s,
                   **locals()
               )


class BQSR(Tool):
    name = "Base Quality Score Recalibration"
    cpu_req = 4
    mem_req = 12 * 1024
    time_req = 4 * 60
    inputs = [inp(format='bam', n='>0')]
    outputs = [out(format='bam'), out(format='bai')]  # Need to check the persistence option

    # no -nt, -nct = 4
    def cmd(self, i, s, p):
        cmd_main = r"""

            {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx{self.mem_req}M
            -jar {s[gatk]}
            -T BaseRecalibrator
            -R {s[reference_fasta]}
            -o $tmpDir/{p[chrom]}.grp
            -knownSites {s[dbsnp_vcf]}
            -knownSites {s[1komni_vcf]}
            -knownSites {s[1kindel_vcf]}
            -knownSites {s[mills_vcf]}
            -nct {self.cpu_req}
            -L {p[chrom]}
            {inputs};

            printf "\n%s BaseRecalibrator ended\n" "{s[date]}" | tee -a /dev/stderr;

            {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx{self.mem_req}M
            -jar {s[gatk]}
            -T PrintReads
            -R {s[reference_fasta]}
            -o $tmpDir/out.bam
            -compress 0
            -BQSR $tmpDir/{p[chrom]}.grp
            -nct {self.cpu_req}
            -L {p[chrom]}
            {inputs};

            """
        return (cmd_init + cmd_main + cmd_out), \
               dict(
                   inputs=_list2input(i['bam'], "-I "),
                   s=s,
                   **locals()
               )


# Mean to be used per sample
class Haplotype_Caller(Tool):
    name = "Haplotype Caller"
    cpu_req = 4
    mem_req = 16 * 1024
    time_req = 12 * 60
    inputs = [inp(format='bam', n='>0')]
    outputs = [out(format='vcf'), out(format='vcf.idx')]

    # -nct available
    def cmd(self, i, s, p):
        cmd_main = r"""

            {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx{self.mem_req}M
            -jar {s[gatk]}
            -T HaplotypeCaller
            -R {s[reference_fasta]}
            -D {s[dbsnp_vcf]}
            -o $tmpDir/out.vcf
            -pairHMM VECTOR_LOGLESS_CACHING
            -L {p[chrom]}
            -nct 1
            --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000
            -A DepthPerAlleleBySample
            -stand_call_conf 30
            -stand_emit_conf 10
            {inputs};

            """
        return (cmd_init + cmd_main + cmd_out_vcf), \
               dict(
                   inputs=_list2input(i['vcf'], "-I "),
                   s=s,
                   **locals()
               )

# Joint Genotyping
class Genotype_GVCFs(Tool):
    name = "Genotype GVCFs"
    cpu_req = 4
    mem_req = 12 * 1024
    time_req = 12 * 60
    inputs = [inp(format='vcf', n='>0')]
    outputs = [out(format='vcf'), out(format='vcf.idx')]

    # -nt available
    def cmd(self, i, s, p):
        cmd_main = r"""

            {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx{self.mem_req}M
            -jar {s[gatk]}
            -T GenotypeGVCFs
            -R {s[reference_fasta]}
            -D {s[dbsnp_vcf]}
            -o $tmpDir/out.vcf
            -L {p[chrom]}
            -nt {self.cpu_req}
            -A Coverage
            -A GCContent
            -A HaplotypeScore
            -A MappingQualityRankSumTest
            -A InbreedingCoeff -A FisherStrand -A QualByDepth -A ChromosomeCounts
            {inputs};

            """
        return (cmd_init + cmd_main + cmd_out_vcf), \
               dict(
                   inputs=_list2input(i['vcf'], "-V "),
                   s=s,
                   **locals()
               )

class VQSR(Tool):
    """
    VQSR
    Note that HaplotypeScore is no longer applicable to indels
    see http://gatkforums.broadinstitute.org/discussion/2463/unified-genotyper-no-haplotype-score-annotated-for-indels

    """
    name = "Variant Quality Score Recalibration"
    cpu_req = 4
    mem_req = 12 * 1024
    time_req = 12 * 60
    inputs = [inp(format='vcf', n='>0')]
    outputs = [out(format='vcf'), out(format='vcf.idx'), out(format='R')]


    # -nt available, -nct not available
    def cmd(self, i, s, p):
        """
        Check gatk forum: http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
        --maxGaussians         8 (default), set    1 for small-scale test
        --minNumBadVariants 1000 (default), set 3000 for small-scale test
        """
        cmd_VQSR = r"""

            {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx{self.mem_req}M
            -jar {s[gatk]}
            -T VariantRecalibrator
            -R {s[reference_fasta]}
            -recalFile    $tmpDir/out.recal
            -tranchesFile $tmpDir/out.tranches
            -rscriptFile  $tmpDir/out.R
            -nt {self.cpu_req}
            -an MQRankSum -an ReadPosRankSum -an DP -an FS -an QD
            -mode {p[glm]}
            -L {p[chrom]}
            --maxGaussians 1
            --minNumBadVariants 3000
            {inputs}
            """
        cmd_SNP = r"""
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {s[hapmap_vcf]}
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0  {s[dbsnp_vcf]}
            -resource:omni,known=false,training=true,truth=true,prior=12.0   {s[1komni_vcf]}
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {s[1ksnp_vcf]};
            """
        cmd_INDEL = r"""
            -resource:mills,known=false,training=true,truth=true,prior=12.0 {s[mills_vcf]}
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {s[dbsnp_vcf]};
            """
        cmd_apply_VQSR = r"""

            printf "\n%s\n" "{s[date]}" | tee -a /dev/stderr;

            {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx{self.mem_req}M
            -jar {s[gatk]}
            -T ApplyRecalibration
            -R {s[reference_fasta]}
            -recalFile    $tmpDir/out.recal
            -tranchesFile $tmpDir/out.tranches
            -o            $tmpDir/out.vcf
            --ts_filter_level 99.9
            -mode {p[glm]}
            -nt {self.cpu_req}
            -L {p[chrom]}
            {inputs};

            mv -f $tmpDir/out.R $OUT.R;

            """
        if p['glm'] == 'SNP':
            cmd_rc = cmd_SNP
        else:
            cmd_rc = cmd_INDEL

        if p['skip_VQSR']:
            return " cp {i[vcf][0]} $OUT.vcf; cp {i[vcf][0]}.idx $OUT.vcf.idx; touch $OUT.R"
        else:
            return (cmd_init + cmd_VQSR + cmd_rc + cmd_apply_VQSR + cmd_out_vcf), \
               dict(
                   inputs=_list2input(i['vcf'], "-input "),
                   s=s,
                   **locals()
               )

class Combine_Variants(Tool):
    name = "Combine Variants"
    cpu_req = 4  # max CPU here
    mem_req = 12 * 1024
    time_req = 2 * 60
    inputs = [inp(format='vcf', n='>0')]
    outputs = [out(format='vcf'), out(format='vcf.idx')]

    # -nt available, -nct not available
    # Too many -nt (20?) will cause write error
    def cmd(self, i, s, p):
        """
        :param genotypemergeoptions: select from the following:
            UNIQUIFY       - Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
            PRIORITIZE     - Take genotypes in priority order (see the priority argument).
            UNSORTED       - Take the genotypes in any order.
            REQUIRE_UNIQUE - Require that all samples/genotypes be unique between all inputs.
        """
        cmd_main = r"""

            {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx{self.mem_req}M
            -jar {s[gatk]}
            -T CombineVariants
            -R {s[reference_fasta]}
            -o $tmpDir/out.vcf
            -genotypeMergeOptions UNSORTED
            -nt {self.cpu_req}
            {inputs};

        """
        return (cmd_init + cmd_main + cmd_out_vcf), \
               dict(
                   inputs=_list2input(i['vcf'], "-V "),
                   s=s,
                   **locals()
               )

#######################################
## General tools   ##


class Sleep(Tool):
    def cmd(self, i, o, time=10):
        return 'sleep {time}'


class Echo(Tool):
    outputs = [out('echo', 'txt')]

    def cmd(self, _, outputs, word):
        return '{s[echo_path]} {word} > {outputs[0]}'.format(s=s, **locals())


class Cat(Tool):
    inputs = [inp(format='txt', n='>=1')]
    outputs = [out('cat', 'txt', 'cat_out.txt', )]

    def cmd(self, inputs, outputs):
        input_txts = inputs[0]  # it's easier to just unpack in the method signature, see examples below.
        out_txt = outputs[0]  # it's easier to just unpack in the method signature, see examples below.
        return 'cat {input} > {out_txt}'.format(
            input=' '.join(map(str, (input_txts,))),
            **locals()
        )


class Paste(Tool):
    inputs = [inp(format='txt')]
    outputs = [out('paste', 'txt', 'paste.txt')]

    def cmd(self, (input_txts, ), (out_txt, )):
        return 'paste {input} > {out_txt}'.format(
            input=' '.join(map(str, (input_txts,))),
            **locals()
        )


class WordCount(Tool):
    inputs = [inp(format='txt')]
    outputs = [out('wc', 'txt')]

    def cmd(self, (input_txts, ), (out_txt, )):
        return 'wc {input} > {out_txt}'.format(
            input=' '.join(map(str, (input_txts,))),
            **locals()
        )


class Fail(Tool):
    def cmd(self, inputs, outputs):
        return '__fail__'


class MD5Sum(Tool):
    inputs = [inp(format='*', n=1)]
    outputs = [out(name='md5', format='md5')]

    def cmd(self, (in_file, ), _, out_md5):
        out_md5.basename = in_file.basename + '.md5'
        return 'md5sum {in_file}'.format(**locals())
