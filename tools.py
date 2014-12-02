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

###################
# Genomekey tools #
###################

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

    class Bam_To_BWA(Tool):
        name = "BAM to BWA"
        cpu_req = 4  # orchestra: 4
        mem_req = 12 * 1024  # max mem used (RSS): 12.4 Gb?
        time_req = 2 * 60
        inputs = ['bam']
        outputs = ['bam', 'bai']

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
            return (cmd_init + cmd_rg + cmd_main + cmd_out)


### GATK ###

# Indel Realigner

class IndelRealigner(Tool):
    name = "IndelRealigner"
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


###################
# General tools   #
###################

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
