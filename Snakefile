import pandas as pd

configfile: "config.yaml"




samples_df = pd.read_table(config["sample_tsv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])


rule all:
    input:
        expand("trimmed_merged/{sample}_R1.fastq.gz", sample=sample_list),
        expand("trimmed_merged/{sample}_R2.fastq.gz", sample=sample_list),
        #expand("alignment/{sample}.srt.bam", sample=sample_list),
        #expand("alignment/{sample}.srt.bam.bai", sample=sample_list),
#        expand("alignment/stats/{sample}.srt.bam.flagstat", sample=sample_list),
        expand("markdup/{sample}.mkdup.bam", sample=sample_list),
#        expand("metrics/picard/{sample}.picard.rna.metrics.txt", sample=sample_list),
#        "featurecounts/featurecounts.readcounts.tsv"
        expand("freebayes/{sample}.raw.vcf", sample=sample_list),
        expand("assembly/{sample}/scaffolds.fasta", sample=sample_list),
        expand("assembly/{sample}/scaffolds_to_ref.bam", sample=sample_list),
        expand("lumpy/{sample}.raw.vcf", sample=sample_list),
    output:
        "multiqc_report.html"
    shell: """
        multiqc fastqc alignment markdup metrics featurecounts
"""


rule trim_merge:
    output: "trimmed_merged/{sample}_R1.fastq.gz",
            "trimmed_merged/{sample}_R2.fastq.gz"
    params:
        sample = lambda wildcards:  wildcards.sample,
        trimmomatic = config["trimmomatic_path"],
        trimmomatic_adapters = config["trimmomatic_adapters"],
        fastq_list_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq_list_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
    resources: threads="6", maxtime="2:00:00", memory="6gb",
    shell: """
        java -jar {params.trimmomatic}  PE -threads 4 {params.fastq_list_1} {params.fastq_list_2} trimmed_merged/{params.sample}_R1.fastq.gz trimmed_merged/{params.sample}_fwdtmp.fastq.gz  trimmed_merged/{params.sample}_R2.fastq.gz trimmed_merged/{params.sample}_revtmp.fastq.gz    ILLUMINACLIP:{params.trimmomatic_adapters}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
        rm -f trimmed_merged/{params.sample}_fwdtmp.fastq.gz trimmed_merged/{params.sample}_revtmp.fastq.gz


"""

rule alignment:
    input: "trimmed_merged/{sample}_R1.fastq.gz",
            "trimmed_merged/{sample}_R2.fastq.gz"
    output: "alignment/{sample}.srt.bam",
            "alignment/{sample}.srt.bam.bai"
    params:
#        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
#        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
        layout = config["layout"],
        sample = lambda wildcards:  wildcards.sample,
        aligner_name = config["aligner_name"],
        aligner = config["aligner_path"],
        aligner_index = config["aligner_index"],
        samtools = config["samtools_path"],
#        rg = "@RG\tID:{sample}\tSM:{sample}"
    resources: threads="6", maxtime="8:00:00", memory="8gb",
    shell: """

           if [ "{params.layout}" = "single" ]
            then
           {params.aligner} -x {params.aligner_index} --rg ID:{params.sample} --rg SM:{params.sample} --rg LB:{params.sample}  -U trimmed_merged/{params.sample}_R1.fastq.gz -p 8  --summary-file alignment/{params.sample}.hisat.summary.txt | {params.samtools} view -@ 2 -b | {params.samtools} sort -T /scratch/samtools_{params.sample} -@ 8 - 1> alignment/{params.sample}.srt.bam
           {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam
            else
            {params.aligner} mem  {params.aligner_index}   trimmed_merged/{params.sample}_R1.fastq.gz  trimmed_merged/{params.sample}_R2.fastq.gz  -t 8 | {params.samtools} view -@ 2 -b | {params.samtools} sort -T /scratch/samtools_{params.sample} -@ 4 -m 512M  - 1> alignment/{params.sample}.srt.bam
           {params.samtools} index -@ 4 alignment/{params.sample}.srt.bam

fi
""" 


rule picard_markdup:
    input: "alignment/{sample}.srt.bam"
    output: "markdup/{sample}.mkdup.bam"
    params:
        sample = lambda wildcards:  wildcards.sample,
        picard = config['picard_path'],
        java = config['java_path']
    resources: threads="8", maxtime="24:00:00", memory="12gb",
    shell: """
            {params.java} -Xmx8G -Xms8G   -jar {params.picard} MarkDuplicates I=alignment/{params.sample}.srt.bam O=markdup/{params.sample}.mkdup.bam M=markdup/{params.sample}.mkdup.log.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=true  MAX_RECORDS_IN_RAM=4000000 ASSUME_SORTED=true MAX_FILE_HANDLES=768
            rm -f alignment/{params.sample}.srt.bam
"""

rule freebayes:
    input: "markdup/{sample}.mkdup.bam"
    output: "freebayes/{sample}.raw.vcf"
    params:
        sample = lambda wildcards:  wildcards.sample,
        freebayes = config['freebayes_path'],
        ref_fa = config['reference_fasta'],
    resources: threads="4", maxtime="24:00:00", memory="4gb",
    shell: """
            {params.freebayes} --min-coverage 5 --limit-coverage 100 --min-alternate-fraction .2 --min-mapping-quality 15 --min-alternate-count 2 -f {params.ref_fa}  markdup/{params.sample}.mkdup.bam  > freebayes/{params.sample}.raw.vcf
#            rm -f alignment/{params.sample}.srt.bam
"""

rule lumpy:
    input: "markdup/{sample}.mkdup.bam"
    output: "lumpy/{sample}.raw.vcf"
    params:
        sample = lambda wildcards:  wildcards.sample,
        lumpyexp = config['lumpyexp_path'],
        lumpy_scripts = config['lumpy_scripts'],
        ref_fa = config['reference_fasta'],
        samtools = config['samtools_path'],
    resources: threads="4", maxtime="4:00:00", memory="4gb",
    shell: """
        {params.samtools} addreplacerg -r 'ID:1' -r 'LB:1' -r 'SM:1' -o lumpy/{params.sample}.rg.bam markdup/{params.sample}.mkdup.bam
        {params.samtools} index lumpy/{params.sample}.rg.bam

        {params.samtools} view -b -F 1294 lumpy/{params.sample}.rg.bam > lumpy/{params.sample}.discordants.unsorted.bam
        {params.samtools} view -h lumpy/{params.sample}.rg.bam | {params.lumpy_scripts}/extractSplitReads_BwaMem -i stdin | {params.samtools} view -Sb - > lumpy/{params.sample}.splitters.unsorted.bam

        {params.samtools} sort -T /scratch/samtools_lumpyd_{params.sample} -m 512M lumpy/{params.sample}.discordants.unsorted.bam > lumpy/{params.sample}.discordants.bam
        {params.samtools} sort -T /scratch/samtools_lumpys_{params.sample} -m 512M lumpy/{params.sample}.splitters.unsorted.bam > lumpy/{params.sample}.splitters.bam

        {params.lumpyexp} -B lumpy/{params.sample}.rg.bam -S lumpy/{params.sample}.splitters.bam -D lumpy/{params.sample}.discordants.bam -o lumpy/{params.sample}.raw.vcf
        
"""


rule spades:
    input: "trimmed_merged/{sample}_R1.fastq.gz",
            "trimmed_merged/{sample}_R2.fastq.gz"
    output: "assembly/{sample}/scaffolds.fasta",
    params:
#        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
#        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
        layout = config["layout"],
        sample = lambda wildcards:  wildcards.sample,
        spades = config["spades_path"],
#        aligner_index = config["aligner_index"],
        samtools = config["samtools_path"],
    resources: threads="8", maxtime="16:00:00", memory="12gb",
    shell: """
        {params.spades} -t 8 --only-assembler -1 trimmed_merged/{params.sample}_R1.fastq.gz  -2 trimmed_merged/{params.sample}_R2.fastq.gz -o assembly/{params.sample}
"""

rule minimap:
    input: "assembly/{sample}/scaffolds.fasta"
    output: "assembly/{sample}/scaffolds_to_ref.bam"
    params:
        sample = lambda wildcards:  wildcards.sample,
        ref =  config['reference_fasta'],
        minimap = config["minimap2_path"],
        samtools = config["samtools_path"],
    resources: threads="2", maxtime="1:00:00", memory="2gb",
    shell: """
        {params.minimap} -x asm5 -a {params.ref} assembly/{params.sample}/scaffolds.fasta | {params.samtools} sort > "assembly/{params.sample}/scaffolds_to_ref.bam"
        {params.samtools} index assembly/{params.sample}/scaffolds_to_ref.bam
"""

