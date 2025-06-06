## Quality check of trimmed reads
rule fastqc_after_trimming:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        unpack(get_fastq2)
    output:
        r1 = os.path.join(config["outdir"],"qc", "fastqc", "{sample}.1P_fastqc.zip"),
        r2 = os.path.join(config["outdir"],"qc", "fastqc", "{sample}.2P_fastqc.zip"),
    params:
        outdir = os.path.join(config["outdir"], "qc", "fastqc")
    threads:
        2
    shell:
        """
        fastqc --quiet --outdir {params.outdir} --noextract -f fastq {input} -t {threads}
        """

## Mapping rates etc.
rule bamstats:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"], "bam", "{sample}.bam")
    output:
        os.path.join(config["outdir"], "qc", "bamtools", "{sample}_bamtools.stats")
    threads:
        1
    shell:
        """
        bamtools stats -in {input} | grep -v "*" > {output}
        """

## Make bed windows for window-based metrics: coverage
rule make_ref_window:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"ref","ref.fasta.fai")
    output:
        os.path.join(config["outdir"],"ref","ref.window.bed")
    params:
        w_size = config["w_size"]
    threads:
        1
    shell:
        """
        bedtools makewindows -g {input} -w {params.w_size} > {output}
        """

## Coverage
rule genomeCov:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        bam = os.path.join(config["outdir"], "bam", "{sample}.bam"),
        bed = os.path.join(config["outdir"],"ref","ref.window.bed")
    output:
        os.path.join(config["outdir"], "qc", "coverage", "{sample}_coverage.txt")
    threads:
        20
    shell:
        """
        bedtools coverage \
            -a {input.bed} -b {input.bam} \
            > {output}
        """

## Qualimap: insert size, GC, coverage...
rule Qualimap:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        bam=os.path.join(config["outdir"],"bam","{sample}.bam")
    output:
        os.path.join(config["outdir"],"qc","qualimap","{sample}", "qualimapReport.html")
    params:
        outdir =  os.path.join(config["outdir"],"qc","qualimap","{sample}")
    threads:
        10
    log:
        os.path.join(config["outdir"],"logs","Qualimap","{sample}.log")
    shell:
        """
        qualimap bamqc \
            -bam {input} \
            -outdir {params.outdir}\
            -outformat HTML \
            -nt {threads} \
            --java-mem-size=50G \
            > {log} \
            2>{log}
        """

## vcf stats using bcftools
rule vcf_stats:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" + ".HQ" + ".vcf.gz")
    output:
        os.path.join(config["outdir"],"qc", "bcftools_stats", config["project"] + ".allSite" + ".lcm" + ".HQ" + ".vcf.stats")
    threads:
        1
    shell:
        """
        bcftools stats {input} > {output}
        """

## Visualize vcf stats
rule plot_vcfstats:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"qc", "bcftools_stats", config["project"] + ".allSite" + ".lcm" + ".HQ" + ".vcf.stats")
    output:
        os.path.join(config["outdir"],"qc","bcftools_stats", "plot-vcfstats.log")
    params:
        outdir = os.path.join(config["outdir"],"qc","bcftools_stats")
    threads:
        1
    shell:
        """
        plot-vcfstats \
            -p {params.outdir} \
            --no-PDF \
            {input}
        """

## Quality check of vcf files using vcftools
rule qc_vcf:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".SNPonly" + ".lcm" + ".HQ" + ".vcf.gz")
    output:
        os.path.join(config["outdir"], "qc", "vcftools", config["project"]+".het")
    params:
        outdir = os.path.join(config["outdir"], "qc", "vcftools"),
        proj = config["project"]
    threads:
        1
    log:
        os.path.join(config["outdir"],"logs","vcftools","vcftools.log")
    shell:
        """
        vcftools --gzvcf {input} --freq2 --max-alleles 2 --out {params.outdir}/{params.proj} >{log} 2>{log}
        vcftools --gzvcf {input} --depth --out {params.outdir}/{params.proj} >{log} 2>{log}
        vcftools --gzvcf {input} --site-mean-depth --out {params.outdir}/{params.proj} >{log} 2>{log}
        vcftools --gzvcf {input} --site-quality --out {params.outdir}/{params.proj} >{log} 2>{log}
        vcftools --gzvcf {input} --missing-indv --out {params.outdir}/{params.proj} >{log} 2>{log}
        vcftools --gzvcf {input} --missing-site --out {params.outdir}/{params.proj} >{log} 2>{log}
        vcftools --gzvcf {input} --het --out {params.outdir}/{params.proj} >{log} 2>{log}
        """

## Summarize all qc files using multiqc
rule multiqc:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        expand(os.path.join(config["outdir"],"qc","qualimap","{sample}", "qualimapReport.html"), sample= sample_sheet.index),
        expand(os.path.join(config["outdir"],"qc", "fastqc", "{sample}.{R}_fastqc.zip"), sample= sample_sheet.index, R=["1P", "2P"]),
        expand(os.path.join(config["outdir"], "qc", "bamtools","{sample}_bamtools.stats"), sample= sample_sheet.index),
        expand(os.path.join(config["outdir"],"qc","fastp","{sample}.fastp.json"), sample= sample_sheet.index),
        os.path.join(config["outdir"],"qc", "bcftools_stats", config["project"] + ".allSite" + ".lcm" + ".HQ" + ".vcf.stats")
    output:
        os.path.join(config["outdir"],"qc","multiqc", config["project"]+"_multiqc_report.html")
    params:
        input_dir = os.path.join(config["outdir"], "qc"),
        output_dir = os.path.join(config["outdir"], "qc", "multiqc"),
        original_output = os.path.join(config["outdir"],"qc", "multiqc", "multiqc_report.html")
    threads:
        1
    log:
        os.path.join(config["outdir"],"logs","multiqc","multiqc.log")
    shell:
        """
        rm -rf {params.output_dir}/* ; \
        multiqc \
        -o {params.output_dir} \
        {params.input_dir} \
        >{log} 2>{log}; \
        mv {params.original_output} {output}
        """
