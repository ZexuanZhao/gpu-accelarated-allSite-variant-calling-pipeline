rule lowComplexityAcclist:
    conda:
        os.path.join(workflow.basedir,"envs/blast.yaml")
    input:
        os.path.join(config["outdir"], "ref", "ref.fasta")
    output:
        os.path.join(config["outdir"], "ref", "ref.lcm.acclist")
    threads:
        1
    shell:
        """
        dustmasker \
            -in {input} \
            -outfmt acclist \
            > {output}
        """
rule lowComplexityBed:
    input:
        acclist = os.path.join(config["outdir"],"ref","ref.lcm.acclist")
    output:
        bed = os.path.join(config["outdir"],"ref","ref.lcm.bed")
    threads:
        1
    run:
        def acclist_to_bed(input_file_path, output_file_path):
            with open(input_file_path,"r") as infile, open(output_file_path,"w") as outfile:
                for line in infile:
                    line = line.strip()  # Remove leading/trailing whitespace
                    if not line:  # Skip empty lines.
                        continue
                    parts = line.split("\t")
                    seq_name = parts[0].replace(">","")
                    seq_name = seq_name.split(" ")[0]
                    start = int(parts[1]) - 1
                    end = int(parts[2])

                    outfile.write(f"{seq_name}\t{start}\t{end}\n")
        acclist_to_bed(input.acclist, output.bed)

rule removeLowComplexityVariant:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        vcf = os.path.join(config["outdir"],"vcf",config["project"]+".allSite"+".vcf.gz"),
        bed = os.path.join(config["outdir"],"ref","ref.lcm.bed")
    output:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" +".vcf.gz")
    threads:
        1
    shell:
        """
            bedtools subtract \
                -a {input.vcf} \
                -b {input.bed} \
                -A -header | \
                gzip \
            > {output}
        """

rule getNonVariant:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" +".vcf.gz")
    output:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" + ".nonVariant" + ".vcf.gz")
    threads:
        1
    shell:
        """
            bcftools view \
                --max-ac=1 \
                {input} | \
                gzip \
                >{output}
        """


rule QualityFiltering:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" +".vcf.gz")
    output:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" + ".variantHQ" + ".vcf.gz")
    threads:
        1
    shell:
        """
            vcftools \
            --gzvcf {input} \
            --stdout --recode --recode-INFO-all \
            --minQ 30 \
            --max-alleles 2 \
            --max-meanDP 50 \
            --minGQ 20 \
            --max-missing 0.5 | \
            gzip \
            > {output}
        """

rule mergeHQandNonVariant:
    conda:
        os.path.join(workflow.basedir,"envs/gatk4.yaml")
    input:
        non_variant=os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" + ".nonVariant" + ".vcf.gz"),
        variant=os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" + ".variantHQ" + ".vcf.gz")
    output:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" + ".HQ" + ".vcf.gz")
    threads:
        1
    log:
        os.path.join(config["outdir"],"logs","mergeVcfs","mergeHQvariantAndNonVariant.log")
    resources:
        mem_gb = int(round(get_total_memory_gb() * 0.8, 0))
    shell:
        """
            gatk \
            --java-options "-Xmx{resources.mem_gb}g" \
            MergeVcfs \
            I={input.non_variant} I={input.variant} \
            O={output} \
            > {log} \
            2>{log}
        """

rule SNPonly:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".lcm" + ".HQ" + ".vcf.gz")
    output:
        os.path.join(config["outdir"],"vcf",config["project"] + ".allSite" + ".SNPonly" + ".lcm" + ".HQ" + ".vcf.gz")
    threads:
        1
    shell:
        """
            vcftools \
            --gzvcf {input} \
            --stdout --recode --recode-INFO-all \
            --remove-indels | \
            gzip \
            > {output}
        """
