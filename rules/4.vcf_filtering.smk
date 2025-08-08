rule removeLowQual:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        os.path.join(config["outdir"], "vcf", config["project"]+".raw"+".vcf.gz")
    output:
        os.path.join(config["outdir"], "vcf_filtered", config["project"]+".removeLowQual"+".vcf.gz")
    threads:
        1
    shell:
        """
        bcftools filter --exclude FILTER="LowQual" {input} -o {output}
        """


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
                    start = int(parts[1])
                    end = int(parts[2])

                    outfile.write(f"{seq_name}\t{start}\t{end}\n")
        acclist_to_bed(input.acclist, output.bed)

rule removeLowComplexityVariant:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        vcf = os.path.join(config["outdir"], "vcf_filtered", config["project"]+".removeLowQual"+".vcf.gz"),
        bed = os.path.join(config["outdir"],"ref","ref.lcm.bed")
    output:
        os.path.join(config["outdir"], "vcf_filtered", config["project"]+".removeLowQual"+".lcm"+".vcf.gz")
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

rule getRefCalls:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        os.path.join(config["outdir"], "vcf_filtered", config["project"]+".removeLowQual"+".lcm"+".vcf.gz")
    output:
        os.path.join(config["outdir"], "vcf_filtered", config["project"]+".removeLowQual"+".lcm"+".refCalls"+".vcf.gz")
    threads:
        1
    shell:
        """
            bcftools view \
                --max-ac=1 \
                {input} \
                -Oz -o {output}
        """


rule getHQSNPs:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        os.path.join(config["outdir"], "vcf_filtered", config["project"]+".removeLowQual"+".lcm"+".vcf.gz")
    output:
        os.path.join(config["outdir"], "vcf_filtered", config["project"]+".removeLowQual"+".lcm"+".HQSNPs"+".vcf.gz")
    threads:
        1
    shell:
        """
            bcftools filter \
                -S . -e 'FMT/GQ<20' \
                {input} | \
            bcftools filter \
                -s 'FAIL' \
                -e 'QUAL<30 || N_ALLELES>2 || AVG(FMT/DP)>50 || TYPE!="snp" || F_MISSING>0.5' | \
            bcftools view \
                -f 'PASS' \
                -Oz -o {output} 
        """

rule mergeRefCallsAndHQSNPs:
    conda:
        os.path.join(workflow.basedir,"envs/gatk4.yaml")
    input:
        non_variant=os.path.join(config["outdir"], "vcf_filtered", config["project"]+".removeLowQual"+".lcm"+".refCalls"+".vcf.gz"),
        variant=os.path.join(config["outdir"], "vcf_filtered", config["project"]+".removeLowQual"+".lcm"+".HQSNPs"+".vcf.gz")
    output:
        os.path.join(config["outdir"], "vcf_final", config["project"]+".removeLowQual"+".lcm"+".allSite"+".vcf.gz")
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
