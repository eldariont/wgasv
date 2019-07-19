rule run_minimap2:
    input:
        amb="assemblies_{chr}/{sample}_paf_{chr}.fa",
        ref="assemblies_{chr}/hg38_{chr}.fa"
    output:
        "minimap2_output/{sample}.{chr}.{setting}.hg38.paf"
    threads: 10
    shell:
        "minimap2 -t {threads} -c --cs -x {wildcards.setting} {input.ref} {input.amb} > {output}"


rule run_paftools_call:
    input:
        paf="minimap2_output/{sample}.{chr}.{setting}.hg38.paf",
        ref="assemblies_{chr}/hg38_{chr}.fa"
    output:
        "paftools_output/{sample}.{chr}.{setting}.hg38.vcf"
    params:
        min_length_coverage = 10000,
        min_length_call = 50000,
        min_mapq = 5
    shell:
        "sort -k6,6 -k8,8n {input.paf} | paftools.js call -f {input.ref} -s {wildcards.sample} -l {params.min_length_coverage} -L {params.min_length_call} -q {params.min_mapq} - | awk 'OFS=\"\\t\" {{ if($1 ~ /^#/) {{ print $0 }} else {{ NUM+=1; printf \"%s\\t%s\\t{wildcards.sample}.%s\", $1, $2, NUM; for(i=4;i<=NF;i++) {{ printf \"\\t%s\", $i }}; printf \"\\n\" }} }}' > {output}"


rule select_svs:
    input:
        "paftools_output/{sample}.{chr}.{setting}.hg38.vcf"
    output:
        "paftools_output/{sample}.{chr}.{setting}.hg38.svs.vcf"
    shell:
        "bcftools view -i 'STRLEN(REF)>10 | STRLEN(ALT)>10' {input} > {output}"


rule merge_svs:
    input:
        ref="assemblies_{chr}/hg38_{chr}.fa",
        svs=expand("paftools_output/{sample}.{{chr}}.{{setting}}.hg38.svs.vcf", sample=SAMPS)
    output:
        filelist=temp("paftools_output/svs.{chr}.{setting}.fofn"),
        merged="paftools_output/merged.{chr}.{setting}.hg38.svs.clustered.vcf"
    params:
        prefix="paftools_output/merged.{chr}.{setting}.hg38.svs"
    run:
        filelist_open = open(output.filelist, "w")
        for f in input.svs:
            print(f, file=filelist_open)
        filelist_open.close()
        shell("svanalyzer merge --ref {input.ref} --fof {output.filelist} --prefix {params.prefix}")


rule select_small_variants:
    input:
        vcf="paftools_output/{sample}.{chr}.{setting}.hg38.vcf",
        ref="assemblies_{chr}/hg38_{chr}.fa"
    output:
        "paftools_output/{sample}.{chr}.{setting}.hg38.small.vcf.gz"
    shell:
        "bcftools view -i 'STRLEN(REF)<=10 & STRLEN(ALT)<=10' --output-type z {input.vcf} | bcftools norm --fasta-ref {input.ref} --rm-dup none --output-type z > {output}"


rule tabix:
    input:
        "{name}.vcf.gz"
    output:
        "{name}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"


rule merge_small_variants:
    input:
        vcfs=expand("paftools_output/{sample}.{{chr}}.{{setting}}.hg38.small.vcf.gz", sample=SAMPS),
        tbis=expand("paftools_output/{sample}.{{chr}}.{{setting}}.hg38.small.vcf.gz.tbi", sample=SAMPS)
    output:
        "paftools_output/merged.{chr}.{setting}.hg38.small.clustered.vcf"
    shell:
        "bcftools merge --merge none {input.vcfs} > {output}"


rule merge_svs_and_small_variants:
    input:
        small="paftools_output/merged.{chr}.{setting}.hg38.small.clustered.vcf",
        svs="paftools_output/merged.{chr}.{setting}.hg38.svs.clustered.vcf"
    output:
        temp_small=temp("paftools_output/merged.{chr}.{setting}.hg38.small.clustered.nogt.vcf.gz"),
        temp_svs=temp("paftools_output/merged.{chr}.{setting}.hg38.svs.clustered.nogt.vcf.gz"),
        temp_small_tbi=temp("paftools_output/merged.{chr}.{setting}.hg38.small.clustered.nogt.vcf.gz.tbi"),
        temp_svs_tbi=temp("paftools_output/merged.{chr}.{setting}.hg38.svs.clustered.nogt.vcf.gz.tbi"),
        final="paftools_output/merged.{chr}.{setting}.final.vcf.gz"
    run:
        #Remove GT fields
        shell("awk 'OFS=\"\\t\" {{ if($1 ~ /^##/) {{ print $0 }} else {{ print $1, $2, $3, $4, $5, $6, $7, $8 }} }}' {input.small} | bcftools sort --output-type z > {output.temp_small}")
        shell("awk 'OFS=\"\\t\" {{ if($1 ~ /^##/) {{ print $0 }} else {{ print $1, $2, $3, $4, $5, $6, $7, $8 }} }}' {input.svs} | bcftools sort --output-type z > {output.temp_svs}")
        shell("tabix -p vcf {output.temp_small} && tabix -p vcf {output.temp_svs}")
        shell("bcftools concat --allow-overlaps --output-type z {output.temp_small} {output.temp_svs} > {output.final}")