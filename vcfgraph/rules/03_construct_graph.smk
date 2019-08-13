rule construct_graph:
    input:
        fasta="assemblies_chr10/hg38_chr10.fa",
        vcf="paftools_output/merged.chr10.asm20.final.vcf.gz",
        tbi="paftools_output/merged.chr10.asm20.final.vcf.gz.tbi"
    output:
        vg="graph/chr10.asm20.vg"
    shell:
        "vg construct -r {input.fasta} -v {input.vcf} -S -a -f -p > {output}"


rule split_nodes:
    input:
        "graph/chr10.asm20.vg"
    output:
        "graph/chr10.asm20.chop32.vg"
    shell:
        "vg mod -X 32 {input} > {output}"


rule prune_graph:
    input:
        "graph/chr10.asm20.chop32.vg"
    output:
        "graph/chr10.asm20.chop32.pruned.vg"
    shell:
        "vg prune -r {input} > {output}"


rule index_gcsa:
    input:
        "graph/chr10.asm20.chop32.pruned.vg"
    output:
        "graph/chr10.asm20.gcsa"
    params:
        temp_dir="/data1/temp"
    shell:
        "vg index -g {output} -k 11 -p -b {params.temp_dir} {input}"


rule index_xg:
    input:
        "graph/chr10.asm20.chop32.vg"
    output:
        "graph/chr10.asm20.xg"
    shell:
        "vg index -x {output} {input}"
