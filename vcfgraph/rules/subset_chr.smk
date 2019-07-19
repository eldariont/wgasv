rule download_assemblies:
    output:
        "assemblies/{sample}_shasta_marginpolish_helen_consensus.fa"
    shell:
        "wget -P assemblies https://storage.googleapis.com/kishwar-helen/polished_genomes/london_calling_2019/CHM13_shasta_marginpolish_helen_consensus.fa"


rule align_assemblies_to_ref:
    input:
        amb="assemblies/{sample}_shasta_marginpolish_helen_consensus.fa",
        ref="assemblies/hg38.fa"
    output:
        "assemblies_aligned_hg38/{sample}_hg38.paf"
    threads: 3
    shell:
        "minimap2 -t {threads} -c --cs -x asm5 {input.ref} {input.amb} > {output}"


rule subset_paf:
    input:
        paf="assemblies_aligned_hg38/{sample}_hg38.paf",
    output:
        pdf="assemblies_{chr}/{sample}_paf_{chr}.pdf",
        contigs="assemblies_{chr}/contigs_{sample}_paf_{chr}.txt",
    shell:
        "Rscript scripts/pafFilterAmb.R {input.paf} {output.pdf} {wildcards.chr} {output.contigs}"


rule subset_assemblies:
    input:
        fasta="assemblies/{sample}_shasta_marginpolish_helen_consensus.fa",
        ids="assemblies_{chr}/contigs_{sample}_paf_{chr}.txt"
    output:
        "assemblies_{chr}/{sample}_paf_{chr}.fa"
    shell:
        "python scripts/filterFastaByIds.py -f {input.fasta} -i {input.ids} -o {output}"


rule subset_reference:
    input:
        "assemblies/hg38.fa"
    output:
        "assemblies_{chr}/hg38_{chr}.fa"
    shell:
        "samtools faidx {input} {chr} > {output}"