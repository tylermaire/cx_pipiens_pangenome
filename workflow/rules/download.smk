rule download_genome:
    """Download genome assembly from NCBI."""
    output:
        fasta="resources/genomes/{sample}.fasta",
        gff="resources/annotations/{sample}.gff3"
    params:
        accession=lambda wc: samples.loc[wc.sample, "accession"]
    conda: "../envs/datasets.yaml"
    shell:
        """
        datasets download genome accession {params.accession} \
            --include genome,gff3 --filename {wildcards.sample}.zip
        unzip -o {wildcards.sample}.zip -d tmp_{wildcards.sample}
        find tmp_{wildcards.sample} -name "*.fna" | head -1 | \
            xargs -I{} mv {} {output.fasta}
        find tmp_{wildcards.sample} -name "*.gff" | head -1 | \
            xargs -I{} mv {} {output.gff} || touch {output.gff}
        rm -rf tmp_{wildcards.sample} {wildcards.sample}.zip
        """
