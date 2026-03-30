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
        # Download
        datasets download genome accession {params.accession} \
            --include genome,gff3 --filename {wildcards.sample}.zip

        # Unzip
        unzip -o {wildcards.sample}.zip -d tmp_{wildcards.sample}

        # Move genome FASTA
        find tmp_{wildcards.sample} -name "*.fna" | head -1 | \
            xargs -I PLACEHOLDER mv PLACEHOLDER {output.fasta}

        # Move GFF if it exists, otherwise create empty file
        GFF_FILE=$(find tmp_{wildcards.sample} -name "*.gff" | head -1)
        if [ -n "$GFF_FILE" ]; then
            mv "$GFF_FILE" {output.gff}
        else
            touch {output.gff}
        fi

        # Cleanup
        rm -rf tmp_{wildcards.sample} {wildcards.sample}.zip
        """
