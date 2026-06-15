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
        if [ "{wildcards.sample}" = "Cx_tarsalis" ]; then
            # The CtarK1 assembly (Main et al. 2021) is hosted on OSF, not NCBI.
            # The original NCBI accession GCA_016859205.1 was reassigned to a
            # Fusarium oxysporum genome — do NOT use the datasets CLI for tarsalis.
            curl -L -o /tmp/{wildcards.sample}.fa.gz   https://osf.io/download/3dgxa/
            curl -L -o /tmp/{wildcards.sample}.gff3.gz https://osf.io/download/rj97c/
            gunzip -c /tmp/{wildcards.sample}.fa.gz   > {output.fasta}
            gunzip -c /tmp/{wildcards.sample}.gff3.gz > {output.gff}
            rm /tmp/{wildcards.sample}.fa.gz /tmp/{wildcards.sample}.gff3.gz
            exit 0
        fi
        # Default path: NCBI Datasets CLI
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
