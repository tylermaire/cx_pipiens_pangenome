rule download_genome:
    """Download a genome assembly (and its GFF, when NCBI has one) with the
    NCBI Datasets CLI."""
    output:
        fasta="resources/genomes/{sample}.fasta",
        gff="resources/annotations/{sample}.gff3"
    wildcard_constraints:
        sample=one_of(NCBI_SAMPLES)
    params:
        accession=lambda wc: samples.loc[wc.sample, "accession"]
    conda: "../envs/datasets.yaml"
    shell:
        """
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


rule download_ensembl:
    """Genome and gene set of a sample annotated by Ensembl (the outgroup).
    Both come from the Ensembl FTP so that the GFF sequence names match the
    FASTA. The genome is the unmasked FASTA under genome/ (softmasked if that
    is all there is). The gene set is the genes.gff3.gz of the newest release
    under <provider>/geneset/, the provider being ensembl when present (the
    directory name Ensembl uses for its own gene builds), else the first
    other one listed (for example braker). config outgroup genome_url and
    gff_url override either file. check_seqids.py stops the run if the two
    do not match."""
    output:
        fasta="resources/genomes/{sample}.fasta",
        gff="resources/annotations/{sample}.gff3"
    wildcard_constraints:
        sample=one_of(ENSEMBL_SAMPLES)
    params:
        base=config["outgroup"].get("ensembl_base", ""),
        genome_url=config["outgroup"].get("genome_url", "") or "",
        gff_url=config["outgroup"].get("gff_url", "") or ""
    conda: "../envs/phylo.yaml"
    shell:
        """
        set -euo pipefail
        BASE="{params.base}"
        GENOME_URL="{params.genome_url}"
        GFF_URL="{params.gff_url}"
        links() {{
            curl -fsSL "$1" | grep -oE 'href="[^"?/][^"]*"' | sed -e 's/^href="//' -e 's/"$//'
        }}
        if [ -z "$GENOME_URL" ]; then
            FILES=$(links "$BASE/genome/" || true)
            F=$(printf '%s\n' "$FILES" | grep -E '[.]fa[.]gz$' | grep -m1 'unmasked' || true)
            if [ -z "$F" ]; then
                F=$(printf '%s\n' "$FILES" | grep -E '[.]fa[.]gz$' | grep -m1 'softmasked' || true)
            fi
            if [ -z "$F" ]; then
                echo "No genome FASTA found under $BASE/genome/. Open that address in a" >&2
                echo "browser and put the FASTA's full URL in config.yaml outgroup genome_url." >&2
                exit 1
            fi
            GENOME_URL="$BASE/genome/$F"
        fi
        if [ -z "$GFF_URL" ]; then
            PROVIDERS=$(links "$BASE/" | grep -E '/$' | sed 's#/$##' | grep -v -x -E 'genome|vep|variation' || true)
            ORDERED=$(printf '%s\n' "$PROVIDERS" | grep -x ensembl || true)
            ORDERED="$ORDERED $(printf '%s\n' "$PROVIDERS" | grep -v -x ensembl || true)"
            for P in $ORDERED; do
                REL=$(links "$BASE/$P/geneset/" 2>/dev/null | grep -oE '^[0-9]{{4}}_[0-9]{{2}}' | sort -u | tail -1 || true)
                [ -n "$REL" ] || continue
                F=$(links "$BASE/$P/geneset/$REL/" | grep -E 'genes[.]gff3[.]gz$' | head -1 || true)
                GFF_URL="$BASE/$P/geneset/$REL/${{F:-genes.gff3.gz}}"
                break
            done
            if [ -z "$GFF_URL" ]; then
                echo "No gene set release found under $BASE/<provider>/geneset/. Open $BASE/" >&2
                echo "in a browser and put the full URL of the gene set GFF3 in config.yaml" >&2
                echo "outgroup gff_url." >&2
                exit 1
            fi
        fi
        echo "genome: $GENOME_URL"
        echo "genes:  $GFF_URL"
        mkdir -p resources/genomes resources/annotations
        curl -fsSL "$GENOME_URL" | gunzip -c > {output.fasta}
        curl -fsSL "$GFF_URL" | gunzip -c > {output.gff}
        printf '%s\n%s\n' "$GENOME_URL" "$GFF_URL" > resources/annotations/{wildcards.sample}.source.txt
        python workflow/scripts/check_seqids.py {output.fasta} {output.gff}
        """
