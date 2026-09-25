"""Tests for the V5 outgroup analyses (run: python3 tests/test_rooted_analyses.py [--keep DIR]).

Unit tests cover rooting.py, extract_sco5.py, codon_align.py, d_statistics.py,
format_cafe_input.py, normalize_gff.py, check_seqids.py and filter_cds.py.

When snakemake, mafft, trimal and iqtree are on the PATH, the rules
extract_sco5, rooted_alignments, rooted_tree, rooted_summary, d_statistics
and prepare_cafe_input also run on a simulated data set: 80 loci evolved on
((pallens, quinquefasciatus), (molestus, pipiens)) with the outgroup, where
pipiens received a quarter of its loci from the pallens lineage. The rooted
tree must recover the species tree with the root between the two pairs, and
the D statistics must find the pallens and pipiens excess (D < 0 in the
first planned test, D > 0 in the third) and none in the second.
"""
import csv
import os
import random
import runpy
import shutil
import subprocess
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SCRIPTS = os.path.join(ROOT, "workflow", "scripts")
sys.path.insert(0, SCRIPTS)

import rooting  # noqa: E402
import extract_sco5  # noqa: E402
import codon_align  # noqa: E402
import d_statistics  # noqa: E402

QUI, PAL, MOL, PIP, OUT = ("Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus",
                           "Cx_pipiens", "Cx_perexiguus")
INGROUP = [QUI, PAL, MOL, PIP]
SAMPLES = INGROUP + [OUT]


# ---------------------------------------------------------------- unit tests
def test_rooting():
    # IQ-TREE writes unrooted trees with a basal trifurcation on any taxon
    t = rooting.parse("(Cx_molestus:0.01,Cx_pipiens:0.02,((Cx_pallens:0.01,"
                      "Cx_quinquefasciatus:0.01)95/60.1/40.2:0.003,Cx_perexiguus:0.1)100:0.004);")
    clades = rooting.ingroup_clades(t, OUT)
    got = {frozenset(c): (lab, ln) for c, lab, ln in clades}
    assert got == {frozenset([PAL, QUI]): ("95/60.1/40.2", 0.003),
                   frozenset([MOL, PIP]): ("100", 0.004)}, got
    topo = rooting.rooted_topology(list(got), INGROUP)
    assert topo == "((Cx_molestus,Cx_pipiens),(Cx_pallens,Cx_quinquefasciatus));", topo
    assert rooting.sister_pair(list(got), (PAL, QUI, PIP)) == frozenset([PAL, QUI])
    assert rooting.quartet_split(list(got), INGROUP) == \
        "Cx_molestus+Cx_pipiens | Cx_pallens+Cx_quinquefasciatus"
    # the same tree written rooted on the outgroup
    t2 = rooting.parse("(Cx_perexiguus:0.1,(Cx_pallens:0.01,Cx_quinquefasciatus:0.01)95:0.003,"
                       "(Cx_molestus:0.01,Cx_pipiens:0.02)100:0.004);")
    assert rooting.rooted_topology([c for c, _, _ in rooting.ingroup_clades(t2, OUT)],
                                   INGROUP) == topo
    # a caterpillar: the outgroup joins beside pipiens
    t3 = rooting.parse("(Cx_perexiguus:0.1,Cx_pipiens:0.02,((Cx_pallens:0.01,"
                       "Cx_quinquefasciatus:0.01)90:0.003,Cx_molestus:0.02)80:0.002);")
    c3 = [c for c, _, _ in rooting.ingroup_clades(t3, OUT)]
    assert rooting.rooted_topology(c3, INGROUP) == \
        "(((Cx_pallens,Cx_quinquefasciatus),Cx_molestus),Cx_pipiens);"
    assert rooting.sister_pair(c3, (MOL, PIP, QUI)) == frozenset([MOL, QUI])
    assert rooting.nested(c3, INGROUP) == (((PAL, QUI), MOL), PIP)


def test_extract_frames():
    # frame 1 CDS (a partial model starting mid codon) with an internal stop codon
    cds = "G" + "ATGAAA" + "TAA" + "CCCGGG" + "TGA"
    prot, codons = extract_sco5.in_frame(cds, 1)
    assert prot == "MKPG" and codons == ["ATG", "AAA", "CCC", "GGG"]
    ident, frame, trans, cod = extract_sco5.best_frame(cds, "MKPG")
    assert (ident, frame, trans) == (1.0, 1, "MKPG")
    assert extract_sco5.translate_codon("GCN") == "A"
    assert extract_sco5.translate_codon("NNN") == "X"
    assert extract_sco5.translate_codon("RAY") == "X"      # AAC/AAT/GAC/GAT: N or D


def test_codon_align():
    with tempfile.TemporaryDirectory() as tmp:
        sco, aln = os.path.join(tmp, "sco"), os.path.join(tmp, "aln")
        os.makedirs(sco)
        os.makedirs(aln)
        with open(os.path.join(sco, "OG1.faa"), "w") as fh:
            fh.write(">A\nMKP\n>B\nMP\n")
        with open(os.path.join(sco, "OG1.fna"), "w") as fh:
            fh.write(">A\nATGAAACCC\n>B\nATGCCA\n")
        with open(os.path.join(aln, "OG1.aln"), "w") as fh:
            fh.write(">A\nMKP\n>B\nM-P\n")
        with open(os.path.join(aln, "OG1.cols"), "w") as fh:
            fh.write("#ColumnsMap\t0, 2\n")
        # a locus whose alignment is not the translation of its CDS
        with open(os.path.join(sco, "OG2.faa"), "w") as fh:
            fh.write(">A\nMK\n>B\nMK\n")
        with open(os.path.join(sco, "OG2.fna"), "w") as fh:
            fh.write(">A\nATGAAA\n>B\nATGAAA\n")
        with open(os.path.join(aln, "OG2.aln"), "w") as fh:
            fh.write(">A\nMR\n>B\nMK\n")
        with open(os.path.join(aln, "OG2.cols"), "w") as fh:
            fh.write("#ColumnsMap\t0, 1\n")
        out = os.path.join(tmp, "codon")
        tally = codon_align.run(sco, aln, out, os.path.join(tmp, "summary.tsv"))
        assert tally["written"] == 1 and tally["alignment_not_translation"] == 1, tally
        seqs, _ = codon_align.read_fasta(os.path.join(out, "OG1.fna"))
        assert seqs == {"A": "ATGCCC", "B": "ATGCCA"}, seqs


def test_d_counts_and_jackknife():
    # columns (P1, P2, P3, other, O); O = A is ancestral
    cols = ["AGGAA",   # ABBA
            "GAGAA",   # BABA
            "GGAAA",   # BBAA
            "GGGAA",   # all three derived: not counted
            "AGTAA",   # two derived states: not counted
            "AG-AA",   # gap: not a site
            "AAAAA",   # invariant site
            "CGGCC",   # ABBA
            "GAGAA"]   # BABA; third positions are columns 2, 5 (a gap) and 8
    names = ["P1", "P2", "P3", "X", "O"]
    seqs = {n: "".join(c[i] for c in cols) for i, n in enumerate(names)}
    tests = [("T1", "P1", "P2", "P3", True, "planned")]
    counts, n_sites = d_statistics.count_locus(seqs, names, tests, "O")
    allc = counts[("T1", "all_sites")]
    assert (allc["ABBA"], allc["BABA"], allc["BBAA"]) == (2, 2, 1), allc
    assert n_sites["all_sites"] == 8 and n_sites["third_positions"] == 2, n_sites
    third = counts[("T1", "third_positions")]
    assert (third["ABBA"], third["BABA"], third["BBAA"]) == (0, 1, 1), third
    # equal blocks reduce to the ordinary delete one jackknife
    blocks = [(5, 1, 10), (3, 2, 10), (4, 4, 10), (6, 1, 10)]
    d, se, g = d_statistics.block_jackknife(blocks)
    A, B = 18, 8
    loo = [((A - a) - (B - b)) / ((A - a) + (B - b)) for a, b, _ in blocks]
    mean = sum(loo) / 4
    se_ref = ((3 / 4) * sum((x - mean) ** 2 for x in loo)) ** 0.5
    assert abs(d - 10 / 26) < 1e-12 and abs(se - se_ref) < 1e-12 and g == 4, (d, se, se_ref)


def test_plan_tests():
    with tempfile.TemporaryDirectory() as tmp:
        balanced = os.path.join(tmp, "b.nwk")
        with open(balanced, "w") as fh:
            fh.write("(Cx_perexiguus:0.1,(Cx_pallens:0.01,Cx_quinquefasciatus:0.01)95:0.003,"
                     "(Cx_molestus:0.01,Cx_pipiens:0.02)100:0.004);")
        planned = [[PAL, QUI, PIP], [PAL, QUI, MOL], [MOL, PIP, PAL], [MOL, PIP, QUI]]
        tests = d_statistics.plan_tests(balanced, OUT, INGROUP, planned)
        assert [t[1:5] for t in tests] == [(PAL, QUI, PIP, True), (PAL, QUI, MOL, True),
                                           (MOL, PIP, PAL, True), (MOL, PIP, QUI, True)]
        derived = d_statistics.plan_tests(balanced, OUT, INGROUP, [])
        assert sorted(t[1:4] for t in derived) == sorted(tuple(p) for p in planned)
        cater = os.path.join(tmp, "c.nwk")
        with open(cater, "w") as fh:
            fh.write("(Cx_perexiguus:0.1,Cx_pipiens:0.02,((Cx_pallens:0.01,"
                     "Cx_quinquefasciatus:0.01)90:0.003,Cx_molestus:0.02)80:0.002);")
        tests = d_statistics.plan_tests(cater, OUT, INGROUP, planned)
        flags = [(t[1], t[2], t[3], t[4], t[5]) for t in tests]
        assert (MOL, PIP, PAL, False, "planned") in flags
        assert (MOL, PAL, PIP, True, "implied by the rooted tree") in flags


def test_format_cafe_input():
    script = os.path.join(SCRIPTS, "format_cafe_input.py")
    with tempfile.TemporaryDirectory() as tmp:
        of = os.path.join(tmp, "of", "Results_X", "Orthogroups")
        os.makedirs(of)
        with open(os.path.join(of, "Orthogroups.GeneCount.tsv"), "w") as fh:
            fh.write("Orthogroup\t" + "\t".join(SAMPLES) + "\tTotal\n")
            fh.write("OG1\t1\t1\t1\t1\t1\t5\nOG2\t0\t0\t0\t0\t3\t3\n"
                     "OG3\t101\t1\t1\t1\t0\t104\nOG4\t2\t0\t1\t0\t0\t3\n")
        rooted = os.path.join(tmp, "rooted.nwk")
        unrooted = os.path.join(tmp, "concat.nwk")
        with open(unrooted, "w") as fh:
            fh.write("(Cx_molestus:0.01,Cx_pipiens:0.01,(Cx_pallens:0.01,"
                     "Cx_quinquefasciatus:0.01)100:0.004);")
        for text, expect in (
                ("(Cx_perexiguus:0.1,(Cx_pallens:0.01,Cx_quinquefasciatus:0.01)95:0.003,"
                 "(Cx_molestus:0.01,Cx_pipiens:0.02)100:0.004);",
                 "((Cx_molestus:1.0,Cx_pipiens:1.0):0.5,(Cx_pallens:1.0,"
                 "Cx_quinquefasciatus:1.0):0.5);"),
                ("(Cx_perexiguus:0.1,Cx_pipiens:0.02,((Cx_pallens:0.01,"
                 "Cx_quinquefasciatus:0.01)90:0.003,Cx_molestus:0.02)80:0.002);",
                 "(((Cx_pallens:1.0,Cx_quinquefasciatus:1.0):1.0,Cx_molestus:2.0):1.0,"
                 "Cx_pipiens:3.0);")):
            with open(rooted, "w") as fh:
                fh.write(text)
            sm = types.SimpleNamespace(
                input=types.SimpleNamespace(counts=os.path.join(tmp, "of"), tree=unrooted,
                                            rooted=rooted),
                params=types.SimpleNamespace(ingroup=INGROUP, outgroup=[OUT]),
                output=types.SimpleNamespace(counts=os.path.join(tmp, "counts.tsv"),
                                             tree=os.path.join(tmp, "tree.nwk")))
            runpy.run_path(script, init_globals={"snakemake": sm}, run_name="__main__")
            tree = open(sm.output.tree).read().strip()
            assert tree == expect, tree
        rows = list(csv.DictReader(open(sm.output.counts), delimiter="\t"))
        assert [r["Family ID"] for r in rows] == ["OG1", "OG4"], rows
        assert list(rows[0]) == ["Family ID", "Desc"] + INGROUP


def test_gff_tools():
    with tempfile.TemporaryDirectory() as tmp:
        gff = os.path.join(tmp, "in.gff3")
        with open(gff, "w") as fh:
            fh.write("##gff-version 3\n"
                     "chr1\tens\tgene\t1\t90\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n"
                     "chr1\tens\tmRNA\t1\t90\t.\t+\t.\tID=transcript:T1.1;Parent=gene:G1\n"
                     "chr1\tens\texon\t1\t90\t.\t+\t.\tParent=transcript:T1.1;Name=E1\n"
                     "chr1\tens\tCDS\t1\t90\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1.1\n"
                     "chrX\tens\tgene\t1\t9\t.\t+\t.\tID=gene:G2\n")
        out = os.path.join(tmp, "out.gff3")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "normalize_gff.py"), gff, out],
                       check=True, capture_output=True)
        text = open(out).read()
        assert "ID=T1.1;Parent=G1" in text and "ID=P1;Parent=T1.1" in text, text
        assert "gene:" not in text and "transcript:" not in text
        fa = os.path.join(tmp, "g.fa")
        with open(fa, "w") as fh:
            fh.write(">chr1 dna:chromosome\n" + "A" * 100 + "\n")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "check_seqids.py"), fa, out],
                           capture_output=True, text=True)
        assert r.returncode == 1 and "differently" in r.stderr, r   # chrX is 1 of 5 features
        with open(out, "a") as fh:
            for i in range(200):
                fh.write(f"chr1\tens\texon\t1\t9\t.\t+\t.\tParent=T1.1;Name=X{i}\n")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "check_seqids.py"), fa, out],
                           capture_output=True, text=True)
        assert r.returncode == 0, r
        # filter_cds keeps the transcripts of the protein file, periods stripped
        cds = os.path.join(tmp, "cds.fa")
        with open(cds, "w") as fh:
            fh.write(">T1.1 gene=G1\natgaaa\n>T9.1\nATGCCC\n")
        prot = os.path.join(tmp, "p.fa")
        with open(prot, "w") as fh:
            fh.write(">T11\nMK\n")
        res = os.path.join(tmp, "cds_out.fa")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "filter_cds.py"), "--cds", cds,
                        "--proteins", prot, "--out", res, "--strip-periods"],
                       check=True, capture_output=True)
        assert open(res).read() == ">T11\nATGAAA\n"


# ---------------------------------------------------------------- simulation
CODONS = [c for c, a in extract_sco5.CODON_TABLE.items() if a != "*"]


def evolve(seq, rate, rng):
    out = []
    for codon in seq:
        c = list(codon)
        for i in range(3):
            if rng.random() < rate:
                c[i] = rng.choice([b for b in "ACGT" if b != c[i]])
        new = "".join(c)
        out.append(codon if extract_sco5.CODON_TABLE[new] == "*" else new)
    return out


def simulate_locus(rng, n_codons, introgressed):
    anc = ["ATG"] + [rng.choice(CODONS) for _ in range(n_codons - 1)]
    out = evolve(anc, 0.08, rng)
    ing = evolve(anc, 0.004, rng)
    a1, a2 = evolve(ing, 0.006, rng), evolve(ing, 0.006, rng)
    qui, mol = evolve(a1, 0.008, rng), evolve(a2, 0.008, rng)
    if introgressed:
        p = evolve(a1, 0.004, rng)
        pal, pip = evolve(p, 0.004, rng), evolve(p, 0.004, rng)
    else:
        pal, pip = evolve(a1, 0.008, rng), evolve(a2, 0.008, rng)
    return {QUI: qui, PAL: pal, MOL: mol, PIP: pip, OUT: out}


def translate(codons):
    return "".join(extract_sco5.CODON_TABLE[c] for c in codons)


def build_workspace(ws, n_loci=80, seed=7):
    """A workflow directory with the inputs of the outgroup rules."""
    rng = random.Random(seed)
    for name in ("Snakefile", "workflow"):
        os.symlink(os.path.join(ROOT, name), os.path.join(ws, name))
    os.makedirs(os.path.join(ws, "config"))
    shutil.copy(os.path.join(ROOT, "config", "config.yaml"), os.path.join(ws, "config"))
    shutil.copy(os.path.join(ROOT, "config", "samples.tsv"), os.path.join(ws, "config"))
    for d in ("results/proteins", "results/cds", "results/annotation", "results/phylo",
              "results/orthofinder/output/Results_Sim/Orthogroups"):
        os.makedirs(os.path.join(ws, d))
    prot = {s: [] for s in SAMPLES}
    cds = {s: [] for s in SAMPLES}
    gff = {s: ["##gff-version 3"] for s in INGROUP}
    members, counts = [], []
    for i in range(n_loci):
        og = f"OG{i:07d}"
        loc = simulate_locus(rng, rng.randint(150, 300), introgressed=i % 4 == 0)
        ids = {}
        for s in SAMPLES:
            if s == OUT:
                tid = f"CPERX{i:05d}T1"
            else:
                tid = f"rna-XM_{i:09d}1"
            ids[s] = tid
            nt = "".join(loc[s])
            aa = translate(loc[s])
            if i == 1 and s == PAL:        # partial model: CDS starts mid codon
                nt = "C" + nt
            if i == 2 and s == MOL:        # broken model: in frame stop codon
                nt = nt[:30] + "TAG" + nt[30:]
            if i == 3 and s == PIP:        # protein that is not the CDS translation
                aa = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in aa)
            prot[s].append(f">{tid}\n{aa}")
            cds[s].append(f">{tid}\n{nt}")
            if s != OUT:
                chrom = "NC_1" if i < n_loci // 2 else "NC_2"
                start = 1 + (i % (n_loci // 2)) * 1_300_000
                valid = "False" if (i in (2, 5, 9) and s != QUI) else "True"
                attrs = f"ID={tid[:-1]}.{tid[-1]};Parent=gene-G{i}"
                if s != QUI:
                    attrs += f";valid_ORF={valid}"
                gff[s].append(f"{chrom}\tsim\tmRNA\t{start}\t{start + 900}\t.\t+\t.\t{attrs}")
        members.append([og] + [ids[s] for s in SAMPLES])
        counts.append([og] + ["1"] * len(SAMPLES))
    # a two copy orthogroup and an outgroup only one: neither is used
    members.append(["OGX1", "rna-XM_9990000011, rna-XM_9990000021", "", "", "", ""])
    counts.append(["OGX1", "2", "0", "0", "0", "0"])
    ofd = os.path.join(ws, "results/orthofinder/output/Results_Sim/Orthogroups")
    for name, rows in (("Orthogroups.tsv", members), ("Orthogroups.GeneCount.tsv", counts)):
        with open(os.path.join(ofd, name), "w") as fh:
            fh.write("Orthogroup\t" + "\t".join(SAMPLES) + ("\tTotal" if "Count" in name else "")
                     + "\n")
            for r in rows:
                extra = [str(sum(int(x) for x in r[1:]))] if "Count" in name else []
                fh.write("\t".join(r + extra) + "\n")
    for s in SAMPLES:
        with open(os.path.join(ws, f"results/proteins/{s}.fa"), "w") as fh:
            fh.write("\n".join(prot[s]) + "\n")
        with open(os.path.join(ws, f"results/cds/{s}.fa"), "w") as fh:
            fh.write("\n".join(cds[s]) + "\n")
    for s in INGROUP:
        with open(os.path.join(ws, f"results/annotation/{s}_liftoff.gff3"), "w") as fh:
            fh.write("\n".join(gff[s]) + "\n")
    with open(os.path.join(ws, "results/phylo/concat_tree.treefile"), "w") as fh:
        fh.write("(Cx_molestus:0.01,Cx_pipiens:0.01,(Cx_pallens:0.01,"
                 "Cx_quinquefasciatus:0.01)100:0.004);\n")


def run_rules(ws):
    targets = ["results/phylo/dstat/d_statistics.tsv", "results/phylo/rooted/rooted_summary.tsv",
               "results/cafe/ultrametric_tree.nwk"]
    rules = ["extract_sco5", "rooted_alignments", "rooted_tree", "rooted_summary",
             "d_statistics", "prepare_cafe_input"]
    cmd = ["snakemake", "--cores", "2", "--allowed-rules", *rules, "--rerun-triggers", "mtime",
           "--nolock", *targets]
    r = subprocess.run(cmd, cwd=ws, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-5000:], r.stderr[-8000:])
        raise SystemExit("snakemake failed on the simulated data")
    return r


def test_workflow(keep=None):
    tools = ["snakemake", "mafft", "trimal", "iqtree"]
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        print(f"workflow test skipped: {', '.join(missing)} not on the PATH")
        return
    ws = keep or tempfile.mkdtemp(prefix="v5sim_")
    if keep:
        os.makedirs(ws, exist_ok=True)
    build_workspace(ws)
    run_rules(ws)
    acc = {(r["item"], r["sample"]): r["value"] for r in
           csv.DictReader(open(os.path.join(ws, "results/phylo/rooted/sco5_accounting.tsv")),
                          delimiter="\t")}
    assert acc[("single_copy_in_all_taxa", "")] == "80", acc
    assert acc[("removed_cds_mismatch", "")] == "1", acc
    assert acc[("written", "")] == "79", acc
    assert acc[("cds_translation_frame_not_0", PAL)] == "1", acc
    summ = dict(csv.reader(open(os.path.join(ws, "results/phylo/rooted/rooted_summary.tsv")),
                           delimiter="\t"))
    assert summ["rooted ingroup topology"] == \
        "((Cx_molestus,Cx_pipiens),(Cx_pallens,Cx_quinquefasciatus));", summ
    assert summ["root position"] == ("between (Cx_molestus,Cx_pipiens) and "
                                     "(Cx_pallens,Cx_quinquefasciatus)"), summ
    assert summ["(Cx_pallens,Cx_quinquefasciatus) gCF"] not in ("NA", ""), summ
    rows = {r["analysis"]: r for r in
            csv.DictReader(open(os.path.join(ws, "results/phylo/dstat/d_statistics.tsv")),
                           delimiter="\t")}
    t1, t2, t3 = (rows[f"T{i} all_loci all_sites"] for i in (1, 2, 3))
    assert t1["statistic"] == "D(Cx_pallens,Cx_quinquefasciatus;Cx_pipiens,Cx_perexiguus)"
    assert float(t1["D"]) < 0 and float(t1["Z"]) <= -3, t1
    assert float(t3["D"]) > 0 and float(t3["Z"]) >= 3, t3
    assert abs(float(t2["Z"])) < 3, t2
    assert t1["excess_derived_sharing"] == "Cx_pallens+Cx_pipiens"
    assert "T1 intact_loci third_positions" in rows
    intact = rows["T1 intact_loci all_sites"]
    assert int(intact["n_loci"]) == 76, intact          # loci 2, 5 and 9 are not intact
    tree = open(os.path.join(ws, "results/cafe/ultrametric_tree.nwk")).read().strip()
    assert tree == ("((Cx_molestus:1.0,Cx_pipiens:1.0):0.5,"
                    "(Cx_pallens:1.0,Cx_quinquefasciatus:1.0):0.5);"), tree
    print(f"workflow test passed in {ws}")
    if not keep:
        shutil.rmtree(ws, ignore_errors=True)


def main():
    keep = sys.argv[sys.argv.index("--keep") + 1] if "--keep" in sys.argv else None
    test_rooting()
    test_extract_frames()
    test_codon_align()
    test_d_counts_and_jackknife()
    test_plan_tests()
    test_format_cafe_input()
    test_gff_tools()
    print("rooted analysis unit tests passed")
    test_workflow(keep)


if __name__ == "__main__":
    sys.exit(main())
