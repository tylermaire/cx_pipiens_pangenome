"""Tests for workflow/scripts/classify_absence.py (run: python3 tests/test_classify_absence.py)."""
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "classify_absence", os.path.join(HERE, "..", "workflow", "scripts", "classify_absence.py"))
ca = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ca)

GFF = {
    # reference: genes X, Y, W, Z
    "Q": [("chr1", 100, 900, "rna-XM_1.1", "gene-X"),
          ("chr1", 2000, 2900, "rna-XM_2.1", "gene-Y"),
          ("chr1", 4000, 4900, "rna-XM_4.1", "gene-W")],
    # molestus: broken transfer of X; paralog Z where Y's probe lands;
    # model V near identical to W where W's probe lands; a non coding transcript
    "M": [("chr1", 100, 900, "rna-XM_1.1", "gene-X"),
          ("chr1", 2000, 2900, "rna-XM_3.1", "gene-Z"),
          ("chr1", 4000, 4900, "rna-XM_5.1", "gene-V"),
          ("chr1", 6000, 6900, "rna-XR_9.1", "gene-N")],
    # pallens: faithful transfer of X; nothing at Y's locus; nothing for W
    "P": [("chr1", 100, 900, "rna-XM_1.1", "gene-X")],
}
PROT = {
    "Q": {"rna-XM_11": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ",
          "rna-XM_21": "MSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSEFELRRQACGRTRAPPPPPLRSGC",
          "rna-XM_41": "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLLVKHSQ"},
    "M": {"rna-XM_11": "MKTAYIAKQRQISFVKSHFSRQLEWGVRGQNLASKHLSASSLQRDADGNTQDKDLAVEKALRCQ",
          "rna-XM_31": "MDNLKQVVEQALLGRSKLTEAHWPDISYVAGNRGKIVMLDNGQQTFEPVTHAQ",
          "rna-XM_51": "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLLVKHSQ"},
    "P": {"rna-XM_11": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"},
}
MEMBERS = {("M", "rna-XM_11"): "OG1",
           ("Q", "rna-XM_11"): "OG2", ("P", "rna-XM_11"): "OG2",
           ("Q", "rna-XM_21"): "OG3", ("M", "rna-XM_31"): "OG5",
           ("Q", "rna-XM_41"): "OG4", ("M", "rna-XM_51"): "OG6"}
# probe hits: (identity, coverage, seqid, start, end)
PROT_HITS = {("OG3", "M"): (0.85, 0.9, "chr1", 2050, 2850),
             ("OG3", "P"): (0.90, 0.9, "chr1", 2050, 2850),
             ("OG4", "M"): (0.97, 1.0, "chr1", 4010, 4890)}
DNA_HITS = {}


def write_gffs(tmp):
    out = {}
    for form, feats in GFF.items():
        path = os.path.join(tmp, f"{form}.gff3")
        with open(path, "w") as fh:
            fh.write("##gff-version 3\n")
            for seqid, s, e, tid, gid in feats:
                kind = "transcript" if tid.startswith("rna-XR") else "mRNA"
                fh.write(f"{seqid}\tLiftoff\t{kind}\t{s}\t{e}\t.\t+\t.\t"
                         f"ID={tid};Parent={gid};gbkey=mRNA\n")
        out[form] = ca.read_gff_transcripts(path)
    return out


def main():
    with tempfile.TemporaryDirectory() as tmp:
        gffs = write_gffs(tmp)
    assert gffs["M"][1]["rna-XM_11"] == "gene-X", "periods must be stripped from IDs"
    clf = ca.Classifier(["Q", "M", "P"], PROT, gffs, MEMBERS, PROT_HITS, DNA_HITS)

    # broken transfer of X in M: the reference and pallens hold the same gene
    r = clf.classify("OG1", "cloud", "rna-XM_11", "M", "Q")
    assert (r["call"], r["match_basis"], r["other_orthogroup"]) == \
        ("clustered_elsewhere", "same_gene_id", "OG2"), r
    assert r["best_protein_identity"] < 0.9, "divergent transfer must still count as the same gene"
    r = clf.classify("OG1", "cloud", "rna-XM_11", "M", "P")
    assert (r["call"], r["match_basis"]) == ("clustered_elsewhere", "same_gene_id"), r

    # Y absent from M: the probe lands on a distant paralog
    r = clf.classify("OG3", "cloud", "rna-XM_21", "Q", "M")
    assert r["call"] == "paralog_only" and r["best_overlapping_gene"] == "rna-XM_31", r
    # Y in P: sequence found, no coding model there
    r = clf.classify("OG3", "cloud", "rna-XM_21", "Q", "P")
    assert r["call"] == "unannotated_locus", r

    # W: probe lands on a near identical model of another gene ID
    r = clf.classify("OG4", "cloud", "rna-XM_41", "Q", "M")
    assert (r["call"], r["match_basis"], r["other_orthogroup"]) == \
        ("clustered_elsewhere", "overlap_identity", "OG6"), r
    # W in P: nothing found at all
    r = clf.classify("OG4", "cloud", "rna-XM_41", "Q", "P")
    assert r["call"] == "absent", r

    # a hit on a non coding transcript only is not annotation
    assert clf.overlapping("M", "chr1", 6100, 6800) == []

    import pandas as pd
    calls = pd.DataFrame([
        clf.classify("OG1", "cloud", "rna-XM_11", "M", "Q"),
        clf.classify("OG3", "cloud", "rna-XM_21", "Q", "M"),
        clf.classify("OG3", "cloud", "rna-XM_21", "Q", "P"),
        clf.classify("OG4", "cloud", "rna-XM_41", "Q", "P")])
    s = ca.summarise(calls).iloc[0]
    assert s.n_events == 4 and s.absent == 1 and s.paralog_only == 1
    assert s.clustered_elsewhere_same_gene_id == 1 and s.unannotated_locus == 1
    assert s.pct_supported == 50.0 and s.pct_artifact == 50.0
    print("classify_absence: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
