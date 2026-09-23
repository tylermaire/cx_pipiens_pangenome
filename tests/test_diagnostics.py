"""Tests for divergence_diagnostics.py and transfer_quality.py (run: python3 tests/test_diagnostics.py)."""
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, "..", "workflow", "scripts", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


dd = load("divergence_diagnostics")
tq = load("transfer_quality")


def main():
    concat = "(A:0.04,(B:0.03,C:0.05)100:0.013,D:0.038);"
    genes = ["(A:0.004,(B:0.003,C:0.006)97:0.005,D:0.0000010000);",
             "(A:0.30,(B:0.002,D:0.006)40:0.002,C:0.005);",
             "(A:0.002,(C:0.004,D:0.004)88:0.001,B:0.001);"]
    rows, n = dd.branch_table(concat, genes)
    by = {r["branch"]: r for r in rows}
    assert n == 3 and by["A"]["concat_tree"] == 0.04
    assert by["A"]["gene_tree_median"] == 0.004 and by["A"]["share_above_0.1"] == round(1 / 3, 4)
    assert by["D"]["share_minimum"] == round(1 / 3, 4)
    assert by["internal"]["concat_tree"] == 0.013 and by["internal"]["n_gene_trees"] == 3

    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "OG1.trim"), "w") as fh:
            fh.write(">A\nMKLV-A\n>B\nMKLVQA\n>C\nMRLVQA\n")
        with open(os.path.join(tmp, "OG2.trim"), "w") as fh:
            fh.write(">A\nMKK\n>B\nMKK\n>C\nMKK\n")
        pair_rows, invariant, few = dd.pairwise_table(sorted(
            os.path.join(tmp, f) for f in os.listdir(tmp)))
    ab = next(r for r in pair_rows if (r["taxon_a"], r["taxon_b"]) == ("A", "B"))
    assert ab["n_loci"] == 2 and ab["median_identity"] == 1.0
    ac = next(r for r in pair_rows if (r["taxon_a"], r["taxon_b"]) == ("A", "C"))
    assert ac["share_below_0.95"] == 0.5          # OG1: 4 of 5 shared columns
    assert invariant == {"OG2"} and few == {"OG1", "OG2"}

    import tempfile as _t
    with _t.TemporaryDirectory() as tmp:
        gff = os.path.join(tmp, "M.gff3")
        with open(gff, "w") as fh:
            fh.write("##gff-version 3\n")
            for tid, attrs in (("rna-1.1", "valid_ORF=True;matches_ref_protein=True"),
                               ("rna-2.1", "valid_ORF=False;inframe_stop_codon=True;matches_ref_protein=False"),
                               ("rna-3.1", "valid_ORF=False;missing_start_codon=True;missing_stop_codon=True"),
                               ("rna-4.1", "valid_ORF=True")):
                fh.write(f"c1\tLiftoff\tmRNA\t1\t9\t.\t+\t.\tID={tid};Parent=g;{attrs}\n")
        refgff = os.path.join(tmp, "Q.gff3")
        with open(refgff, "w") as fh:
            for tid, attrs in (("rna-1.1", "gbkey=mRNA"),
                               ("rna-2.1", "exception=unclassified transcription discrepancy"),
                               ("rna-3.1", "partial=true;start_range=.,1"),
                               ("rna-4.1", "gbkey=mRNA")):
                fh.write(f"c1\tGnomon\tmRNA\t1\t9\t.\t+\t.\tID={tid};Parent=g;{attrs}\n")
        flags = tq.liftoff_flags(gff)
        ref_flags = tq.liftoff_flags(refgff)
    kept = ["rna-11", "rna-21", "rna-31", "rna-41"]
    assert flags["rna-21"]["inframe_stop_codon"] == "True"
    row = tq.summarise_sample("M", kept, flags, False, ref_flags)
    assert row["n_invalid_orf"] == 2 and row["pct_invalid_orf"] == 50.0
    assert row["n_inframe_stop"] == 1 and row["n_missing_start"] == 1
    assert row["n_mismatch_ref_protein"] == 1 and row["liftoff_flags"] == "yes"
    # both broken transfers come from reference models that are not clean ORFs
    assert row["n_ref_model_not_clean"] == 2 and row["n_clean_ref_models"] == 2
    assert row["n_invalid_orf_clean_ref"] == 0 and row["pct_invalid_orf_clean_ref"] == 0.0
    ref = tq.summarise_sample("Q", kept, ref_flags, True, ref_flags)
    assert ref["liftoff_flags"].startswith("no (reference") and ref["pct_invalid_orf"] == ""
    assert ref["n_ref_model_not_clean"] == 2 and ref["pct_invalid_orf_clean_ref"] == ""
    where = {("M", "rna-11"): "core", ("M", "rna-21"): "cloud", ("M", "rna-31"): "cloud"}
    comp = {r["compartment"]: r for r in tq.by_compartment("M", kept, flags, where)}
    assert comp["cloud"]["pct_invalid_orf"] == 100.0 and comp["core"]["n_invalid_orf"] == 0
    assert comp["unassigned"]["n_models"] == 1 and "shell" not in comp
    where[("M", "rna-41")] = "outgroup_only"
    comp = {r["compartment"]: r for r in tq.by_compartment("M", kept, flags, where)}
    assert comp["outgroup_only"]["n_models"] == 1 and "unassigned" not in comp
    print("diagnostics: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
