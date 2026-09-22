"""Tests for workflow/scripts/cloud_composition.py (run: python3 tests/test_cloud_composition.py)."""
import importlib.util
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "cloud_composition", os.path.join(HERE, "..", "workflow", "scripts", "cloud_composition.py"))
cc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cc)

INGROUP = ["Q", "M", "P"]
# OG1: M only plus outgroup (transfer of gene X split from its source in OG2)
# OG3: Q only plus outgroup, gene Y has no model anywhere else
# OG4: M only, two in-paralogs, no outgroup, gene Z also modelled in P (OG5)
MEMBERS = pd.DataFrame({
    "Q": ["", "rna-XM_11", "rna-XM_21", "", ""],
    "M": ["rna-XM_11", "", "", "rna-XM_31, rna-XM_32", ""],
    "P": ["", "rna-XM_11", "", "", "rna-XM_31"],
    "T": ["t1", "t2", "t3", "", ""],
}, index=["OG1", "OG2", "OG3", "OG4", "OG5"])
TABLE = pd.DataFrame({
    "Q": [0, 1, 1, 0, 0], "M": [1, 0, 0, 2, 0], "P": [0, 1, 0, 0, 1],
    "T": [1, 1, 1, 0, 0],
    "compartment": ["cloud", "shell", "cloud", "cloud", "cloud"],
}, index=MEMBERS.index)
TX = {"Q": {"rna-XM_11": "gene-X", "rna-XM_21": "gene-Y"},
      "M": {"rna-XM_11": "gene-X", "rna-XM_31": "gene-Z", "rna-XM_32": "gene-Z2"},
      "P": {"rna-XM_11": "gene-X", "rna-XM_31": "gene-Z"}}
CODING = {"Q": {"rna-XM_11", "rna-XM_21"},
          "M": {"rna-XM_11", "rna-XM_31", "rna-XM_32"},
          "P": {"rna-XM_11", "rna-XM_31"}}


def main():
    per = cc.compose(TABLE, MEMBERS, TX, CODING, INGROUP, "Q").set_index("orthogroup")
    og1 = per.loc["OG1"]
    assert og1.form == "M" and og1.n_outgroup_genes == 1 and og1.n_ingroup_genes == 1
    assert og1.status == "same_gene_elsewhere" and og1.in_reference_elsewhere
    assert og1.same_gene_orthogroups == "OG2" and og1.same_gene_forms == "P,Q"
    assert og1.same_gene_compartments == "shell"
    og3 = per.loc["OG3"]
    assert og3.is_reference_form and og3.status == "no_counterpart"
    og4 = per.loc["OG4"]
    assert og4.n_outgroup_genes == 0 and og4.status == "same_gene_elsewhere"
    assert not og4.in_reference_elsewhere and og4.same_gene_forms == "P"
    s = cc.summarise(per.reset_index(), INGROUP).set_index("form")
    assert s.loc["total", "n_cloud_orthogroups"] == 4
    assert s.loc["total", "with_outgroup_gene"] == 2
    assert s.loc["M", "same_gene_elsewhere"] == 2 and s.loc["M", "same_gene_in_reference"] == 1
    assert s.loc["Q", "no_counterpart"] == 1
    print("cloud_composition: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
