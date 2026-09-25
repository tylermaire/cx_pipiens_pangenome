#!/usr/bin/env python3
"""rooting.py - small Newick helpers for the rooted (outgroup) analyses.

A tree from IQ-TREE is unrooted, however it is written. Each internal branch
splits the taxa in two, and with an outgroup the side without the outgroup is
an ingroup clade of the rooted tree. Working from these splits, rather than
re-rooting a parsed tree, keeps branch labels (bootstrap values, concordance
branch ids) on the branch they belong to.

    parse(text)                  -> Node (children, name, label, length)
    ingroup_clades(tree, og)     -> [(frozenset clade, label, length)] per internal branch
    nested(clades, ingroup)      -> the rooted ingroup tree as nested tuples
    rooted_topology(clades, ing) -> canonical Newick of the rooted ingroup tree
    sister_pair(clades, triple)  -> the two taxa of triple that are sisters, or None
    quartet_split(clades, ing)   -> 'a+b | c+d', the unrooted ingroup split
"""
import re


class Node:
    __slots__ = ("children", "name", "label", "length")

    def __init__(self):
        self.children, self.name, self.label, self.length = [], None, None, None

    def tips(self):
        if not self.children:
            return [self.name]
        out = []
        for c in self.children:
            out.extend(c.tips())
        return out

    def internal(self):
        """Every internal node below this one (not this one)."""
        for c in self.children:
            if c.children:
                yield c
                yield from c.internal()


TOKEN = re.compile(r"\s*([(),:;])\s*|\s*([^(),:;]+)\s*")


def parse(text):
    """Parse one Newick tree. Tip names go to name; internal labels to label."""
    tokens = [(m.group(1), (m.group(2) or "").strip()) for m in TOKEN.finditer(text.strip())]
    pos = 0

    def node():
        nonlocal pos
        n = Node()
        if tokens[pos][0] == "(":
            pos += 1
            n.children.append(node())
            while tokens[pos][0] == ",":
                pos += 1
                n.children.append(node())
            if tokens[pos][0] != ")":
                raise ValueError(f"expected ')' in {text[:80]}")
            pos += 1
            if pos < len(tokens) and tokens[pos][0] is None:
                n.label = tokens[pos][1]
                pos += 1
        else:
            if tokens[pos][0] is not None:
                raise ValueError(f"expected a name in {text[:80]}")
            n.name = tokens[pos][1]
            pos += 1
        if pos < len(tokens) and tokens[pos][0] == ":":
            pos += 1
            n.length = float(tokens[pos][1])
            pos += 1
        return n

    root = node()
    return root


def ingroup_clades(tree, outgroup):
    """(clade, label, length) for every internal branch: clade is the side of
    the branch without the outgroup. Branches whose ingroup side is a single
    taxon or the whole ingroup are not splits of the ingroup and are left out."""
    taxa = frozenset(tree.tips())
    if outgroup not in taxa:
        raise ValueError(f"outgroup {outgroup} not in the tree ({sorted(taxa)})")
    ingroup = taxa - {outgroup}
    out = []
    for n in tree.internal():
        below = frozenset(n.tips())
        side = (taxa - below) if outgroup in below else below
        if 1 < len(side) < len(ingroup):
            out.append((side, n.label, n.length))
    return out


def newick(node):
    """Newick text (no lengths) of a nested tuple tree."""
    return node if isinstance(node, str) else "(" + ",".join(newick(c) for c in node) + ")"


def nested(clades, ingroup):
    """The rooted ingroup tree defined by a set of compatible clades, as
    nested tuples in canonical order (clades before single taxa, then by
    name); unresolved nodes stay multifurcating."""
    ingroup = frozenset(ingroup)
    sets = sorted({frozenset(c) for c in clades} | {ingroup}, key=len)

    def build(members):
        inner = [c for c in sets if c < members]
        maximal = [c for c in inner if not any(c < d for d in inner)]
        covered = frozenset().union(*maximal) if maximal else frozenset()
        parts = [build(c) for c in maximal] + sorted(members - covered)
        parts.sort(key=lambda p: (isinstance(p, str), newick(p)))
        return tuple(parts)

    return build(ingroup)


def rooted_topology(clades, ingroup):
    """Canonical Newick (no lengths) of the rooted ingroup tree."""
    return newick(nested(clades, ingroup)) + ";"


def sister_pair(clades, triple):
    """The two taxa of triple that form a clade excluding the third, or None
    when the tree does not resolve the triple."""
    triple = frozenset(triple)
    for c in clades:
        inside = triple & frozenset(c)
        if len(inside) == 2:
            return inside
    return None


def quartet_split(clades, ingroup):
    """The ingroup split of a four taxon ingroup: the size two clade and its
    complement, as 'a+b | c+d' with names sorted; None if unresolved."""
    ingroup = frozenset(ingroup)
    for c in clades:
        if len(c) == 2:
            first, second = sorted([sorted(c), sorted(ingroup - frozenset(c))])
            return f"{'+'.join(first)} | {'+'.join(second)}"
    return None
