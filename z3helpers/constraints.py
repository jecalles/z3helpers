import itertools

from synbio.codes import Code
from synbio.utils import get_codons
from z3 import Or, PbEq, PbLe

from z3helpers.definitions import *


# general constraints
def this_many_booleans(bool_list, threshold, weighting=None):
    if weighting is None:
        weighting = [1 for _ in range(len(bool_list))]

    return PbLe(list(zip(bool_list, weighting)), k=threshold)


# constrain nucleotide -> codon mapping
def f_codon_true_mapping(f_codon):
    """
    Adds constraints on f_nuc_to_codon such that triplets of nucleotides map to
    appropriate codons

    :return:
    List[Bool]: constraints to be satisfied
    """
    nucleotide_triplets = list(itertools.product(z3nucleotides, repeat=3))
    return [
        f_codon(n1, n2, n3) == c
        for (n1, n2, n3), c in zip(nucleotide_triplets, z3codons)
    ]


# hard constraints on Genetic Codes


def at_least_one_codon_per_amino(T, exclude=(NULL,)):
    return [
        Or([T(c) == aa for c in z3codons])
        for aa in z3aminos if aa not in exclude
    ]


def at_most_one_codon_per_amino(T, exclude=(NULL, STOP)):
    def _any_codon_in_exclude(codons, T, exclude):
        return Or([
            Or([T(c) == exc for c in codons])
            for exc in exclude
        ])

    return [
        Or(
            T(c1) != T(c2),  # check that codons aren't identical
            _any_codon_in_exclude([c1, c2], T, exclude)
            # or that codons are in exclude
        ) for c1, c2 in itertools.combinations(z3codons, r=2)
    ]


def exactly_one_codon_per_amino(T, exclude=(NULL, STOP)):
    return [
        PbEq([(T(c) == aa, 1) for c in z3codons], k=1)
        for aa in z3aminos if aa not in exclude
    ]


def compatible_with_standard_code(T):
    sc = Code()
    return [
        Or(T(rna_to_z3codon[rna_codon]) == amino_to_z3amino[sc[rna_codon]],
           T(rna_to_z3codon[rna_codon]) == NULL)
        for rna_codon in triplet_rna_codons
    ]


# define sequence based constraints
def translates_same(T, f_codon, seq_variables, part, offset):
    start = part.location.start - offset
    end = part.location.end - offset
    codon_list = get_codons(seq_variables[start:end])

    protein_seq = str(part.seq.translate())

    return [
        T(f_codon(n1, n2, n3)) == amino_to_z3amino[aa]
        for (n1, n2, n3), aa in zip(codon_list, protein_seq)
    ]


# bundled constraints
def FS20(T=f_codon_to_amino, f_codon=f_nuc_to_codon):
    return f_codon_true_mapping(f_codon) + exactly_one_codon_per_amino(T)


def RED20(T=f_codon_to_amino, f_codon=f_nuc_to_codon):
    return FS20(T, f_codon) + compatible_with_standard_code(T)
