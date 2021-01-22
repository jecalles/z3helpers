import itertools

from z3 import Or

from z3helpers.definitions import *


# constrain nucleotide -> codon mapping
def f_Codons_true_mapping():
    """
    Adds constraints on f_Codon such that triplets of nucleotides map to
    appropriate codons

    :return:
    List[Bool]: constraints to be satisfied
    """
    triples = list(itertools.product(z3nucleotides, repeat=3))
    return [
        f_Codon(n1, n2, n3) == z3codon
        for (n1, n2, n3), z3codon in zip(triples, z3codons)
    ]


# hard constraints on Genetic Codes


def at_least_one_codon_per_amino(T, exclude=(NULL)):
    return [
        Or([T(z3codon) == z3amino for z3codon in z3codons])
        for z3amino in z3aminos if z3amino not in exclude
    ]


def at_most_one_codon_per_amino(T, exclude=(NULL, STOP)):
    def _any_codon_in_exclude(z3codons, T, exclude):
        return Or([
            Or([T(z3codon) == exc for z3codon in z3codons])
            for exc in exclude
        ])

    return [
        Or(
            T(c1) != T(c2),  # check that codons aren't identical
            _any_codon_in_exclude([c1, c2], T, exclude)
            # or that codons are in exclude
        ) for c1, c2 in itertools.combinations(z3codons, r=2)
    ]
