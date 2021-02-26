from typing import List, Sequence

from synbio.utils import get_codons
from synbio.annotations import Location, Part

from z3 import Const

from z3helpers.typing import *
from z3helpers.definitions import *
from z3helpers.utils import decode

__all__ = [
    # z3 variable generating functions
    "dna_variables", "protein_variables"
]


def dna_variables(
        loc: Location,
        nucleotide_sort: NucleotideSort = NucleotideEnumSort
) -> List[NucleotideRef]:
    """
    Creates dna variables for each nucleotide at Location loc:

    :param loc:
    :param nucleotide_sort:
    :return dna_variables:
    """
    return [Const(f"dna_{i}", nucleotide_sort)
            for i in range(loc.start, loc.end)]


def protein_variables(
        T: CodeRef,
        seq_variables: Sequence[NucleotideRef],
        part: Part,
        offset: int = 0,
) -> List[AminoRef]:
    """
    A function that generates z3 variables corresponding to the amino acid
    sequence of the translated Part

   """
    begin = part.location.start - offset
    end = part.location.end - offset

    return [
        decode(T, codon)
        for codon in get_codons(seq_variables[begin:end])
    ]
