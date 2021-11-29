from typing import List, Sequence

from synbio.utils import get_codons
from synbio.annotations import Location, Part

from z3 import Const

from z3helpers.typing import *
from z3helpers.definitions import *
from z3helpers.utils import decode

__all__ = [
    # z3 variable generating functions
    "dna_variables", "protein_variables", "code_dict"
]


def dna_variables(
        loc: Location,
        nucleotide_sort: NucleotideSort = NucleotideBitVecSort
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
        amino_sort: AminoSort = AminoBitVecSort,
) -> List[AminoRef]:
    """
    A function that generates z3 variables corresponding to the amino acid
    sequence of the translated Part

    """
    begin = part.location.start - offset
    end = part.location.end - offset

    codons = get_codons(seq_variables[begin:end])
    if isinstance(T, dict):
        prot_seq = [
            Const(f"{part.name}_{i}", amino_sort)
            for i, _ in enumerate(codons)
        ]

    elif isinstance(T, FuncDeclRef):
        prot_seq = [
            decode(T, codon)
            for codon in codons
        ]

    else:
        raise TypeError(f"T is not of type CodeRef")

    return prot_seq


def code_dict(codons: Sequence[CodonRef] = triplet_dna_codons,
              amino_sort: AminoSort = AminoBitVecSort) -> CodeRef:
    """
    A function that returns a python Dict[str -> AminoRef] mapping DNA codons
    to amino acids

    :param codons:
    :param amino_sort:
    :return code:
    """
    return {
        codon: Const(f"T({codon})", amino_sort)
        for codon in codons
    }
