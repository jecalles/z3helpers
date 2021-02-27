import itertools
from typing import List, Sequence

from synbio.utils import aminoacids, dNTPs, rNTPs
from z3 import BitVec, BitVecSort, EnumSort, Function

from z3helpers.typing import *

__all__ = [
    # str values
    "triplet_dna_codons", "triplet_rna_codons", "aminoacids",
    # z3 Sorts and possible values
    "NucleotideEnumSort", "z3nucleotides", "triplet_z3nucleotides",
    "CodonEnumSort", "z3codons",
    "AminoEnumSort", "z3_enum_aminos", "STOP", "NULL",
    "AminoBitVecSort", "z3_bitvec_aminos",
    # mappings from str to z3 values
    "dna_to_z3codon", "dna_to_z3nucleotide", "rna_to_z3codon",
    "amino_to_z3_enum_amino", "amino_to_z3_bitvec_amino",
    # z3 Functions
    "f_nuc_to_codon", "f_codon_to_amino", "f_nuc_to_amino",
    # z3 constant generating functions
    "get_start", "get_stop", "get_null",
]

# useful str collections
triplet_dna_codons: List[CodonRef] = [
    ''.join(codon_list)
    for codon_list in itertools.product(dNTPs, repeat=3)
]
triplet_rna_codons: List[CodonRef] = [
    ''.join(codon_list)
    for codon_list in itertools.product(rNTPs, repeat=3)
]
aminoacids += ["0"]  # adds NULL

# useful z3 Sorts
NucleotideEnumSort, z3nucleotides = EnumSort("Nucleotides",
                                             [f"d{n}" for n in dNTPs])
triplet_z3nucleotides = list(itertools.product(z3nucleotides, repeat=3))
CodonEnumSort, z3codons = EnumSort("Codons", triplet_dna_codons)

AminoEnumSort, z3_enum_aminos = EnumSort("Amino Acids", aminoacids)
AminoBitVecSort = BitVecSort(21)
z3_bitvec_aminos = [
    BitVec((2 ** i) - 1, AminoBitVecSort)
    for i in range(22)
]


# z3amino definitions
def get_start(aminos: Sequence[AminoRef] = z3_enum_aminos) -> AminoRef:
    return aminos[6]


START = z3_enum_aminos[6]
assert str(START) == "M"


def get_stop(aminos: Sequence[AminoRef] = z3_enum_aminos) -> AminoRef:
    return aminos[-2]


def get_null(aminos: Sequence[AminoRef] = z3_enum_aminos) -> AminoRef:
    return aminos[-1]


STOP, NULL = z3_enum_aminos[-2:]
assert str(STOP) == "*"
assert str(NULL) == "0"

# python dicts to convert between str and z3 sorts
# TODO: convert static dicts into functions that generate dicts based on
#  which sorts to use (unary, EnumSort, raw Bools, etc)
dna_to_z3codon = {
    str_codon: z3_codon
    for str_codon, z3_codon in zip(triplet_dna_codons, z3codons)
}
dna_to_z3nucleotide = {
    str_nuc: z3nuc
    for str_nuc, z3nuc in zip(dNTPs, z3nucleotides)
}
rna_to_z3codon = {
    str_codon: z3_codon
    for str_codon, z3_codon in zip(triplet_rna_codons, z3codons)
}
amino_to_z3_enum_amino = {
    str_amino: z3_amino
    for str_amino, z3_amino in zip(aminoacids, z3_enum_aminos)
}
amino_to_z3_bitvec_amino = {
    str_amino: z3_bitvec_amino
    for str_amino, z3_bitvec_amino in zip(aminoacids, z3_bitvec_aminos)
}

# z3 Functions
f_nuc_to_codon = Function(
    "nucleotides -> codons",
    NucleotideEnumSort, NucleotideEnumSort, NucleotideEnumSort,
    CodonEnumSort
)
f_codon_to_amino = Function("codons -> aminos", CodonEnumSort, AminoEnumSort)
f_nuc_to_amino = Function(
    "nucleotides -> aminos",
    NucleotideEnumSort, NucleotideEnumSort, NucleotideEnumSort,
    AminoEnumSort
)


