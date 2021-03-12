import itertools
from typing import List, Sequence

from synbio.utils import aminoacids, dNTPs, rNTPs
from z3 import BitVec, BitVecSort, EnumSort, Function

from z3helpers.typing import *

__all__ = [
    # str values
    "triplet_dna_codons", "triplet_rna_codons", "aminoacids",
    # z3 Sorts and possible values
    "NucleotideEnumSort", "z3_enum_nucleotides", "triplet_z3_enum_nucleotides",
    "NucleotideBitVecSort", "z3_bitvec_nucleotides",
    "triplet_z3_bitvec_nucleotides",
    "CodonEnumSort", "z3_enum_codons",
    "AminoEnumSort", "z3_enum_aminos",
    "AminoBitVecSort", "z3_bitvec_aminos",
    # mappings from str to z3 values
    "dna_to_z3_bv_nuc", "dna_to_z3_enum_nuc",
    "dna_to_z3_enum_codon", "rna_to_z3_enum_codon",
    "amino_to_z3_bv_amino", "amino_to_z3_enum_amino",
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
NucleotideEnumSort, z3_enum_nucleotides = EnumSort(
    "Nucleotides", [f"d{n}" for n in dNTPs]
)
triplet_z3_enum_nucleotides = list(
    itertools.product(z3_enum_nucleotides, repeat=3)
)
NucleotideBitVecSort = BitVecSort(3)
z3_bitvec_nucleotides = [
    BitVec((2 ** i) - 1, NucleotideBitVecSort)
    for i in range(4)
]
triplet_z3_bitvec_nucleotides = list(
    itertools.product(z3_bitvec_nucleotides, repeat=3)
)

CodonEnumSort, z3_enum_codons = EnumSort("Codons", triplet_dna_codons)

AminoEnumSort, z3_enum_aminos = EnumSort("Amino Acids", aminoacids)
AminoBitVecSort = BitVecSort(21)
z3_bitvec_aminos = [
    BitVec((2 ** i) - 1, AminoBitVecSort)
    for i in range(22)
]


# z3amino definitions
def get_start(aminos: Sequence[AminoRef] = z3_enum_aminos) -> AminoRef:
    return aminos[6]


def get_stop(aminos: Sequence[AminoRef] = z3_enum_aminos) -> AminoRef:
    return aminos[-2]


def get_null(aminos: Sequence[AminoRef] = z3_enum_aminos) -> AminoRef:
    return aminos[-1]


# python dicts to convert between str and z3 sorts
# TODO: convert static dicts into functions that generate dicts based on
#  which sorts to use (unary, EnumSort, raw Bools, etc)
dna_to_z3_bv_nuc = {
    str_nuc: z3nuc
    for str_nuc, z3nuc in zip(dNTPs, z3_bitvec_nucleotides)
}
dna_to_z3_enum_codon = {
    str_codon: z3_codon
    for str_codon, z3_codon in zip(triplet_dna_codons, z3_enum_codons)
}
dna_to_z3_enum_nuc = {
    str_nuc: z3nuc
    for str_nuc, z3nuc in zip(dNTPs, z3_enum_nucleotides)
}
rna_to_z3_enum_codon = {
    str_codon: z3_codon
    for str_codon, z3_codon in zip(triplet_rna_codons, z3_enum_codons)
}
amino_to_z3_enum_amino = {
    str_amino: z3_amino
    for str_amino, z3_amino in zip(aminoacids, z3_enum_aminos)
}
amino_to_z3_bv_amino = {
    str_amino: z3_bitvec_amino
    for str_amino, z3_bitvec_amino in zip(aminoacids, z3_bitvec_aminos)
}

# z3 Functions
f_nuc_to_codon = Function(
    "nucleotides -> codons",
    NucleotideBitVecSort, NucleotideBitVecSort, NucleotideBitVecSort,
    CodonEnumSort
)
f_codon_to_amino = Function(
    "codons -> aminos",
    CodonEnumSort,
    AminoBitVecSort
)
f_nuc_to_amino = Function(
    "nucleotides -> aminos",
    NucleotideBitVecSort, NucleotideBitVecSort, NucleotideBitVecSort,
    AminoBitVecSort
)
