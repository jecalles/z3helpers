import itertools

from synbio.utils import aminoacids, dNTPs, rNTPs
from z3 import EnumSort, Function

__all__ = [
    # str values
    "triplet_dna_codons", "triplet_rna_codons", "aminoacids",
    # z3 Sorts and possible values
    "N", "z3nucleotides", "C", "z3codons", "A", "z3aminos",
    "STOP", "NULL",
    # mappings from str to z3 values
    "dna_to_z3codon", "rna_to_z3codon", "amino_to_z3amino",
    # z3 Functions
    "f_Code", "f_Codon"
]

# useful str collections
triplet_dna_codons = [
    ''.join(codon_list)
    for codon_list in itertools.product(dNTPs, repeat=3)
]
triplet_rna_codons = [
    ''.join(codon_list)
    for codon_list in itertools.product(rNTPs, repeat=3)
]
aminoacids += ["0"]                                     # adds NULL

# useful z3 Sorts
N, z3nucleotides = EnumSort("Nucleotides", dNTPs)       # N = set of nucleotides
C, z3codons = EnumSort("Codons", triplet_dna_codons)    # triplet codons
A, z3aminos = EnumSort("Amino Acids", aminoacids)       # translation signals

# z3amino definitions
STOP, NULL = z3aminos[-2:]

# python dicts to convert between str and z3 sorts
dna_to_z3codon = {
    str_codon: z3_codon
    for str_codon, z3_codon in zip(triplet_dna_codons, z3codons)
}
rna_to_z3codon = {
    str_codon: z3_codon
    for str_codon, z3_codon in zip(triplet_rna_codons, z3codons)
}
amino_to_z3amino = {
    str_amino: z3_amino
    for str_amino, z3_amino in zip(aminoacids, z3aminos)
}

# z3 Functions
f_Codon = Function("nucleotides -> codons", N, N, N, C)     # F(N,N,N) -> C
f_Code = Function("Genetic Code", C, A)                     # F(C) -> A
