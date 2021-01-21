import itertools

from synbio.utils import aminoacids, dNTPs, rNTPs
from z3 import EnumSort

# useful str collections
triplet_dna_codons = [
    ''.join(codon_list)
    for codon_list in itertools.product(dNTPs, repeat=3)
]
triplet_rna_codons = [
    ''.join(codon_list)
    for codon_list in itertools.product(rNTPs, repeat=3)
]
aminoacids += ["0"]

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
