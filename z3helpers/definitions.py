import itertools

from synbio.annotations import Location
from synbio.utils import aminoacids, dNTPs, rNTPs
from z3 import Const, EnumSort, Function

__all__ = [
    # str values
    "triplet_dna_codons", "triplet_rna_codons", "aminoacids",
    # z3 Sorts and possible values
    "NucleotideSort", "z3nucleotides", "CodonSort", "z3codons",
    "AminoSort", "z3aminos", "STOP", "NULL",
    # mappings from str to z3 values
    "dna_to_z3codon", "rna_to_z3codon", "amino_to_z3amino",
    # z3 Functions
    "f_nuc_to_codon", "f_codon_to_amino", "f_nuc_to_amino",
    # z3 Constant generating functions
    "sequence_variables"
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
aminoacids += ["0"]  # adds NULL

# useful z3 Sorts
NucleotideSort, z3nucleotides = EnumSort("Nucleotides", dNTPs)
CodonSort, z3codons = EnumSort("Codons", triplet_dna_codons)
AminoSort, z3aminos = EnumSort("Amino Acids", aminoacids)

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
f_nuc_to_codon = Function("nucleotides -> codons",
                          NucleotideSort, NucleotideSort, NucleotideSort,
                          CodonSort)
f_codon_to_amino = Function("codons -> aminos", CodonSort, AminoSort)
f_nuc_to_amino = Function("nucleotides -> aminos",
                          NucleotideSort, NucleotideSort, NucleotideSort,
                          AminoSort)


# z3 Constant generating functions
def sequence_variables(loc: Location):
    return [Const(f"x{i}", NucleotideSort)
            for i in range(loc.start, loc.end)]
