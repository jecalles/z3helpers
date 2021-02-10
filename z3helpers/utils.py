from typing import Dict, Iterable, Optional, Sequence

from z3helpers.definitions import *
from z3helpers.typing import *


def rmap(dict_):
    return {
        value: key for key, value in dict_.items()
    }


def z3nuc_to_str(
        nucleotide_list: Sequence[NucleotideRef],
        mapping: Dict[str, NucleotideRef] = dna_to_z3nucleotide
) -> str:
    rmap_ = rmap(mapping)

    return ''.join(rmap_[nuc] for nuc in nucleotide_list)


def z3codon_to_str(codon: CodonRef,
                   mapping: Dict[str, CodonRef] = dna_to_z3codon) -> str:
    rmap_ = rmap(mapping)
    return rmap_[codon]


def z3amino_to_str(amino: AminoRef,
                   mapping: Optional[Dict[str, CodonRef]] = None):
    if mapping is None:
        if isinstance(amino, DatatypeRef):
            mapping = amino_to_z3_enum_amino
        elif isinstance(amino, BitVecRef):
            mapping = amino_to_z3_bitvec_amino

    rmap_ = rmap(mapping)
    return rmap_[amino]


def decode(T: CodeRef, key: CodonRef) -> AminoRef:
    """
    A method used to get amino acids from a CodeRef type object, given a
    CodonRef key obj

    :param T: genetic code
    :param key: codons or nucleotides
    :return amino: amino acid corresponding to "key"

    >>> decode(f_codon_to_amino, z3codons[3])
    codons -> aminos(TTG)
    >>> decode(f_nuc_to_amino, (z3nucleotides[0], z3nucleotides[1], z3nucleotides[2]))
    nucleotides -> aminos(dT, dC, dA)
    >>> decode(Code(), "AUG")
    'M'
    """
    if isinstance(T, dict):
        if isinstance(key, DatatypeRef):  # of CodonRef type
            key = z3codon_to_str(key)

        elif isinstance(key[0], DatatypeRef):  # list of NucleotideRef
            key = z3nuc_to_str(key)

        amino = T[key]  # try direct indexing with key

    elif isinstance(T, FuncDeclRef):
        if T.arity() == 1:  # accepts codons
            amino = T(key)
        elif T.arity() == 3:  # accepts nucleotides
            n1, n2, n3 = key
            amino = T(n1, n2, n3)
        else:
            raise ValueError(f"T has incorrect arity {T.arity()}")
    else:
        raise TypeError("T is not of type CodeRef")

    return amino


def add_constraints(solver: SolverType,
                    constraints: Iterable,
                    hard: bool = True) -> None:
    method = solver.add if hard else solver.add_soft

    for constr in constraints:
        method(constr)

    return