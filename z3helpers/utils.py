from typing import Dict, Optional, Sequence
import itertools

from z3helpers.definitions import *
from z3helpers.typing import *

__all__ = [
    "rmap", "z3nuc_to_str", "z3codon_to_str", "z3amino_to_str",
    "decode", "add_constraints"
]


def rmap(dict_):
    return {
        value: key for key, value in dict_.items()
    }


def z3nuc_to_str(
        nucleotide_list: Sequence[NucleotideRef],
        mapping: Optional[Dict[str, NucleotideRef]] = None
) -> str:
    if mapping is None:
        sort = "bv" if isinstance(nucleotide_list[0], BitVecRef) else "enum"
        mapping = dna_to_z3_bv_nuc if sort == "bv" else dna_to_z3_enum_nuc

    rmap_ = rmap(mapping)

    return ''.join(rmap_[nuc] for nuc in nucleotide_list)


def z3codon_to_str(
        codon: CodonRef,
        mapping: Optional[Dict[str, CodonRef]] = None
) -> str:
    if mapping is None:
        mapping = dna_to_z3_enum_codon

    rmap_ = rmap(mapping)

    return rmap_[codon]


def z3amino_to_str(amino: AminoRef,
                   mapping: Optional[Dict[str, CodonRef]] = None):
    if mapping is None:
        sort = "bv" if isinstance(amino, BitVecRef) else "enum"
        mapping = amino_to_z3_bv_amino if sort == "bv" \
            else amino_to_z3_enum_amino

    rmap_ = rmap(mapping)
    return rmap_[amino]


def decode(T: CodeRef, key: CodonRef) -> AminoRef:
    """
    A method used to get amino acids from a CodeRef type object, given a
    CodonRef key obj

    :param T: genetic code
    :param key: codons or nucleotides
    :return amino: amino acid corresponding to "key"

    >>> decode(f_codon_to_amino, z3_enum_codons[3])
    codons -> aminos(TTG)
    >>> decode(
    ...     f_nuc_to_amino,
    ...     (z3_enum_nucleotides[0], z3_enum_nucleotides[1], z3_enum_nucleotides[2])
    ... )
    nucleotides -> aminos(dT, dC, dA)
    >>> decode(Code(), "AUG")
    'M'
    """
    if isinstance(T, dict):
        if isinstance(key, str):
            pass
        elif isinstance(key, DatatypeRef):  # of CodonRef type
            key = z3codon_to_str(key)
        elif isinstance(key[0], DatatypeRef):  # list of NucleotideRef
            key = z3nuc_to_str(key)
        else:
            raise TypeError("key is not of type CodonRef")

        amino = T[key]  # try direct indexing with key

    elif isinstance(T, FuncDeclRef):
        # convert string keys into z3 defs
        if isinstance(key, str):
            key = dna_to_z3_enum_codon[key]
        elif isinstance(key, Sequence) and isinstance(key[0], str):
            sort = "bv" if isinstance(T.domain(0), BitVecSortRef) else "enum"
            mapping = dna_to_z3_bv_nuc if sort == "bv"  else dna_to_z3_enum_nuc
            key = (mapping[n] for n in key)

        if T.arity() == 1:  # accepts codons
            if isinstance(key, Sequence):
                amino = T(f_nuc_to_codon(*key))
            else:
                amino = T(key)

        elif T.arity() == 3:  # accepts nucleotides
            amino = T(*key)
        else:
            raise ValueError(f"T has incorrect arity {T.arity()}")
    else:
        raise TypeError("T is not of type CodeRef")

    return amino


def add_constraints(solver: SolverType,
                    constraints: Iterable[ConstraintRef],
                    weights: Optional[Iterable[int]] = None,
                    hard: bool = True) -> None:

    if hard:
        for constr in constraints:
            solver.add(constr)
    else:
        if weights is None:
            weights = itertools.repeat(1)

        for constr, w in zip(constraints, weights):
            solver.add_soft(constr, weight=w)

    return
