import itertools
from typing import Dict, List, Optional, Sequence

from synbio.annotations import Location, Part
from synbio.utils import get_codons
from z3 import And, BitVecVal, Extract, Implies, Or, PbEq, PbLe

from z3helpers.definitions import *
from z3helpers.typing import *
from z3helpers.utils import decode


# general constraints
def this_many_true(bool_list: Sequence[BoolRef],
                   num_true: int,
                   weighting=Optional[Sequence[int]]) -> ConstraintRef:
    """
    A function that takes a list of z3 Bool objects (bool_list) and an integer
    representing the number (num_true) of these booleans that should be TRUE.
    Returns a z3 Bool object which is:
        TRUE
            if _exactly_ num_true booleans in bool_list are TRUE, and
        FALSE
            otherwise.

    Optionally, specify the weighting of each boolean in bool_list (default
    is 1, or that every boolean is equally valued)

    :param bool_list: z3 Bool objects to compare
    :param num_true: number of bools that should be true
    :param weighting: optional, how important is each bool (default=1)
    :return: BoolRef
    """
    if weighting is None:
        weighting = [1 for _ in range(len(bool_list))]

    return PbLe(list(zip(bool_list, weighting)), k=num_true)


# constrain nucleotide -> codon mapping
def f_codon_true_mapping(f_codon: FuncDeclRef) -> List[ConstraintRef]:
    """
    Adds constraints on f_nuc_to_codon such that triplets of nucleotides map to
    appropriate codons

    :return: list of constraints
    """
    nucleotide_triplets = list(itertools.product(z3nucleotides, repeat=3))
    return [
        f_codon(n1, n2, n3) == c
        for (n1, n2, n3), c in zip(nucleotide_triplets, z3codons)
    ]


# hard constraints on unary encoding
def amino_bitvec_unary_restriction(
        amino_list: List[AminoRef] = z3_bitvec_aminos
) -> List[ConstraintRef]:
    amino_sort_size = amino_list[0].sort().size() - 1

    return [
        Implies(
            Extract(i + 1, i + 1, amino) == BitVecVal(1, 1),
            Extract(i, i, amino) == BitVecVal(1, 1),
        ) for i in range(amino_sort_size)
        for amino in amino_list
    ]


# hard constraints on Genetic Codes


def at_least_one_codon_per_amino(
        T: CodeRef,
        codons: Sequence[CodonRef] = z3codons,
        aminos: Sequence[AminoRef] = z3_enum_aminos,
        exclude: Optional[AminoRef] = None
) -> List[ConstraintRef]:
    """

    :param T: genetic code
    :param codons: domain of T (list or list-like of codons)
    :param aminos: range of T (list or list-like of aminos)
    :param exclude: range of T to ignore (list or list-like of aminos) (
    default=NULL)
    :return: list of constraints
    """
    if exclude is None:
        exclude = (get_null(aminos),)

    return [
        Or([decode(T, codon) == aa for codon in codons])
        for aa in aminos if aa not in exclude
    ]


def at_most_one_codon_per_amino(
        T: CodeRef,
        codons: Sequence[CodonRef] = z3codons,
        aminos: Sequence[AminoRef] = z3_enum_aminos,
        exclude: Optional[AminoRef] = None
) -> List[ConstraintRef]:
    """
    :param aminos:
    :param T: genetic code
    :param codons: domain of T (list or list-like of codons)
    :param exclude: range of T to ignore (list or list-like of aminos)
    :return: list of constraints
    """
    if exclude is None:
        exclude = (get_null(aminos), get_stop(aminos))

    def _any_codon_in_exclude(T, codons, exclude):
        return Or([
            Or([decode(T, codon) == exc for codon in codons])
            for exc in exclude
        ])

    return [
        Or(
            decode(T, c1) != decode(T, c2),
            # check that codons aren't identical
            _any_codon_in_exclude(T, (c1, c2), exclude)
            # or that codons are in exclude
        ) for c1, c2 in itertools.combinations(codons, r=2)
    ]


def exactly_one_codon_per_amino(
        T: CodeRef,
        codons: Sequence[CodonRef] = z3codons,
        aminos: Sequence[AminoRef] = z3_enum_aminos,
        exclude: Optional[AminoRef] = None
) -> List[ConstraintRef]:
    """

    :param T: genetic code
    :param codons: list of CodonRef objects
    :param aminos: list of AminoRef objects
    :param exclude: range of T to ignore (list or list-like of aminos)
    :return: list of constraints
    """
    if exclude is None:
        exclude = (get_null(aminos), get_stop(aminos))

    return [
        PbEq([(decode(T, c) == aa, 1) for c in codons], k=1)
        for aa in aminos if aa not in exclude
    ]


def compatible_with_standard_code(
        T: CodeRef,
        dna_codons: Sequence[CodonRef] = triplet_dna_codons,
        aminos: Sequence[AminoRef] = z3_enum_aminos,
        amino_dict: Dict[str, AminoRef] = amino_to_z3_enum_amino
) -> List[ConstraintRef]:
    null = get_null(aminos)

    sc_constraints = standard_code(T, dna_codons, amino_dict)

    return [
        Or(sc_const, decode(T, codon) == null)
        for codon, sc_const in zip(dna_codons, sc_constraints)
    ]


def standard_code(
        T: CodeRef,
        dna_codons: Sequence[CodonRef] = triplet_dna_codons,
        # codon_dict: Dict[str, NucleotideRef] = dna_to_z3codon,
        amino_dict: Dict[str, AminoRef] = amino_to_z3_enum_amino
) -> List[ConstraintRef]:
    dna_to_rna = dict(zip(dna_codons, triplet_rna_codons))
    sc = {
        codon: amino_dict[Code()[dna_to_rna[codon]]]
        for codon in dna_codons
    }

    return [
        decode(T, codon) == sc[codon]
        for codon in dna_codons
    ]


# define sequence based constraints
def translation_constraints(
        T: CodeRef,
        dna_variables: Sequence[NucleotideRef],
        prot_variables: Sequence[AminoRef],
        location: Location,
        offset: int = 0,
        nucleotides: Sequence[NucleotideRef] = z3nucleotides,
        aminos: Sequence[AminoRef] = z3_enum_aminos,
        start_flag: bool = False,
        stop_flag: bool = False
) -> List[ConstraintRef]:
    """
    A function that generates constraints on DNA variables that are involved
    in encoding proteins (must start with "M", no NULL codons in middle,
    must end with STOP)

    :param T: genetic code
    :param dna_variables:
    :param prot_variables:
    :param location:
    :param offset:
    :param nucleotides:
    :param aminos:
    :param start_flag:
    :param stop_flag:
    :return:
    """
    # translation signals
    START = get_start(aminos)
    STOP = get_stop(aminos)
    NULL = get_null(aminos)

    # locations
    begin = location.start - offset
    end = location.end - offset

    codon_variable_list = get_codons(dna_variables[begin:end])

    # ensure proteins start with Met
    start_constraints = []
    if start_flag:
        start_constraints += [prot_variables[0] == START]

    # ensure no NULL in middle of protein and proper termination
    null_constraints = [
        aa != NULL for aa in prot_variables
    ]

    stop_constraints = [
        aa != STOP for aa in prot_variables[:-1]
    ]

    if stop_flag:
        stop_constraints += [prot_variables[-1] == STOP]
    else:
        stop_constraints += [prot_variables[-1] != STOP]

    code_implications = [
        # TODO: try refactoring to create intermediate variables
        # TODO: refactor to generalize over Sorts
        Implies(And(codon_variable[0] == z3codon[0],
                    codon_variable[1] == z3codon[1],
                    codon_variable[2] == z3codon[2]),
                z3amino == decode(T, z3codon))
        for codon_variable, z3amino in zip(codon_variable_list, prot_variables)
        for z3codon in itertools.product(nucleotides, repeat=3)
    ]

    return start_constraints + stop_constraints + null_constraints + \
           code_implications


def same_sequence(
        dna_variables: Sequence[NucleotideRef],
        wt_sequence: str
) -> List[ConstraintRef]:
    seq_in_z3nucleotides = (
        dna_to_z3nucleotide[n] for n in wt_sequence
    )
    return [
        variable == wt_value
        for variable, wt_value in zip(dna_variables, seq_in_z3nucleotides)
    ]


def translates_same(
        prot_variables: Sequence[AminoRef],
        part: Part,
        amino_dict: Dict[str, AminoRef] = amino_to_z3_enum_amino
) -> List[ConstraintRef]:
    prot_seq = str(part.seq.translate())

    return [
        var_aa == amino_dict[wt_aa]
        for var_aa, wt_aa in zip(prot_variables, prot_seq)
    ]


# bundled constraints
def FS20(
        T: CodeRef,
        f_codon: FuncDeclRef = f_nuc_to_codon,
        codons: Sequence[CodonRef] = z3codons,
        aminos: Sequence[AminoRef] = z3_enum_aminos,
) -> List[ConstraintRef]:
    constraints = exactly_one_codon_per_amino(T, codons, aminos)
    # if T is a Function, add true nucleotide -> codon mapping
    if isinstance(T, FuncDeclRef) and T.arity() == 1:
        constraints += f_codon_true_mapping(f_codon)

    return constraints


def RED20(
        T: CodeRef,
        f_codon: FuncDeclRef = f_nuc_to_codon,
        codons: Sequence[CodonRef] = z3codons,
        aminos: Sequence[AminoRef] = z3_enum_aminos,
        amino_dict: Dict[str, AminoRef] = amino_to_z3_enum_amino
) -> List[ConstraintRef]:
    return FS20(T, f_codon, codons, aminos) + compatible_with_standard_code(
        T, codons, aminos, amino_dict)
