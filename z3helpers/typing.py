from typing import Union, Dict, Sequence
from typing import TypeVar

from z3 import (
    FuncDeclRef, BitVecRef, DatatypeRef, BoolRef,
    DatatypeSortRef, BitVecSortRef,
    Solver, Optimize,
)
from synbio.codes import Code
from synbio.polymers import DNA, Protein

__all__ = [
    # z3 defined classes
    "FuncDeclRef", "BitVecRef", "DatatypeRef", "BoolRef",
    "DatatypeSortRef", "BitVecSortRef",
    # user defined types
    "NucleotideRef", "AminoRef", "CodonRef", "CodeRef", "ConstraintRef",
    "NucleotideSort", "AminoSort", "CodonSort",
    "SolverType"
]
NucleotideRef = TypeVar("NucleotideRef", str, DNA, DatatypeRef, BitVecRef)
NucleotideSort = TypeVar("NucleotideSort",
                         str, DNA, DatatypeSortRef, BitVecSortRef)

AminoRef = TypeVar("AminoRef", str, Protein, DatatypeRef, BitVecRef)
AminoSort = TypeVar("AminoSort", str, Protein, DatatypeSortRef, BitVecSortRef)

CodonRef = TypeVar("CodonRef", NucleotideRef, Sequence[NucleotideRef])
CodonSort = TypeVar("CodonSort", NucleotideSort, Sequence[NucleotideSort])

# codes are either dict[codon --> amino] or z3.Function
CodeRef = TypeVar("CodeRef", Dict[CodonRef, AminoRef], Code, FuncDeclRef)

# constraints are either bool or z3.BoolRef
ConstraintRef = TypeVar("ConstraintRef", bool, BoolRef)

# SolverType is either z3.Solver or z3.Optimize
SolverType = TypeVar("SolverType", Solver, Optimize)
