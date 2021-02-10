from typing import Union, Dict, Sequence
from typing import NewType

from z3 import (
    FuncDeclRef, BitVecRef, DatatypeRef, BoolRef,
    DatatypeSortRef, BitVecSortRef,
    Solver, Optimize
)
from synbio.codes import Code
from synbio.polymers import DNA, Protein

__all__ = [
    # z3 defined classes
    "FuncDeclRef", "BitVecRef", "DatatypeRef", "BoolRef",
    "DatatypeSortRef", "BitVecSortRef",
    # synbio defined classes
    "DNA", "Protein", "Code",
    # user defined types
    "NucleotideRef", "AminoRef", "CodonRef", "CodeRef", "ConstraintRef",
    "NucleotideSort", "AminoSort", "CodonSort",
    "SolverType"
]
NucleotideRef = Union[str, DNA, DatatypeRef, BitVecRef]
NucleotideSort = Union[str, DNA, DatatypeSortRef, BitVecSortRef]

AminoRef = Union[str, Protein, DatatypeRef, BitVecRef]
AminoSort = Union[str, Protein, DatatypeSortRef, BitVecSortRef]

CodonRef = Union[NucleotideRef, Sequence[NucleotideRef]]
CodonSort = Union[NucleotideSort, Sequence[NucleotideSort]]

# codes are either dict[codon --> amino] or z3.Function
CodeRef = Union[Dict[CodonRef, AminoRef], Code, FuncDeclRef]

# constraints are either bool or z3.BoolRef
ConstraintRef = Union[bool, BoolRef]

# SolverType is either z3.Solver or z3.Optimize
SolverType = Union[Solver, Optimize]
