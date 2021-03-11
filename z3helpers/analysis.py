from synbio.codes import Code
from synbio.polymers import DNA

from z3helpers.definitions import (
    dna_to_z3_bv_nuc, dna_to_z3_enum_nuc, triplet_dna_codons
)
from z3helpers.utils import rmap, decode
from z3helpers.typing import *

__all__ = [
    # extract variables from model
    "code_from_model", "dna_from_model",
    # objective functions
    "objective_function"
]


def code_from_model(code: CodeRef, model) -> Code:
    def rna(dna): return dna.replace("T", "U")

    if isinstance(code, dict):
        def amino(codon): return str(model[decode(code, codon)])

        dict_ = {
            rna(codon): amino(codon) for codon in code.keys()
        }
    elif isinstance(code, FuncDeclRef):
        code_as_list = model[code].as_list()
        null = str(code_as_list[-1])

        if code.arity() == 3:
            codons_and_aminos_chunked = (
                (tuple(args[:3]), args[-1])
                for args in code_as_list[:-1]
            )
            codon_tuples, aminos = zip(*codons_and_aminos_chunked)

            def codon_from_tuple(tup):
                sort = "bv" if isinstance(tup[0], BitVecRef) else "enum"
                mapping = dna_to_z3_bv_nuc if sort == "bv" \
                    else dna_to_z3_enum_nuc

                return ''.join(
                    rmap(mapping)[n]
                    for n in tup
                )
            codons = map(codon_from_tuple, codon_tuples)
        else:
            codons, aminos = zip(*code_as_list[:-1])

        dict_ = {rna(codon): null for codon in triplet_dna_codons}
        update = {
            rna(str(codon)): str(amino)
            for codon, amino in zip(codons, aminos)
        }
        dict_.update(update)

    else:
        raise ValueError("code is not of type CodeRef")

    return Code(dict_)


def dna_from_model(dna_variables, model):
    class AmbiguousDNA(DNA):
        basepairing = {
            k: v for k, v in DNA.basepairing.items()
        }
        basepairing['X'] = 'X'

        def alphabet(self):
            dna_alphabet = super().alphabet()
            dna_alphabet.append('X')
            return dna_alphabet

    sort = "bv" if isinstance(dna_variables[0], BitVecRef) else "enum"
    mapping = dna_to_z3_bv_nuc if sort == "bv" else dna_to_z3_enum_nuc
    dict_ = rmap(mapping)

    dna_assignments = (
        model[var] for var in dna_variables
    )
    dna_str = ''.join(
        'X' if nuc is None else dict_[nuc]
        for nuc in dna_assignments
    )
    if any(nt == 'X' for nt in dna_str):
        out_seq = AmbiguousDNA(dna_str)
    else:
        out_seq = DNA(dna_str)

    return out_seq


def objective_function(out_seq, wt_seq, weights=None):
    if weights is None:
        weights = [1] * len(out_seq)

    matches = (
        str(a1) == str(a2)
        for a1, a2 in zip(out_seq, wt_seq)
    )

    return sum(
        match * weight
        for match, weight in zip(matches, weights)
    )
