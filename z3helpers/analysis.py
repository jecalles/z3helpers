from synbio.codes import Code
from synbio.polymers import DNA

from z3helpers.definitions import dna_to_z3nucleotide
from z3helpers.utils import rmap

__all__ = [
    # extract variables from model
    "code_from_model", "dna_from_model",
    # objective functions
    "objective_function"
]


def code_from_model(code_dict, model):
    def rna(dna): return dna.replace("T", "U")
    def amino(codon): return str(model[code_dict[codon]])

    dict_ = {
        rna(codon): amino(codon) for codon in code_dict.keys()
    }

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

    dict_ = rmap(dna_to_z3nucleotide)

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
