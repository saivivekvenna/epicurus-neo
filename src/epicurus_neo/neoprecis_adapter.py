from __future__ import annotations

from collections import OrderedDict


def normalize_neoprecis_allele(value: str) -> str:
    allele = str(value).removeprefix("HLA-")
    if "*" not in allele:
        allele = allele[0] + "*" + allele[1:]
    return allele


def map_sequence_to_core(sequence: str, core: str) -> OrderedDict[int, int]:
    mapping: OrderedDict[int, int] = OrderedDict()
    sequence_index = 0
    core_index = 0
    while sequence_index < len(sequence) and core_index < len(core):
        if sequence[sequence_index] == core[core_index]:
            mapping[sequence_index] = core_index
            sequence_index += 1
            core_index += 1
        elif len(sequence) > len(core):
            sequence_index += 1
        elif len(sequence) < len(core):
            core_index += 1
        else:
            break
    return mapping


def wildtype_pseudo_core(wildtype: str, mutant: str, mutant_core: str) -> str:
    if len(wildtype) != len(mutant):
        return ""
    pseudo_core = list(mutant_core)
    for sequence_index, core_index in map_sequence_to_core(mutant, mutant_core).items():
        pseudo_core[core_index] = wildtype[sequence_index]
    return "".join(pseudo_core)
