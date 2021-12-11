from typing import Any, Union
from CONSTANTS import *


def read_seq(fasta_file, probabilities_file):
    with open(fasta_file, "r") as fasta, open(probabilities_file, "r") as prob:
        database = fasta.readlines()[0]
        prob_file = prob.readlines()[0]

    probability = [float(i) for i in prob_file.split()]
    print("database", len(database))
    return database, probability


def get_proba_nuc(database, probability):
    prob_nucleotide: list[dict[str, Union[float, Any]]] = list()

    for i in range(len(database)):
        nuc_dict = {}
        for nuc in NUC_LIST:
            if nuc == database[i]:
                nuc_dict[nuc] = probability[i]
            else:
                nuc_dict[nuc] = (1 - probability[i]) / 3
        prob_nucleotide.append(nuc_dict)
    print(prob_nucleotide[:10])
    return prob_nucleotide
