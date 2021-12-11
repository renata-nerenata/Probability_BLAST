from Data_manupulation.get_init_data import read_seq, get_proba_nuc
from Data_manupulation.get_query import create_new_database, choose_query

from ProbBLAST.HSP import find_words_HSP
from ProbBLAST.extension_NW import Needleman_Wunsch


import time


def main():
    start = time.process_time()
    database, probability = read_seq(
        "Data/chr22.maf.ancestors.42000000.complete.boreo.conf",
        "Data/chr22.maf.ancestors.42000000.complete.boreo.fa",
    )
    prob_nucleotide = get_proba_nuc(database, probability)
    start_pos, query = choose_query(database, prob_nucleotide)
    index_words, new_database = create_new_database(database, prob_nucleotide)
    HSP_dict = find_words_HSP(query, index_words, prob_nucleotide, database)
    Needleman_Wunsch(HSP_dict, query, prob_nucleotide, new_database)
    print(time.process_time() - start)


if __name__ == "__main__":
    main()
