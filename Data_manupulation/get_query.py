import numpy as np
from CONSTANTS import *


def delete_nuc(query):
    for i in range(QUERY_LENGTH):
        if np.random.uniform(0.0, 1.0) < INDEL_FRAQ / 2:
            if i == 0:
                return query[1:]
            elif i == QUERY_LENGTH - 1:
                return query[:-1]
        return query[:i] + query[i + 1:]


def insert_nuc(query):
    nuc = np.random.choice(NUC_LIST)
    for i in range(QUERY_LENGTH + 1):
        if np.random.uniform(0.0, 1.0) < INDEL_FRAQ / 2:
            if i == 0:
                return nuc + query
            elif i == QUERY_LENGTH:
                return query + nuc
        return query[:i] + nuc + query[i:]


def choose_query(database, prob_nucleotide):
    start_pos = np.random.randint(0, len(database) - QUERY_LENGTH)
    query = ""
    for i in range(QUERY_LENGTH):
        query += np.random.choice(
            list(prob_nucleotide[start_pos + i].keys()),
            p=list(prob_nucleotide[start_pos + i].values()),
        )

    query = delete_nuc(query)
    query = insert_nuc(query)
    return start_pos, query


def create_new_database(database, prob_nucleotide):
    new_database = ""
    for i in range(len(database)):
        new_database += np.random.choice(
            list(prob_nucleotide[i].keys()), p=list(prob_nucleotide[i].values())
        )
    index_words = {}
    for i in range(len(new_database) - WORD_SIZE + 1):
        if new_database[i: i + WORD_SIZE] in index_words.keys():
            index_words[new_database[i: i + WORD_SIZE]].append(i)
        else:
            index_words[new_database[i: i + WORD_SIZE]] = [i]
    return index_words, new_database
