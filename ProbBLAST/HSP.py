from CONSTANTS import *


def get_new_score(
    l_extension,
    r_extension,
    prob_nucleotide,
    query,
    start_real,
    start_created,
    l_index,
    l_score,
    r_index,
    r_score,
):
    end_real = start_real + WORD_SIZE - 1
    end_created = start_created + WORD_SIZE - 1

    if not l_extension:
        if r_extension:
            for nuc in prob_nucleotide[end_real + r_index]:
                if nuc == query[end_created + r_index]:
                    r_score += prob_nucleotide[end_real + r_index][nuc]
                else:
                    r_score += -prob_nucleotide[end_real + r_index][nuc] / 3
            return r_score
        return

    for nuc in prob_nucleotide[start_real - l_index]:
        if nuc != query[start_created - l_index]:
            l_score += -prob_nucleotide[start_real - l_index][nuc] / 3
        else:
            l_score += prob_nucleotide[start_real - l_index][nuc]
    return l_score


def check_DELTA(score, score_max, index, extension, highest_index):
    if score + DELTA < score_max:
        extension = False
    if score >= score_max:
        score_max = score
        highest_index = index
    return extension, score_max, highest_index


def find_HSP(start_real, start_created, prob_nucleotide, query, database):

    end_real = start_real + WORD_SIZE - 1
    end_created = start_created + WORD_SIZE - 1

    l_extension, r_extension = True, True

    l_score_max = 0
    l_score = 0
    l_index_max = 0
    l_index = 0

    r_score_max = 0
    r_score = 0
    r_index_max = 0
    r_index = 0

    while l_extension or r_extension:
        if l_extension:
            l_index += 1
            if (start_real - l_index < 0) or (start_created - l_index < 0):
                l_extension = False
            l_score = get_new_score(
                l_extension,
                r_extension,
                prob_nucleotide,
                query,
                start_real,
                start_created,
                l_index,
                l_score,
                end_real,
                end_created,
                r_index,
                r_score,
            )
            l_extension, l_score_max, l_index_max = check_DELTA(
                l_score, l_score_max, l_index, l_extension, l_index_max
            )

        if r_extension:
            r_index += 1
            if (end_real + r_index >= len(database)) or (
                end_created + r_index >= QUERY_LENGTH
            ):
                r_extension = False
            r_score = get_new_score(
                l_extension,
                r_extension,
                prob_nucleotide,
                query,
                start_real,
                start_created,
                l_index,
                l_score,
                end_real,
                end_created,
                r_index,
                r_score,
            )
            r_score_max, r_extension, r_index_max = check_DELTA(
                r_score, r_score_max, r_index, r_extension, r_index_max
            )

    score = l_score_max + r_score_max
    return score, l_index_max, r_index_max


def find_words_HSP(query, index_words, prob_nucleotide, database):
    HSP = dict()
    for i_start_real in range(QUERY_LENGTH - WORD_SIZE + 1):
        for i_start_created in index_words[
            query[i_start_real : i_start_real + WORD_SIZE]
        ]:
            score, l_index, r_index = find_HSP(
                i_start_created, i_start_real, prob_nucleotide, query, database
            )

            if score > HSP_THREASHOLD:
                HSP[
                    (
                        i_start_created + WORD_SIZE + r_index - 1,
                        i_start_created - l_index,
                        i_start_real - l_index,
                        i_start_real + WORD_SIZE + r_index - 1,
                    )
                ] = score
    return HSP
