import numpy as np
from CONSTANTS import *


def forward_trace(score_matrix, path_matrix, query, prob_nucleotide, curr_HSP, left):
    if left:
        index = 2
    else:
        index = 3

    for i in range(len(score_matrix)):
        score_matrix[i][0] = -1 * i
        path_matrix[i][0] = "UP"
    for i in range(len(score_matrix[0])):
        score_matrix[0][i] = -1 * i
        path_matrix[0][i] = "LEFT"
    for i in range(1, len(score_matrix)):
        for j in range(1, len(score_matrix[0])):
            a1 = score_matrix[i - 1][j - 1]
            for nuc in NUC_LIST:
                if nuc == query[curr_HSP[index] - j]:
                    a1 += prob_nucleotide[curr_HSP[index - 2] - i][nuc]
                else:
                    a1 -= prob_nucleotide[curr_HSP[index - 2] - i][nuc]
            a2 = score_matrix[i - 1][j] - 1
            a3 = score_matrix[i][j - 1] - 1
            score_matrix[i][j] = max(a1, a2, a3)
            if a1 == score_matrix[i][j]:
                path_matrix[i][j] = "DIAG"
            elif a2 == score_matrix[i][j]:
                path_matrix[i][j] = "UP"
            else:
                path_matrix[i][j] = "LEFT"
    return score_matrix, path_matrix


def back_trace(score_matrix, path_matrix, query, new_database, curr_HSP, left):
    max_score = max([i[-1] for i in score_matrix])
    j = len(score_matrix[0]) - 1
    i = [i[-1] for i in score_matrix].index(max_score)
    while (i >= 0 and j >= 0) and not (i == 0 and j == 0):
        if left:
            index = 2
            if path_matrix[i][j] == "DIAG":
                seq_real = seq_real + query[curr_HSP[index] - j]
                seq_created = seq_created + new_database[curr_HSP[index - 2] - i]
                i -= 1
                j -= 1
            elif path_matrix[i][j] == "UP":
                seq_real = seq_real + "_"
                seq_created = seq_created + new_database[curr_HSP[index - 2] - i]
                i -= 1
            elif path_matrix[i][j] == "LEFT":
                seq_real = seq_real + query[curr_HSP[index] - j]
                seq_created = seq_created + "_"
                j -= 1
        else:
            index = 3
            if path_matrix[i][j] == "DIAG":
                seq_real = query[curr_HSP[index] + j] + seq_real
                seq_created = new_database[curr_HSP[index - 2] + i] + seq_created
                i -= 1
                j -= 1
            elif path_matrix[i][j] == "UP":
                seq_real = "_" + seq_real
                seq_created = new_database[curr_HSP[index - 2] + i] + seq_created
                i -= 1
            elif path_matrix[i][j] == "LEFT":
                seq_real = query[curr_HSP[index] + j] + seq_real
                seq_created = "_" + seq_created
                j -= 1
    return seq_real, seq_created


def Needleman_Wunsch(HSP, query, prob_nucleotide, new_database):
    for curr_HSP in HSP:
        l_seq_real = ""
        r_seq_real = ""
        l_seq_created = ""
        r_seq_created = ""
        l_score_max = 0
        r_score_max = 0

        if curr_HSP[2] != 0:
            left = True
            index = 2
            if curr_HSP[0] - 4 * (curr_HSP[index] + 1) < 0:
                score_matrix = np.zeros((curr_HSP[index - 2] + 1, curr_HSP[index] + 1))
                path_matrix = np.zeros((curr_HSP[index - 2] + 1, curr_HSP[index] + 1))
            else:
                score_matrix = np.zeros(
                    (4 * curr_HSP[index - 2] + 1, curr_HSP[index] + 1)
                )
                path_matrix = np.zeros(
                    (4 * curr_HSP[index - 2] + 1, curr_HSP[index] + 1)
                )

            score_matrix, path_matrix = forward_trace(
                score_matrix, path_matrix, query, prob_nucleotide, curr_HSP, left
            )
            l_seq_created, l_seq_real = back_trace(
                score_matrix, path_matrix, query, new_database, curr_HSP, left
            )

        if len(query) != curr_HSP[3] + 1:
            left = False
            index = 3
            if 4 * (len(query) - curr_HSP[3]) + curr_HSP[1] < len(new_database):
                score_matrix = np.zeros(
                    (
                        len(new_database) - curr_HSP[index - 2] - 1,
                        len(query) - curr_HSP[index],
                    )
                )
                path_matrix = np.zeros(
                    (
                        len(new_database) - curr_HSP[index - 2] - 1,
                        len(query) - curr_HSP[index],
                    )
                )
            else:
                score_matrix = np.zeros(
                    (4 * (len(query) - curr_HSP[index]), len(query) - curr_HSP[index])
                )
                path_matrix = np.zeros(
                    (4 * (len(query) - curr_HSP[index]), len(query) - curr_HSP[index])
                )

            score_matrix, path_matrix = forward_trace(
                score_matrix, path_matrix, query, prob_nucleotide, curr_HSP, left
            )

            r_seq_created, r_seq_real = back_trace(
                score_matrix, query, path_matrix, new_database, curr_HSP, left
            )

        seq_real = l_seq_real + new_database[curr_HSP[0]: curr_HSP[1] + 1] + r_seq_real
        seq_created = (
                l_seq_created + query[curr_HSP[2]: curr_HSP[3] + 1] + r_seq_created
        )
        score = HSP[curr_HSP] + r_score_max + l_score_max
        return seq_real, seq_created, score
