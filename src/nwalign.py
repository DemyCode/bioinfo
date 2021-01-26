#!/usr/bin/python

import argparse
import re
import numpy as np
from blossum62 import blossum62


def dict_ind(letter, sub_dict):
    for i, cur_lett in enumerate(sub_dict):
        if letter == cur_lett:
            return i


def needlewunschman(cmd, x, y, gamma_e, gamma_o, mode):
    score_matrix = np.zeros((len(x) + 1, len(y) + 1))
    backtrack_matrix = np.zeros((len(x) + 1, len(y) + 1, 2), dtype='int')
    substitution_matrix = np.identity(4) + (np.identity(4) - 1)
    sub_dict = ['A', 'T', 'C', 'G']
    if mode == 'rna':
        sub_dict = ['A', 'U', 'C', 'G']
    if mode == 'protein':
        sub_dict = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P ', 'S',
                    'T', 'W', 'Y', 'V', 'B', 'Z', 'X']
        substitution_matrix = blossum62()
    for i in range(score_matrix.shape[0]):
        score_matrix[i, 0] = i * gamma_e
    for j in range(score_matrix.shape[1]):
        score_matrix[0, j] = j * gamma_e
    for i in range(1, score_matrix.shape[0]):
        for j in range(1, score_matrix.shape[1]):
            # For x and y the indices are one index ahead so x - 1 instead of x, same for y
            substitution = substitution_matrix[dict_ind(x[i - 1], sub_dict), dict_ind(y[j - 1], sub_dict)] + \
                           score_matrix[i - 1, j - 1]
            insertion = gamma_e + score_matrix[i - 1, j]
            deletion = gamma_e + score_matrix[i, j - 1]
            score_matrix[i, j] = substitution
            backtrack_matrix[i, j] = [i - 1, j - 1]
            if score_matrix[i, j] < insertion:
                score_matrix[i, j] = insertion
                backtrack_matrix[i, j] = [i - 1, j]
            if score_matrix[i, j] < deletion:
                score_matrix[i, j] = deletion
                backtrack_matrix[i, j] = [i, j - 1]
    return score_matrix, backtrack_matrix


def matrix_treatment(cmd, score_matrix, backtrack_matrix):
    if cmd == 'score':
        print(score_matrix[-1, -1])
        return
    u, v = score_matrix.shape[0] - 1, score_matrix.shape[1] - 1
    resx = ''
    resy = ''
    while u != 0 and v != 0:
        newu, newv = backtrack_matrix[u, v]
        if newu == u - 1 and newv == v - 1:
            resx = x[u - 1] + resx
            resy = y[v - 1] + resy
        if newu == u - 1 and newv == v:
            resx = x[u - 1] + resx
            resy = '-' + resy
        if newu == u and newv == v - 1:
            resx = '-' + resx
            resy = y[v - 1] + resy
        u, v = newu, newv
    newu, newv = backtrack_matrix[u, v]
    if newu == u - 1 and newv == v - 1:
        resx = x[u - 1] + resx
        resy = y[v - 1] + resy
    if newu == u - 1 and newv == v:
        resx = x[u - 1] + resx
        resy = '-' + resy
    if newu == u and newv == v - 1:
        resx = '-' + resx
        resy = y[v - 1] + resy
    u, v = newu, newv
    print(resx)
    print(resy)


def main(cmd, x, y, gamma_e, gamma_o, mode):
    score_matrix, backtrack_matrix = needlewunschman(cmd, x, y, gamma_e, gamma_o, mode)
    # print(backtrack_matrix)
    matrix_treatment(cmd, score_matrix, backtrack_matrix)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gamma', nargs=2, type=int)
    parser.add_argument('cmd')
    parser.add_argument('x')
    parser.add_argument('y')
    arguments = parser.parse_args()
    cmd = arguments.cmd
    x = arguments.x
    y = arguments.y
    gamma_e = -1
    gamma_o = 0
    if arguments.gamma is not None:
        gamma_e, gamma_o = arguments.gamma

    # First Check
    if cmd not in ['score', 'align']:
        raise RuntimeError('cmd : {} is unrecognized'.format(cmd))

    # Second Check
    mode = None
    pattern = re.compile(r'[ARNDCQEGHILKMFPSTWYVBZX]+')
    if pattern.fullmatch(x) and pattern.fullmatch(y):
        mode = 'protein'
    pattern = re.compile(r'[AUCG]+')
    if pattern.fullmatch(x) and pattern.fullmatch(y):
        mode = 'rna'
    pattern = re.compile(r'[ATCG]+')
    if pattern.fullmatch(x) and pattern.fullmatch(y):
        mode = 'dna'
    if mode is None:
        raise RuntimeError('x y : {} {} doesnt match any mode'.format(x, y))

    # Third Check
    if arguments.gamma is not None and (gamma_e > 0 or gamma_o > 0):
        raise RuntimeError('gamma o : {} {} invalid value'.format(gamma_e, gamma_o))

    main(cmd, x, y, gamma_e, gamma_o, mode)