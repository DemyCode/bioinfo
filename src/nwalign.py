#!/usr/bin/python

import argparse
import re
import numpy as np
from blossum62 import blossum62


def dict_ind(letter, sub_dict):
    for i, cur_lett in enumerate(sub_dict):
        if letter == cur_lett:
            return i


def step1(cmd, x, y, gamma_e, gamma_o, mode):
    score_matrix = np.zeros((len(x), len(y)))
    substitution_matrix = np.identity(4)
    substitution_matrix = substitution_matrix + (substitution_matrix - 1)
    # print(substitution_matrix)
    sub_dict = ['A', 'T', 'C', 'G']
    if cmd == 'rna':
        sub_dict = ['A', 'U', 'C', 'G']
    if cmd == 'protein':
        sub_dict = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P ', 'S',
                    'T', 'W', 'Y', 'V', 'B', 'Z', 'X']
    if cmd == 'protein':
        substitution_matrix = blossum62()
    for i in range(score_matrix.shape[0]):
        score_matrix[i, 0] = i * gamma_e
    for j in range(score_matrix.shape[1]):
        score_matrix[0, j] = j * gamma_e
    for i in range(1, score_matrix.shape[0]):
        for j in range(1, score_matrix.shape[1]):
            substitution = substitution_matrix[dict_ind(x[i], sub_dict), dict_ind(y[j], sub_dict)] + score_matrix[i - 1, j - 1]
            insertion = gamma_e + score_matrix[i - 1, j]
            deletion = gamma_e + score_matrix[i, j - 1]
            score_matrix[i, j] = max(substitution, insertion, deletion)
    return score_matrix


def main(cmd, x, y, gamma_e, gamma_o, mode):
    print(step1(cmd, x, y, gamma_e, gamma_o, mode)[-1, -1] + 1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gamma', nargs=2)
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
    if cmd not in ['score', 'nwalign']:
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