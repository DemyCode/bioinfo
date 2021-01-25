import argparse
import re


def main(cmd, x, y, gamma_e, gamma_o, mode):
    pass


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
    gamma_e = None
    gamma_o = None
    if arguments.gamma is not None:
        gamma_e, gamma_o = arguments.gamma

    # First Check
    if cmd not in ['score', 'nwalign']:
        raise RuntimeError('cmd : {} is unrecognized'.format(cmd))

    # Second Check
    mode = None
    pattern = re.compile('[ARNDCQEGHILKMFPSTWYVBZX]+')
    if pattern.match(x) and pattern.match(y):
        mode = 'protein'
    pattern = re.compile('[AUCG]+')
    if pattern.match(x) and pattern.match(y):
        mode = 'rna'
    pattern = re.compile('[ATCG]+')
    if pattern.match(x) and pattern.match(y):
        mode = 'dna'
    if mode not in ['protein', 'rna', 'dna']:
        raise RuntimeError('x y : {} {} have unrecognized mode'.format(x, y))

    # Third Check
    if gamma_e is not None and (gamma_e > 0 or gamma_o > 0):
        raise RuntimeError('gamma o : {} {} invalid value'.format(gamma_e, gamma_o))

    main(cmd, x, y, gamma_e, gamma_o, mode)
