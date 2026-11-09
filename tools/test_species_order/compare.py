#!/usr/bin/env python3
"""Check that partial RDF and SQ outputs are invariant to trajectory atom order."""
from pathlib import Path
import argparse
import re

import numpy as np


def load(path, names):
    with path.open() as stream:
        header = stream.readline()
    tokens = header.lstrip('#').split()
    values = np.loadtxt(path, comments='#')
    assert values.ndim == 2 and values.shape[1] == len(tokens), path
    return values, {match.group(): index for index, token in enumerate(tokens)
                    if (match := re.search(r'(?:g|S)_\d\d', token))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root', type=Path, nargs='?', default=Path(__file__).resolve().parent)
    parser.add_argument('--expect-swapped', action='store_true', help='expected for the old binary')
    args = parser.parse_args()
    for filename, labels in [('gmixsim.dat', ('g_11', 'g_22')), ('sq.dat', ('S_11', 'S_22'))]:
        a, ai = load(args.root / 'a' / filename, labels)
        b, bi = load(args.root / 'b' / filename, labels)
        assert all(label in ai and label in bi for label in labels), f'{filename}: missing partial headers'
        assert np.allclose(a[:, 0], b[:, 0]), f'{filename}: different grid'
        for label in labels:
            match = np.allclose(a[:, ai[label]], b[:, bi[label]], rtol=1e-3, atol=2e-4)
            if args.expect_swapped:
                assert not match, f'{filename}: {label} unexpectedly agrees with old binary'
            else:
                assert match, f'{filename}: {label} differs between atom orders'
        print(f'{filename}: {"swapped as expected" if args.expect_swapped else "partials agree"}')

    if not args.expect_swapped:
        rdf, idx = load(args.root / 'a/gmixsim.dat', ('g_11', 'g_22'))
        contact = (rdf[:, 0] >= 0.9) & (rdf[:, 0] <= 1.5)
        assert contact.any()
        peak11 = np.nanmax(rdf[contact, idx['g_11']])
        peak22 = np.nanmax(rdf[contact, idx['g_22']])
        assert peak11 > peak22, f'contact peaks: g_11={peak11}, g_22={peak22}'
        print(f'Type 1 contact peak dominates: g_11={peak11:.3f} > g_22={peak22:.3f}')


if __name__ == '__main__':
    main()
