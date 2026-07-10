"""
observables.py

Physics observables for the CERN Quantum Twin challenge.

Compatible with Pulser >=1.8

Works directly from Pulser Results.final_bitstrings
or Results.final_state.
"""

from __future__ import annotations

from collections import Counter

import numpy as np


# ---------------------------------------------------
# Bitstring helpers
# ---------------------------------------------------

def bitstring_to_sz(bitstring):
    """
    Convert a Pulser bitstring to σz values.

    Pulser:
        r = 1
        g = 0

    We use

        σz = +1 (ground)

        σz = -1 (rydberg)
    """

    if isinstance(bitstring, str):
        bits = np.array(
            [int(x) for x in bitstring],
            dtype=int,
        )

    else:
        bits = np.asarray(bitstring)

    return 1 - 2 * bits


# ---------------------------------------------------
# Checkerboard mask
# ---------------------------------------------------

def checkerboard_mask(n_side):

    return np.array(
        [
            [(-1) ** (i + j) for j in range(n_side)]
            for i in range(n_side)
        ]
    )


# ---------------------------------------------------
# Staggered magnetization
# ---------------------------------------------------

def staggered_from_bitstring(
    bitstring,
    n_side,
):

    sz = bitstring_to_sz(bitstring)

    sz = sz.reshape(n_side, n_side)

    mask = checkerboard_mask(n_side)

    return np.mean(sz * mask)


def staggered_magnetization(
    counter,
    n_side,
):
    """
    Average staggered magnetization over all measured bitstrings.

    counter:
        Pulser Results.final_bitstrings
    """

    total = sum(counter.values())

    m = 0

    for bits, count in counter.items():

        m += (
            staggered_from_bitstring(bits, n_side)
            * count
        )

    return m / total


# ---------------------------------------------------
# Bulk magnetization
# ---------------------------------------------------

def bulk_magnetization(counter):

    total = sum(counter.values())

    M = 0

    for bits, count in counter.items():

        sz = bitstring_to_sz(bits)

        M += sz.mean() * count

    return M / total


# ---------------------------------------------------
# AFM structure factor
# ---------------------------------------------------

def structure_factor(counter, n_side):

    total = sum(counter.values())

    S = 0

    mask = checkerboard_mask(n_side)

    for bits, count in counter.items():

        sz = bitstring_to_sz(bits)

        sz = sz.reshape(n_side, n_side)

        value = np.mean(sz * mask)

        S += value ** 2 * count

    return S / total


# ---------------------------------------------------
# Defect density
# ---------------------------------------------------

def defect_density(
    counter,
    bonds,
):

    total = sum(counter.values())

    rho = 0

    for bits, count in counter.items():

        sz = bitstring_to_sz(bits)

        defects = 0

        for i, j in bonds:

            if np.sign(sz[i]) == np.sign(sz[j]):

                defects += 1

        rho += defects / len(bonds) * count

    return rho / total


# ---------------------------------------------------
# Correlation matrix
# ---------------------------------------------------

def correlation_matrix(counter):

    total = sum(counter.values())

    first = next(iter(counter))

    n = len(first)

    C = np.zeros((n, n))

    for bits, count in counter.items():

        sz = bitstring_to_sz(bits)

        C += count * np.outer(sz, sz)

    return C / total


# ---------------------------------------------------
# Binder cumulant
# ---------------------------------------------------

def binder_cumulant(counter, n_side):

    ms = []

    for bits, count in counter.items():

        m = staggered_from_bitstring(
            bits,
            n_side,
        )

        ms.extend([m] * count)

    ms = np.asarray(ms)

    m2 = np.mean(ms ** 2)

    m4 = np.mean(ms ** 4)

    return 1 - m4 / (3 * m2 ** 2)


# ---------------------------------------------------
# Domain size
# ---------------------------------------------------

def average_domain_size(counter):

    sizes = []

    for bits, count in counter.items():

        sz = bitstring_to_sz(bits)

        runs = []

        run = 1

        for i in range(1, len(sz)):

            if sz[i] == sz[i - 1]:

                run += 1

            else:

                runs.append(run)

                run = 1

        runs.append(run)

        sizes.extend(runs * count)

    return np.mean(sizes)