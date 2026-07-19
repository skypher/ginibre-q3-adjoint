#!/usr/bin/env python3
"""Read a moments-file from the C++ moments binary and verify:
1. Values match OEIS b-files where overlapping (sanity check).
2. Chain Inequality Q_3(2m+3) >= 4*Q_3(2m+1) for all m up to the computed range.

Usage:
  python3 verify_chain.py GROUP LOGFILE
"""
import re
import sys
from math import comb


OEIS_FILES = {
    'G2': 'oeis_A227292.txt',
    'F4': 'oeis_A179685.txt',
    'E6': 'oeis_A179684.txt',
    'E7': 'oeis_A179683.txt',
    'E8': 'oeis_A179663.txt',
}


def parse_moments(logfile):
    moms = {}
    with open(logfile) as f:
        for line in f:
            m = re.search(r'\bm_(\d+)\s*=\s*(-?\d+)', line)
            if m:
                index = int(m.group(1))
                value = int(m.group(2))
                if index in moms and moms[index] != value:
                    raise ValueError(
                        f'conflicting values for m_{index} in {logfile}: '
                        f'{moms[index]} and {value}'
                    )
                moms[index] = value
    if not moms:
        raise ValueError(f'no moment lines found in {logfile}')
    maximum = max(moms)
    missing = sorted(set(range(maximum + 1)) - moms.keys())
    if missing:
        raise ValueError(f'missing moment indices in {logfile}: {missing}')
    return [moms[k] for k in range(maximum + 1)]


def parse_oeis(path):
    moms = {}
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            k, v = s.split()
            index = int(k)
            value = int(v)
            if index in moms and moms[index] != value:
                raise ValueError(
                    f'conflicting values for OEIS index {index} in {path}: '
                    f'{moms[index]} and {value}'
                )
            moms[index] = value
    return moms


def Q3(n, mom):
    return 2 * sum(comb(n, k) * (mom[k + 2] * mom[n - k] - mom[k + 1] * mom[n - k + 1])
                   for k in range(n + 1))


def main():
    if len(sys.argv) not in (3, 4):
        print('usage: verify_chain.py GROUP LOGFILE [OEIS_DIR]')
        sys.exit(1)
    group = sys.argv[1]
    logfile = sys.argv[2]
    oeis_dir = sys.argv[3] if len(sys.argv) > 3 else '../references'

    moms = parse_moments(logfile)
    K = len(moms) - 1
    if K < 5:
        raise ValueError(
            f'{logfile} ends at m_{K}; at least m_0..m_5 are required '
            'to verify one Chain inequality'
        )
    print(f'{group}: {K + 1} moments read from {logfile} (m_0..m_{K})')

    # OEIS cross-check
    if group in OEIS_FILES:
        try:
            oeis = parse_oeis(f'{oeis_dir}/{OEIS_FILES[group]}')
            overlap = [k for k in oeis if k <= K]
            if not overlap:
                raise ValueError(
                    f'no overlap with OEIS b-file {OEIS_FILES[group]}'
                )
            mismatches = [k for k in overlap if moms[k] != oeis[k]]
            print(f'  OEIS overlap: m_0..m_{max(overlap)}; {len(mismatches)} mismatches.')
            if mismatches:
                for k in mismatches[:5]:
                    print(f'    m_{k}: computed {moms[k]}, OEIS {oeis[k]}')
                sys.exit(1)
            extended = [k for k in range(K + 1) if k not in oeis]
            if extended:
                print(f'  NEW beyond OEIS: m_{min(extended)}..m_{max(extended)} '
                      f'({len(extended)} new values).')
        except FileNotFoundError:
            print(f'  (OEIS b-file {OEIS_FILES[group]} not found; skipping cross-check)')

    # Chain verification
    max_m = (K - 5) // 2
    print(f'Chain m verifiable up to m = {max_m}:')
    all_ok = True
    for mb in range(0, max_m + 1):
        q1 = Q3(2 * mb + 1, moms)
        q3 = Q3(2 * mb + 3, moms)
        diff = q3 - 4 * q1
        ok = diff >= 0
        if not ok:
            all_ok = False
        q1_ok = q1 >= 0
        if not q1_ok:
            all_ok = False
        if mb % 5 == 0 or mb == max_m or not ok or not q1_ok:
            status = 'OK' if ok else 'FAIL'
            q1status = '' if q1_ok else ' Q3<0!'
            print(f'  m={mb:3d}: Q_3({2*mb+3}) - 4*Q_3({2*mb+1}) = {diff} [{status}]{q1status}')
    print(f'Chain verification: {"ALL PASS" if all_ok else "SOME FAIL"}')
    if not all_ok:
        sys.exit(1)


if __name__ == '__main__':
    main()
