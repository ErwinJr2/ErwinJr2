#! /usr/bin/env python3
# -*- coding:utf-8 -*-
"""
This file is part of fftautocorr

Copyright (C) 2020 CareF
author: CareF
Licensed under a 3-clause BSD style license - see LICENSE.md

The script generates a table for the optimized padded length for FFT.
Running this script is part of pre-compilation. See Makefile for details.
"""

import os

import numpy as np


def tablegen(intmax, ps=(2, 3, 5)):
    """
    Generate a composite number table < intmax with prime composition in [ps]
    """
    logL = np.log2(intmax)
    allp = 1
    for p in ps:
        compp = p**(np.arange(logL/np.log2(p), dtype=np.int))
        allp = np.outer(allp, compp).reshape(-1)
        allp = allp[allp < intmax]
    return np.sort(allp)


def print_as_c_header(filename, array):
    contents = """/*
 * This file is generated by {this}
 */

#ifndef {filename}
#define {filename}
#include <stdlib.h>

static int find_factor(size_t n) {{
    const static int factortable[] = {{
        {numbers}
    }};
    int lo=0, hi=sizeof(factortable)/sizeof(factortable[0]), mid;
    while(lo < hi) {{
        mid = (lo + hi) / 2;
        if(factortable[mid] < n)
            lo = mid+1;
        else
            hi = mid;
        /* factortable[lo:] >= n; factortable[:lo] < n */
    }}
    if(lo == sizeof(factortable)/sizeof(factortable[0]))
        return -1;
    return factortable[lo];
}}
#endif
""".format(filename=filename.replace('.', '_').upper(),
           this=os.path.basename(__file__),
           numbers=', '.join([str(n) for n in array]))
    with open(filename, 'w') as f:
        f.write(contents)


if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    print("Generating %s" % filename)
    print_as_c_header(filename, tablegen(2**31-1))
