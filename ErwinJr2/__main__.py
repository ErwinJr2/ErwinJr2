#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
Make `python -m ErwinJr2` an alias for running `ErwinJr2`.
"""

import os
import sys

import ErwinJr2.gui.erwinjr as erwinjr


def main():
    try:
        filename = os.path.abspath(sys.argv[1])
    except IndexError:
        filename = None
    exit(erwinjr.main(filename))


if __name__ == "__main__":
    main()
