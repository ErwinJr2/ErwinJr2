#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
Make `python -m ErwinJr` an alias for running `ErwinJr`.
"""

import os
import sys
from ErwinJr2.ErwinJr import main


def __main__():
    try:
        fileName = os.path.abspath(sys.argv[1])
    except IndexError:
        fileName = None
    main(fileName)


if __name__ == "__main__":
    __main__()
