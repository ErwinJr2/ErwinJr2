#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
import SaveLoad
import cProfile
import pstats
import io

if __name__ == "__main__":
    cp = cProfile.Profile()
    with open("../example/PQLiu.json") as f:
        qcl = SaveLoad.qclLoad(f)
    qcl.xres = 0.1
    cp.enable()
    qcl._auto_gain()
    cp.disable()
    cp.print_stats(sort="cumulative")

    s = io.StringIO()
    ps = pstats.Stats(cp, stream=s).sort_stats('cumulative')
    ps.print_stats()

    with open('profiling_0.txt', 'w+') as f:
        f.write(s.getvalue())
