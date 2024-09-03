#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import cProfile
import io
import pstats

import SaveLoad
from QCLayers import auto_gain

if __name__ == "__main__":
    cp = cProfile.Profile()
    with open("../example/PQLiu.json") as f:
        qcl = SaveLoad.qclLoad(f)
    qcl.xres = 0.1
    cp.enable()
    auto_gain(qcl)
    cp.disable()
    cp.print_stats(sort="cumulative")

    s = io.StringIO()
    ps = pstats.Stats(cp, stream=s).sort_stats('cumulative')
    ps.print_stats()

    with open('profiling_0.txt', 'w+') as f:
        f.write(s.getvalue())
