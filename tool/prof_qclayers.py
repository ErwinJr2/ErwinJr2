#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import cProfile
import io
import pstats

import ErwinJr2.save_load as save_load
from ErwinJr2.qc_layers import auto_gain

if __name__ == "__main__":
    cp = cProfile.Profile()
    with open("../example/PQLiu.json") as f:
        qcl = save_load.qcl_load(f)
    qcl.x_step = 0.1
    cp.enable()
    auto_gain(qcl)
    cp.disable()
    cp.print_stats(sort="cumulative")

    s = io.StringIO()
    ps = pstats.Stats(cp, stream=s).sort_stats("cumulative")
    ps.print_stats()

    with open("profiling_0.txt", "w+") as f:
        f.write(s.getvalue())
