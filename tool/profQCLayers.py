#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
import SaveLoad
import cProfile
import pstats
import io

if __name__ == "__main__":
    cp = cProfile.Profile()

    cp.enable()
    with open("../example/PQLiu.json") as f:
        qcl = SaveLoad.qclLoad(f)

    qcl.layerSelected = 3
    #  qcl.NonParabolic = False
    qcl.populate_x()
    qcl.solve_whole()
    qcl.dipole(19, 15)
    qcl.calc_FoM(19, 15)
    cp.disable()

    cp.print_stats(sort="cumulative")

    # s = io.StringIO()
    # ps = pstats.Stats(cp, stream=s).sort_stats('cumulative')
    # ps.print_stats()

    # with open('profiling_0.txt', 'w+') as f:
    #     f.write(s.getvalue())
