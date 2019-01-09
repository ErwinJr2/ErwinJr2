#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
need to install pycallgraph package and GraphViz
"""

from context import *
import SaveLoad
import numpy as np
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

if __name__ == "__main__":

    with PyCallGraph(output=GraphvizOutput()):
        
        with open("../example/16um.json") as f:
            qcl = SaveLoad.qclLoad(f)

        qcl.layerSelected = 3
        qcl.NonParabolic = False
        qcl.populate_x()
        qcl.solve_whole()
        qcl.dipole(19, 15)
        qcl.FoM(19, 15)
