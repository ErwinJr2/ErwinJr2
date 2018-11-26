from context import *
import QCLayers as qcl
import SaveLoad

with open("../example/PQLiu.json") as f:
     qcl = SaveLoad.qclLoad(f)
