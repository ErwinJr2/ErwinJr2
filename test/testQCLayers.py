from context import *
import SaveLoad

with open("../example/PQLiu.json") as f:
     qcl = SaveLoad.qclLoad(f)

qcl.populate_x()
