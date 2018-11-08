#!/usr/bin/env python
# -*- coding:utf-8 -*-

from QCLayers import QCLayers
import sys, json
import numpy as np

def qclLoad(fhandle):
    """ Load QCLayers from a json file """
    ldict = json.load(fhandle)
    if ldict["FileType"] != "ErwinJr2 Data File":
        raise TypeError("Wrong file type")
    if ldict["Version"] == "181107":
        o = QCLayers(ldict["Substrate"],
                     ldict["Materials"]["Compostion"], 
                     ldict["Materials"]["Mole Fraction"], 
                     ldict["x resolution"], 
                     ldict["E resolution"], 
                     ldict["QC Layers"]["Width"], 
                     ldict["QC Layers"]["Material"], 
                     ldict["QC Layers"]["Doping"],
                     ldict["QC Layers"]["Active Region"], 
                     ldict["EField"], 
                     ldict["Repeats"], 
                     ldict["Temperature"], 
                     ldict["Solver"], 
                     ldict["Description"])
    else: 
        raise NotImplementedError("Version %s not supported" %
                                  ldict["Version"])
    return o

JSONTemplate = """{
    "FileType": "ErwinJr2 Data File", 
    "Version": "181107", 
    "Description": %s, 
    "Substrate": %s, 
    "EField": %s, 
    "x resolution": %s, 
    "E resolution": %s, 
    "Solver": %s, 
    "Temperature": %s, 
    "Repeats": %s, 
    "Materials": {
        "Compostion": %s, 
        "Mole Fraction": %s
    }, 
    "QC Layers": {
        "Material": %s, 
        "Width": %s, 
        "Doping": %s, 
        "Active Region": %s
    }
}"""
def qclSaveJSON(fhandle, qclayers):
    "Save file with filename as json"
    if not isinstance(qclayers, QCLayers):
        raise TypeError("qclSave: Nothing to save.."
                        "QCLayers not valid type")
    o=qclayers
    parameters = [json.encoder.encode_basestring(o.description)] + [
        repr(s).replace("'","\"") for s in (
            o.substrate, o.EField, o.xres, o.Eres, o.Solver, 
            o.Temperature, o.repeats, o.materials, o.moleFracs, 
            o.layerMaterials.tolist(), o.layerWidths.tolist(),
            o.layerDopings.tolist(), o.layerARs.tolist())]
    fhandle.write(JSONTemplate % tuple(parameters))

# vim: ts=4 sw=4 sts=4 expandtab
