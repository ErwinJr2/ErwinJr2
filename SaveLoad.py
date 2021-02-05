"""
This file defines functions to save and load JSON files from ErwinJr
"""

from QCLayers import QCLayers
from OptStrata import OptStrata
from collections import defaultdict
import json
import typing


def qclLoad(fhandle: typing.TextIO) -> QCLayers:
    """
    Load QCLayers from a json file

    Parameters
    ----------
    fhandle : file handle
        file handle of a json file to read in

    Returns
    -------
    qclayers : QCLayers.QCLayers
        The QCLayers class described in the json file

    Examples
    --------
    >>> import SaveLoad
    >>> with open("path/to/file.json") as f:
    >>>     qcl = SaveLoad.qclLoad(f)

    """
    ldict = json.load(fhandle)
    return parseQcl(ldict)


def loadBoth(fhandle: typing.TextIO) -> typing.Union[QCLayers, OptStrata]:
    """
    Load QCLayers and OptStrata from a json file

    Parameters
    ----------
    fhandle : file handle
        file handle of a json file to read in

    Returns
    -------
    qclayers : QCLayers.QCLayers
        The QCLayers class described in the json file
    strata : OptStrata.OptStrata or None
        The Optical Strata class described in the json file. For older version
        this is None

    Examples
    --------
    >>> import SaveLoad
    >>> with open("path/to/file.json") as f:
    >>>     qcl, stratum = SaveLoad.loadBoth(f)

    """
    ldict = json.load(fhandle)
    qcl = parseQcl(ldict)
    try:
        stratum = parseStrata(ldict)
    except NotImplementedError:
        stratum = None
    return qcl, stratum


def parseQcl(ldict: typing.Dict[str, typing.Any]) -> QCLayers:
    if ldict["FileType"] != "ErwinJr2 Data File":
        raise TypeError("Wrong file type")
    version = int(ldict["Version"])
    if version >= 181107:
        o = QCLayers(substrate=ldict["Substrate"],
                     materials=ldict["Materials"]["Compostion"],
                     moleFracs=ldict["Materials"]["Mole Fraction"],
                     xres=ldict["x resolution"],
                     Eres=ldict["E resolution"],
                     layerWidths=ldict["QC Layers"]["Width"],
                     layerMtrls=ldict["QC Layers"]["Material"],
                     layerDopings=ldict["QC Layers"]["Doping"],
                     layerARs=ldict["QC Layers"]["Active Region"],
                     EField=ldict["EField"],
                     repeats=ldict["Repeats"],
                     T=ldict["Temperature"],
                     solver=ldict["Solver"],
                     description=ldict["Description"])
        o.wl = ldict["Wavelength"] if "Wavelength" in ldict else 1.5
        if version >= 210205:
            o.customIFR = ldict["IFR"]["custom IFR"]
            o.mtrlIFRDelta = ldict["IFR"]["material IFR delta"]
            o.mtrlIFRLambda = ldict["IFR"]["material IFR lambda"]
            o.ifrDelta = ldict["IFR"]["layer IFR delta"]
            o.ifrLambda = ldict["IFR"]["layer IFR lambda"]
    else:
        raise NotImplementedError("Version %s not supported" %
                                  ldict["Version"])
    return o


def parseStrata(ldict: typing.Dict[str, typing.Any]) -> OptStrata:
    if ldict["FileType"] != "ErwinJr2 Data File":
        raise TypeError("Wrong file type")
    if int(ldict["Version"]) >= 200504:
        ldict = ldict["Waveguide"]
        cstidx = {}
        cstprd = {}
        cstgain = {}
        for item in ldict["custom"]:
            cstidx[item] = complex(ldict["custom"][item]["index"])
            if "period" in ldict["custom"][item]:
                cstprd[item] = ldict["custom"][item]["period"]
            if "gain" in ldict["custom"][item]:
                cstgain[item] = ldict["custom"][item]["gain"]
        o = OptStrata(wl=ldict["wavelength"],
                      materials=ldict["materials"],
                      moleFracs=ldict["moleFracs"],
                      dopings=ldict["dopings"],
                      Ls=ldict["width"],
                      mobilities=ldict["mobilities"],
                      cstmIndx=cstidx, cstmPrd=cstprd, cstmGain=cstgain)
    else:
        raise NotImplementedError("Version %s not supported" %
                                  ldict["Version"])
    return o


JSONTemplate = """{
    "FileType": "ErwinJr2 Data File",
    "Version": "210205",
    "Description": %s,
    "Wavelength": %s,
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
    "IFR": {
        "custom IFR": %s,
        "material IFR delta": %s,
        "material IFR lambda": %s,
        "layer IFR delta": %s,
        "layer IFR lambda": %s
    },
    "QC Layers": {
        "Material": %s,
        "Width": %s,
        "Doping": %s,
        "Active Region": %s
    },
    "Waveguide": {
        "wavelength": %s,
        "materials": %s,
        "moleFracs": %s,
        "dopings": %s,
        "width": %s,
        "mobilities": %s,
        "custom": %s
    }
}"""


def EJSaveJSON(fhandle: typing.TextIO, qclayers: QCLayers,
               optstratum: OptStrata):
    """Save QCLayers and OptStratum as a json file

    Parameters
    ----------
    fhandle : file handle
        File handle of a json file to save to.
    qclayers : QCLayers.QCLayers
        The QCLayers class to be saved.
    optstratum : OptStrata.OptStratum
        The OptStratum class to be saced
    """
    if not (isinstance(qclayers, QCLayers)
            and isinstance(optstratum, OptStrata)):
        raise TypeError("qclSave: Nothing to save.."
                        "QCLayers or OptStratum not valid type")
    o = qclayers
    s = optstratum
    cstmtrl = defaultdict(dict)
    for item in s.cstmIndx:
        cstmtrl[item]["index"] = str(s.cstmIndx[item])
        if item in s.cstmPrd:
            cstmtrl[item]["period"] = s.cstmPrd[item]
        if item in s.cstmGain:
            cstmtrl[item]["gain"] = s.cstmGain[item]
    parameters = [json.dumps(s) for s in (o.description, o.wl, o.substrate,
                                          o.EField, o.xres, o.Eres, o.solver,
                                          o.Temperature, o.repeats,
                                          o.materials, o.moleFracs,
                                          o.customIFR,
                                          o.mtrlIFRDelta, o.mtrlIFRLambda,
                                          o.ifrDelta, o.ifrLambda,
                                          o.layerMtrls, o.layerWidths,
                                          o.layerDopings, o.layerARs,
                                          s.wl, s.materials, s.moleFracs,
                                          s.dopings, list(s.Ls), s.mobilities)]
    parameters.append(json.dumps(cstmtrl))
    # parameters.append(json.dumps(cstmtrl, indent=4).replace('\n','\n'+' '*8))
    fhandle.write(JSONTemplate % tuple(parameters))

# vim: ts=4 sw=4 sts=4 expandtab
