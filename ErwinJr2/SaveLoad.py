"""
This file defines functions to save and load JSON files from ErwinJr
"""

from .QCLayers import QCLayers
from .OptStrata import OptStrata
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


def optLoad(fhandle: typing.TextIO) -> OptStrata:
    """
    Load OptStrata from a json file

    Parameters
    ----------
    fhandle : file handle
        file handle of a json file to read in

    Returns
    -------
    strata : OptStrata.OptStrata
        The OptStrata class described in the json file

    Examples
    --------
    >>> import SaveLoad
    >>> with open("path/to/file.json") as f:
    >>>     strata = SaveLoad.optLoad(f)

    """
    ldict = json.load(fhandle)
    return parseStrata(ldict)


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
    if version >= 181107 and version < 210330:
        # Don't touch this for compatibility
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
    elif version >= 210330:
        discription = ldict["Description"]
        ldict = ldict["QCLayers"]
        o = QCLayers(substrate=ldict["Substrate"],
                     materials=ldict["MaterialDefs"]["Compostion"],
                     moleFracs=ldict["MaterialDefs"]["Mole Fraction"],
                     xres=ldict["x resolution"],
                     Eres=ldict["E resolution"],
                     statePerRepeat=ldict["No of states"],
                     layerWidths=ldict["Width"],
                     layerMtrls=ldict["Material"],
                     layerDopings=ldict["Doping"],
                     layerARs=ldict["Active Region"],
                     EField=ldict["EField"],
                     repeats=ldict["Repeats"],
                     T=ldict["Temperature"],
                     solver=ldict["Solver"],
                     description=discription,
                     wl=ldict["Wavelength"])
        if ldict["IFR"]:
            o.includeIFR = True
            o.customIFR = ldict["IFR"]["custom IFR"]
            o.mtrlIFRDelta = ldict["IFR"]["material IFR delta"]
            o.mtrlIFRLambda = ldict["IFR"]["material IFR lambda"]
            o.ifrDelta = ldict["IFR"]["layer IFR delta"]
            o.ifrLambda = ldict["IFR"]["layer IFR lambda"]
        else:
            o.includeIFR = False
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


# Compostion is misspelled for legacy reason
JSONTemplate = """{
    "FileType": "ErwinJr2 Data File",
    "Version": "210330",
    "Description": %s,
    "QCLayers": {
        "Wavelength": %s,
        "Substrate": %s,
        "EField": %s,
        "x resolution": %s,
        "E resolution": %s,
        "No of states": %s,
        "Solver": %s,
        "Temperature": %s,
        "Repeats": %s,
        "MaterialDefs": {
            "Compostion": %s,
            "Mole Fraction": %s
        },
        "Material": %s,
        "Width": %s,
        "Doping": %s,
        "Active Region": %s,
        "IFR": %s
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

IFRSettings = """{
            "custom IFR": %s,
            "material IFR delta": %s,
            "material IFR lambda": %s,
            "layer IFR delta": %s,
            "layer IFR lambda": %s
        }"""


def EJSaveJSON(fhandle: typing.TextIO, qclayers: QCLayers = None,
               optstratum: OptStrata = None):
    """Save QCLayers and OptStratum as a json file

    Parameters
    ----------
    fhandle : file handle
        File handle of a json file to save to.
    qclayers : QCLayers.QCLayers
        The QCLayers class to be saved.
    optstratum : OptStrata.OptStratum
        The OptStratum class to be saved
    """
    if optstratum is None and qclayers is None:
        raise TypeError("Nothing to save")
    if optstratum is None:
        optstratum = OptStrata()
    if qclayers is None:
        qclayers = QCLayers()
    if not (isinstance(qclayers, QCLayers)
            and isinstance(optstratum, OptStrata)):
        raise TypeError("qclSave: Nothing to save.."
                        "QCLayers or OptStratum not valid type")
    o = qclayers
    s = optstratum
    s_cstmtrl = defaultdict(dict)
    for item in s.cstmIndx:
        s_cstmtrl[item]["index"] = str(s.cstmIndx[item])
        if item in s.cstmPrd:
            s_cstmtrl[item]["period"] = s.cstmPrd[item]
        if item in s.cstmGain:
            s_cstmtrl[item]["gain"] = s.cstmGain[item]
    if o.includeIFR:
        ifrParams = IFRSettings % tuple([json.dumps(s) for s in (
            o.customIFR, o.mtrlIFRDelta, o.mtrlIFRLambda,
            o.ifrDelta, o.ifrLambda)])
    else:
        ifrParams = 'false'
    parameters = [json.dumps(s) for s in (o.description, o.wl, o.substrate,
                                          o.EField, o.xres, o.Eres,
                                          o.statePerRepeat, o.solver,
                                          o.temperature, o.repeats,
                                          o.materials, o.moleFracs,
                                          o.layerMtrls, o.layerWidths,
                                          o.layerDopings, o.layerARs)]
    parameters.append(ifrParams)
    parameters += [json.dumps(s) for s in (s.wl, s.materials, s.moleFracs,
                                           s.dopings, list(s.Ls), s.mobilities,
                                           s_cstmtrl)]
    fhandle.write(JSONTemplate % tuple(parameters))
