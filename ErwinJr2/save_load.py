"""
This file defines functions to save and load JSON files from ErwinJr
"""

import json
import typing
from collections import defaultdict

from ErwinJr2.opt_strata import OptStrata
from ErwinJr2.qc_layers import QCLayers


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
        o = QCLayers(
            substrate=ldict["Substrate"],
            materials=ldict["Materials"]["Compostion"],
            mole_fracs=ldict["Materials"]["Mole Fraction"],
            x_res=ldict["x resolution"],
            e_res=ldict["E resolution"],
            layer_widths=ldict["QC Layers"]["Width"],
            layer_matrls=ldict["QC Layers"]["Material"],
            layer_dopings=ldict["QC Layers"]["Doping"],
            layer_ar=ldict["QC Layers"]["Active Region"],
            e_field=ldict["EField"],
            repeats=ldict["Repeats"],
            temp=ldict["Temperature"],
            solver=ldict["Solver"],
            description=ldict["Description"],
        )
        o.wl = ldict["Wavelength"] if "Wavelength" in ldict else 1.5
    elif version >= 210330:
        discription = ldict["Description"]
        ldict = ldict["QCLayers"]
        o = QCLayers(
            substrate=ldict["Substrate"],
            materials=ldict["MaterialDefs"]["Compostion"],
            mole_fracs=ldict["MaterialDefs"]["Mole Fraction"],
            x_res=ldict["x resolution"],
            e_res=ldict["E resolution"],
            state_per_repeat=ldict["No of states"],
            layer_widths=ldict["Width"],
            layer_matrls=ldict["Material"],
            layer_dopings=ldict["Doping"],
            layer_ar=ldict["Active Region"],
            e_field=ldict["EField"],
            repeats=ldict["Repeats"],
            temp=ldict["Temperature"],
            solver=ldict["Solver"],
            description=discription,
            wl=ldict["Wavelength"],
        )
        if ldict["IFR"]:
            o.include_ifr = True
            o.custom_ifr = ldict["IFR"]["custom IFR"]
            o.mtrl_ifr_delta = ldict["IFR"]["material IFR delta"]
            o.mtrl_ifr_lambda = ldict["IFR"]["material IFR lambda"]
            o.ifr_deltas = ldict["IFR"]["layer IFR delta"]
            o.ifr_lambdas = ldict["IFR"]["layer IFR lambda"]
        else:
            o.include_ifr = False
    else:
        raise NotImplementedError("Version %s not supported" % ldict["Version"])
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
        o = OptStrata(
            wl=ldict["wavelength"],
            materials=ldict["materials"],
            moleFracs=ldict["moleFracs"],
            dopings=ldict["dopings"],
            Ls=ldict["width"],
            mobilities=ldict["mobilities"],
            cstmIndx=cstidx,
            cstmPrd=cstprd,
            cstmGain=cstgain,
        )
    else:
        raise NotImplementedError("Version %s not supported" % ldict["Version"])
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


def EJSaveJSON(
    fhandle: typing.TextIO, qclayers: QCLayers = None, optstratum: OptStrata = None
):
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
    if not (isinstance(qclayers, QCLayers) and isinstance(optstratum, OptStrata)):
        raise TypeError(
            "qclSave: Nothing to save.." "QCLayers or OptStratum not valid type"
        )
    o = qclayers
    s = optstratum
    s_cstmtrl = defaultdict(dict)
    for item in s.cstmIndx:
        s_cstmtrl[item]["index"] = str(s.cstmIndx[item])
        if item in s.cstmPrd:
            s_cstmtrl[item]["period"] = s.cstmPrd[item]
        if item in s.cstmGain:
            s_cstmtrl[item]["gain"] = s.cstmGain[item]
    if o.include_ifr:
        ifrParams = IFRSettings % tuple(
            [
                json.dumps(s)
                for s in (
                    o.custom_ifr,
                    o.mtrl_ifr_delta,
                    o.mtrl_ifr_lambda,
                    o.ifr_deltas,
                    o.ifr_lambdas,
                )
            ]
        )
    else:
        ifrParams = "false"
    parameters = [
        json.dumps(s)
        for s in (
            o.description,
            o.wl,
            o.substrate,
            o.e_field,
            o.x_step,
            o.e_step,
            o.state_per_repeat,
            o.solver,
            o.temperature,
            o.repeats,
            o.materials,
            o.mole_fracs,
            o.layer_mtrls,
            o.layer_widths,
            o.layer_dopings,
            o.layer_ar,
        )
    ]
    parameters.append(ifrParams)
    parameters += [
        json.dumps(s)
        for s in (
            s.wl,
            s.materials,
            s.moleFracs,
            s.dopings,
            list(s.Ls),
            s.mobilities,
            s_cstmtrl,
        )
    ]
    fhandle.write(JSONTemplate % tuple(parameters))
