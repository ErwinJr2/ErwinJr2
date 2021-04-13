"""
This contains the plotter functions for QCLayers
"""
from typing import Iterable, Union
import numpy as np
from matplotlib.axes import Axes
from matplotlib import cm
# the normalization used for state population
# from matplotlib.colors import LogNorm as cmNorm
from matplotlib.colors import Normalize as cmNorm
from .QCLayers import QCLayers

config = {
    "wf_scale": 0.3,
    "mode_scale": 3,
    # This number not necessarily but is chosen to be consistent with the
    # default value of tol for QCLayers.periodRecognize as
    # wf_almost_zero / modescale = tol
    "wf_almost_zero": 1.5e-4,
    "wf_colors": ((0.584, 0.450, 0.701), (0.431, 0.486, 0.745),
                  (0.576, 0.694, 0.517), (0.682, 0.780, 0.321),
                  (0.501, 0.501, 0.509), (0.854, 0.741, 0.247),
                  (0.874, 0.607, 0.290), (0.823, 0.341, 0.278),
                  (0.725, 0.321, 0.623), (0.411, 0.741, 0.270),
                  (0.078, 0.078, 0.078), (0.431, 0.803, 0.870),
                  (0.223, 0.321, 0.643))
}  #: The configuration dictionary for plotting.


def plotPotential(axes: Axes, qcl: QCLayers,
                  plotVL: bool = False, plotVX: bool = False,
                  plotLH: bool = False, plotSO: bool = False):
    """Plot the potentials of qcl.

    Parameters
    ----------
    axes :
        The axes to plot the figure on.
    qcl :
        The QCLayers to plot
    plotVL, plotVX, plotLH, plotSO :
        flags wether to plot L-point, X-point, LH band, SO band

    Returns
    -------
    A list of plotted data

    """
    ys = [qcl.xVc]
    axes.plot(qcl.xPoints, qcl.xVc, 'k', linewidth=1)
    # highlight selected layer & make AR layers bold
    xNonARs = np.bitwise_and.reduce(
        [qcl.xLayerMask(n) for n, ar in enumerate(qcl.layerARs) if ar])
    xARVc = np.ma.masked_where(xNonARs, qcl.xVc)
    axes.plot(qcl.xPoints, xARVc, 'k', linewidth=1.5)
    # plot Conduction Band L-Valley/X-Valley, Light Hole Valence Band and
    # Spin-Orbit coupling Valence Band
    for flag, xv, conf in (
                (plotVL, qcl.xVL, 'g--'), (plotVX, qcl.xVX, 'm-.'),
                (plotLH, qcl.xVLH, 'k--'), (plotSO, qcl.xVSO, 'r--')):
        if flag:
            axes.plot(qcl.xPoints, xv, conf, linewidth=1)
        ys.append(xv)
    return ys


def scaleWF(qcl: QCLayers, plotType: str = 'mode'):
    r"""Helper function to scale the wave function for plotting. The scale
    factor is decided by :data:`config`.

    Parameters
    ----------
    qcl :
        The QCLayers to plot
    plotType :
        Can be 'mode' or 'wf', to determine it's plotting mode
        (:math:`\psi^2`) or wavefunction itself (:math:`\psi`).

    Returns
    -------
    The scaled

    """
    if plotType == "mode":
        return qcl.psis**2 * config["mode_scale"]
    elif plotType == "wf":
        return qcl.psis * config["wf_scale"]
    else:
        raise ValueError('Undefined wavefunction time ', plotType)


def plotWF(axes: Axes, qcl: QCLayers, plotType: str = 'mode',
           fillPlot: Union[bool, float] = False,
           pickedStates: Iterable = set()):
    r"""Plot the wavefunctions of qcl.
    The wavefunctions are scaled by :func:`scaleWF`.

    Parameters
    ----------
    axes :
        The axes to plot the figure on.
    qcl :
        The QCLayers to plot
    plotType :
        Can be 'mode' or 'wf', to determine it's plotting mode
        (:math:`\psi^2`) or wavefunction itself (:math:`\psi`).
    fillPlot :
        Wether to fill up the wavefunctions. If it's `False` or `None` it will
        not fill, otherwise it should be a float number smaller than 1, meaning
        the transparency of the fill color.
    pickedStates :
        A set of state indices that should be plotted in thick black color.

    Returns
    -------
    A list of plotted data

    """
    colors = config['wf_colors']
    wfs = scaleWF(qcl, plotType)
    # filter almost zero part
    starts = np.argmax(abs(wfs) > config["wf_almost_zero"], axis=1)
    ends = np.argmax(abs(wfs[:, ::-1]) > config["wf_almost_zero"], axis=1)
    qcl.period_recognize()
    if qcl.status == 'solved-full':
        qcl.period_map_build()
        vmin = 0
        vmax = np.ceil(np.max(qcl.population)*10)/10
        popMap = cm.ScalarMappable(
            cmNorm(vmin=vmin, vmax=vmax), 'plasma')
    for n in range(len(qcl.eigenEs)):
        ls = '-'
        if n in pickedStates:
            color = 'k'
            lw = 2
        else:
            if qcl.status == 'solved-full':
                if qcl.periodMap[n] is not None:
                    color = popMap.to_rgba(
                        qcl.state_population(n))
                else:
                    color = 'g'
            else:
                color = colors[n % len(colors)]
            if qcl.status == 'basis':
                lw = 1.5
            else:
                if n in qcl.periodIdx:
                    # lw = 1 if n in qcl.unBound else 1.5
                    lw = 1.5
                    if n in qcl.unBound:
                        ls = (0, (0.5, 0.5))
                else:
                    lw = 0.5
        x = qcl.xPoints[starts[n]:-ends[n]]
        y = wfs[n, starts[n]:-ends[n]] + qcl.eigenEs[n]
        axes.plot(x, y, lw=lw, ls=ls, color=color)
        if fillPlot:
            axes.fill_between(x, y, qcl.eigenEs[n],
                              facecolor=color, alpha=fillPlot)
    if qcl.status == 'solved-full':
        colorbar_axes = axes.inset_axes([0.02, 0.01, 0.5, 0.02])
        axes.figure.colorbar(
            popMap, cax=colorbar_axes, orientation='horizontal',
            label='electron population')
        colorbar_axes.xaxis.set_ticks_position('top')
        colorbar_axes.xaxis.set_label_position('top')
    return wfs
