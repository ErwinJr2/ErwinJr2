"""
This contains the plotter functions for QCLayers
"""

from typing import Iterable, Union

import numpy as np
from matplotlib import cm
from matplotlib.axes import Axes
from matplotlib.colors import Normalize as cmNorm
from matplotlib.pyplot import gca

from ErwinJr2.qc_layers import QCLayers

config = {
    "wf_scale": 0.3,
    "mode_scale": 3,
    # This number not necessarily but is chosen to be consistent with the
    # default value of tol for QCLayers.periodRecognize as
    # wf_almost_zero / modescale = tol
    "wf_almost_zero": 1.5e-4,
    "default_lw": 1.0,
    "wf_colors": (
        (0.584, 0.450, 0.701),
        (0.431, 0.486, 0.745),
        (0.576, 0.694, 0.517),
        (0.682, 0.780, 0.321),
        (0.501, 0.501, 0.509),
        (0.854, 0.741, 0.247),
        (0.874, 0.607, 0.290),
        (0.823, 0.341, 0.278),
        (0.725, 0.321, 0.623),
        (0.411, 0.741, 0.270),
        (0.078, 0.078, 0.078),
        (0.431, 0.803, 0.870),
        (0.223, 0.321, 0.643),
    ),
}  #: The configuration dictionary for plotting.


def plot_potential(
    qcl: QCLayers,
    axes: Axes = None,
    plot_vl: bool = False,
    plot_vx: bool = False,
    plot_lh: bool = False,
    plot_so: bool = False,
):
    """Plot the potentials of qcl.

    Parameters
    ----------
    qcl :
        The QCLayers to plot
    axes :
        The axes to plot the figure on.
    plotVL, plotVX, plotLH, plotSO :
        flags wether to plot L-point, X-point, LH band, SO band

    Returns
    -------
    A list of plotted data

    """
    if axes is None:
        axes = gca()
    ys = [qcl.x_vc]
    axes.plot(qcl.x_points, qcl.x_vc, "k", linewidth=config["default_lw"])
    # highlight selected layer & make AR layers bold
    x_non_ars = np.bitwise_and.reduce(
        [qcl.x_layer_mask(n) for n, ar in enumerate(qcl.layer_ar) if ar]
    )
    x_ar_vc = np.ma.masked_where(x_non_ars, qcl.x_vc)
    axes.plot(qcl.x_points, x_ar_vc, "k", linewidth=config["default_lw"] + 0.5)
    # plot Conduction Band L-Valley/X-Valley, Light Hole Valence Band and
    # Spin-Orbit coupling Valence Band
    for flag, xv, conf in (
        (plot_vl, qcl.x_vl, "g--"),
        (plot_vx, qcl.x_vx, "m-."),
        (plot_lh, qcl.x_vlh, "k--"),
        (plot_so, qcl.x_vso, "r--"),
    ):
        if flag:
            axes.plot(qcl.x_points, xv, conf, linewidth=config["default_lw"])
        ys.append(xv)
    return ys


def scale_wf(qcl: QCLayers, plot_type: str = "mode"):
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
    The scaled wavefunction.

    """
    if plot_type == "mode":
        return qcl.psis**2 * config["mode_scale"]
    elif plot_type == "wf":
        return qcl.psis * config["wf_scale"]
    else:
        raise ValueError("Undefined wavefunction time ", plot_type)


def plot_wf(
    qcl: QCLayers,
    plot_type: str = "mode",
    fill_plot: Union[bool, float] = False,
    picked_states: Iterable = tuple(),
    show_period: bool = True,
    axes: Axes = None,
):
    r"""Plot the wavefunctions of qcl.
    The wavefunctions are scaled by :func:`scaleWF`.

    Parameters
    ----------
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
    showPeriod :
        Flag to whether emphasis the recognized period of wavefunctions.
    axes :
        The axes to plot the figure on.

    Returns
    -------
    A list of plotted data

    """
    if axes is None:
        axes = gca()
    colors = config["wf_colors"]
    wfs = scale_wf(qcl, plot_type)
    # filter almost zero part
    starts = np.argmax(abs(wfs) > config["wf_almost_zero"], axis=1)
    ends = np.argmax(abs(wfs[:, ::-1]) > config["wf_almost_zero"], axis=1)
    show_pop = False
    if show_period:
        qcl.period_recognize()
        show_pop = qcl.status == "solved-full"
    if show_pop:
        qcl.period_map_build()
        vmin = 0
        vmax = np.ceil(np.max(qcl.population) * 10) / 10
        pop_map = cm.ScalarMappable(cmNorm(vmin=vmin, vmax=vmax), "plasma")
    for n in range(len(qcl.eigen_es)):
        ls = "-"
        if n in picked_states:
            color = "k"
            lw = config["default_lw"] * 2
        else:
            if show_pop:
                if qcl.period_map[n] is not None:
                    color = pop_map.to_rgba(qcl.state_population(n))
                else:
                    color = "g"
            else:
                color = colors[n % len(colors)]
            if qcl.status == "basis":
                lw = config["default_lw"]
            else:
                if show_period and n in qcl.period_idx:
                    # lw = 1 if n in qcl.unBound else 1.5
                    lw = config["default_lw"]
                    if n in qcl.un_bound:
                        ls = (0, (0.5, 0.5))
                else:
                    lw = config["default_lw"] / 2
        x = qcl.x_points[starts[n] : -ends[n]]
        y = wfs[n, starts[n] : -ends[n]] + qcl.eigen_es[n]
        axes.plot(x, y, lw=lw, ls=ls, color=color)
        if fill_plot:
            axes.fill_between(x, y, qcl.eigen_es[n], facecolor=color, alpha=fill_plot)
    if show_pop:
        colorbar_axes = axes.inset_axes([0.03, 0.01, 0.5, 0.02])
        axes.figure.colorbar(
            pop_map,
            cax=colorbar_axes,
            orientation="horizontal",
            label="electron population",
        )
        colorbar_axes.xaxis.set_ticks_position("top")
        colorbar_axes.xaxis.set_label_position("top")
    return wfs
