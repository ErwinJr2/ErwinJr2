#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ===========================================================================
# Reference
# [0]Vurgaftman, I., J. áR Meyer, and L. áR Ram-Mohan. "Band parameters for
#    III–V compound semiconductors and their alloys." Journal of applied 
#    physics 89.11 (2001): 5815-5875
# [1]Handbook of Optics, Vol.2, ISBN: 0070479747
# [2]Herve P J L, Vandamme L K J. Empirical temperature dependence of the
#    refractive index of semiconductors[J]. Journal of Applied Physics,
#    1995, 77(10): 5476-5477.
# [3]Joachim Pipre, Semiconductor Optoelectronic Devices: Introduction to
#    Physics and Simulation, SIBN: 0080469787
# [4]Varshni, Y. P. (1967). Temperature dependence of the energy gap in
#    semiconductors. physica, 34(1), 149-154
# ===========================================================================

from warnings import warn
# TODO: add warning for problematic datas
from numpy import sqrt

# flags for what effect to include
BOWING = True
VARSH  = True

# for In0.53Ga0.47As, EcG = 0.22004154
#    use this as a zero point baseline
bandBaseln = 0.22004154

class Material(object):
    """A semiconductor material class that stores material parameters
    
    Parameters
    ----------
    Name : str
        Name of material
    Temperature : int
        Temperature of the material

    """
    def __init__(self, Name, Temperature=300):
        """
        
        Yields
        ------
        type : str
            Type of the material
        
        """
        self.name = Name
        self.parm = MParm[self.name].copy()
        self.type = self.parm.pop("Crystal")
        self.set_temperature(Temperature)

    def set_temperature(self, Temperature):
        """ 
        Set Temperature of the material and update related parameters: 
        lattice consant and band gap.
        
        Yields
        ------
        T : int
            Updated temperature
        parm : dict
            lattice constant and band gap in this dictionary are updated

        """
        self.T = Temperature
        for k in self.parm:
            if k.endswith("lc") and k+"_T" in self.parm:
                # k is about lattice constant
                self.parm[k] = (MParm[self.name][k] + (self.T - 300)
                                * MParm[self.name][k+"_T"])
        # Varshni correction to bandgap [4]
        # major assumption:
        #   Correction only to conduction band (no influence on VBO)
        #   Although in reality the Varshni correction should be part
        #   conduction band, part valence
        if VARSH:
            for pt in ('G', 'X', 'L'): 
                Varsh = -self.parm["al"+pt] * self.T**2 / (
                    self.T + self.parm["be"+pt])
                self.parm['Eg'+pt] = MParm[self.name]['Eg'+pt] + Varsh

    def set_strain(self, a_parallel):
        """
        Update parameters' dependence on strain, according to Pikus-Bir 
        interaction.

        Parameters
        ----------
        a_parallel : float
            lattice constant of the substrate
        

        Yields
        ------
        eps_parallel : float
            Strain tensor within/parallel to the layer plane
        a_perp : float
            Lattice const. perpendicular to the layer plane
        eps_perp : float
            Strain tensor perpendicular to the layer plane
        parm : dict
            Update parameters' dependence on strain
            
        """
        # TODO: crystal with more than one lattice constant? 
        # eps_parallel: strain tensor within/parallel to the layer plane
        self.eps_parallel = a_parallel / self.parm["alc"] - 1
        # self.a_perp: lattice const. perpendicular to the layer plane
        # Ask if MBE growers care about strain on monolayer thickness
        self.a_perp = self.parm["alc"] * (1 - 2 * self.parm["c12"] /
                                     self.parm["c11"] * self.eps_parallel)
        # eps_perp: strain tensor perpendicular to the layer plane
        self.eps_perp = self.a_perp / self.parm["alc"] - 1
        if self.type == "ZincBlende":
            Pec = (2 * self.eps_parallel + self.eps_perp) * self.parm["acG"]
            Pe  = (2 * self.eps_parallel + self.eps_perp) * self.parm["av"]
            Qe  = (- self.parm["b"] * (self.parm["c11"] + 2 
                                       * self.parm["c12"]) 
                   / self.parm["c11"] * self.eps_parallel)
            self.parm["EcG"] = (self.parm["VBO"] + self.parm["EgG"] + Pec 
                                - bandBaseln)
            self.parm["EcL"] = (self.parm["VBO"] + self.parm["EgL"] 
                                + (2 * self.eps_parallel + self.eps_perp) 
                                * (self.parm["acL"] + self.parm["av"]) 
                                - bandBaseln)
            self.parm["EcX"] = (self.parm["VBO"] + self.parm["EgX"] 
                                + (2 * self.eps_parallel + self.eps_perp) 
                                * (self.parm["acX"] + self.parm["av"]) 
                                + 2/3*self.parm["XiX"]*(
                                    self.eps_perp - self.eps_parallel)
                                - bandBaseln)
            self.parm["ESO"] = sqrt(9 * Qe**2 + 2 * Qe * self.parm["DSO"] +
                                    self.parm["DSO"]**2)
            self.parm["EgLH"] = (self.parm["EgG"] + Pec + Pe - 1/2 * (
                Qe - self.parm["DSO"] + self.parm["ESO"]))
            self.parm["EgSO"] = (self.parm["EgG"] + Pec + Pe - 1/2 * (
                Qe - self.parm["DSO"] - self.parm["ESO"]))
            self.parm["EvLH"] = self.parm["EcG"] - self.parm["EgLH"] 
            self.parm["EvSO"] = self.parm["EcG"] - self.parm["EgSO"]
        else:
            self.parm["EcG"] = (self.parm["VBO"] + self.parm["EgG"]
                                - bandBaseln)

class Alloy(Material):
    """
    An alloy of material with mole fraction $x$ 

    Parameters
    ----------
    Name : str
        Name of the alloy
    x : float
        The first composition's Mole fraction defined in AParm
    Temperature : int
        Temperature of the alloy

    """
    def __init__(self, Name, x, Temperature=300):
        self.name = Name
        self.comp = AParm[self.name]["composition"]
        self.A = Material(self.comp[0], Temperature)
        self.B = Material(self.comp[1], Temperature)
        assert(self.A.type == self.B.type)
        self.type = self.A.type
        # so alloy is A_x B_(1-x)
        self.set_temperature(Temperature)
        self.set_molefrac(x)

    def set_temperature(self, Temperature):
        """ 
        Set Temperature of the alloy and update related parameters, 
        lattice consant and band gap, by updating the temperature of the
        materials in the alloy


        Yields
        ------
        T : int
            Updated temperature
        """
        self.T = Temperature
        self.A.set_temperature(self.T)
        self.B.set_temperature(self.T)

    def set_molefrac(self, x):
        """
        Update parameters of the alloy with Mole fraction x

        Parameters
        ----------
        x : float
            The first composition's Mole fraction defined in AParm

        Yields
        ------
        parm : dict
            stores the parameter of the alloy
        """
        self.parm = {}

        # For the Gamma band gap bowing in AlxGa1-xAs and AlxGa1-xSb
        # a linear or even a quadratic interpolation is not sufficient.
        # Here, the bowing parameter is not constant,
        # it depends on the alloy composition x.
        if self.name == 'AlGaAs': 
            AParm[self.name]['EgG'] = -0.127 + 1.310 * x
        if self.name == 'AlGaSb': 
            AParm[self.name]['EgG'] = -0.044 + 1.22 * x

        for k in self.A.parm:
            #  if not k in B.parm:
            #      continue
            self.parm[k] = (x * self.A.parm[k] + (1-x) * self.B.parm[k])
            if BOWING and k in AParm[self.name]:
                # Bowing parameters
                self.parm[k] -= x * (1-x) * AParm[self.name][k]

MParm = {
    "GaAs":{  # from Vurgaftman[0] unless specified
        "Crystal": "ZincBlende", 
        # Lattice constant and thermal expension
        "alc": 5.65325,       # Angstrom, Ref[3] Table 2.6
        "alc_T": 3.88e-5,     # Angstrom / K
        # Elestic stiffness constants
        'c11': 1221,   # GPa
        'c12': 566,    # GPa
        'EgG': 1.519,  # eV, band gap at Gamma at 0 K
        'EgL': 1.815,  # eV, band gap at L
        'EgX': 1.981,  # eV, band gap at X
        'VBO': -0.80,  # eV, valence band offset
        'DSO': 0.341,  # eV, Delta_SO from top of valence band
        # Pikus-Bir interaction parameters: strain correction
        'acG': -7.17,  # eV
        'acL': -4.91,  # eV, NextNano DB
        'acX': -0.16,  # eV, NaxtNano DB
        'av': -1.16,   # eV, valence band deformation potential
        'b': -2.0,     # eV
        'XiG': 0,      # eV
        'XiL': 6.5,    # eV, NaxtNano DB
        'XiX': 14.26,  # eV, NaxtNano DB
        'me0': 0.067,  # m0, effective mass at Gamma
        # Kane model parameters for effective mass
        'Ep': 28.8,    # eV
        'F': -1.94,    # unit 1
        # Varshni correction
        'alG': 0.5405e-3,   # eV/K, Varshni alpha (at Gamma point)
        'beG': 204,         # K, Varshni beta (at Gamma point)
        'alX': 0.460e-3,    # K, Varshni alpha at X
        'beX': 204,         # K, Varshni beta at X
        'alL': 0.605e-3,    # K, Varshni alpha at L
        'beL': 204,         # K, Varshni beta at L
        'epss': 12.9,       # static permitivity
        'epsInf': 10.86,    # high-freq permitivity
        'hwLO': 35.3e-3   #  LO phonon energy eV
    },

    'InAs':{ # from Vurgaftman[0] unless specified
        "Crystal": "ZincBlende", 
        "alc": 6.0583, "alc_T": 2.74e-5,
        'c11': 832.9, 'c12': 452.6,
        'EgG': 0.417, 'EgL': 1.133, 'EgX': 1.433, 'VBO': -0.59, 'DSO': 0.39,
        'acG': -5.08, 'acL': -3.89, # eV, NextNano DB
        'acX': -0.08,  # eV, NextNano DB
        'av': -1.00, 'b': -1.8, 'XiG': 0, 'XiL': 11.35, 'XiX': 3.7,
        'me0': 0.026, 'Ep': 21.5, 'F': -2.9,
        'alG': 0.276e-3, 'beG': 93,
        'alX': 0.276e-3, 'beX': 93,
        'alL': 0.276e-3, 'beL': 93,
        'epss': 14.3, 'epsInf': 11.6,
        'hwLO': 29.93e-3
    },

    'AlAs':{ # from Vurgaftman[0] unless specified
        "Crystal": "ZincBlende", 
        "alc": 5.6611, "alc_T": 2.90e-5,
        'c11': 1250, 'c12': 534,
        'EgG': 3.099, 'EgL': 2.46, 'EgX': 2.24, 'VBO': -1.33, 'DSO': 0.28,
        'acG': -5.64, 'acL': -3.07,  # NextNano DB
        'acX': 2.54,   # NextNano DB
        'av': -2.47, 'b': -2.3, 'XiG': 0, 'XiL': 11.35, 'XiX': 6.11,
        'me0': 0.15, 'Ep': 21.1, 'F': -0.48,
        'alG': 0.855e-3, 'beG': 530,
        'alX': 0.70e-3, 'beX': 530,
        'alL': 0.605e-3, 'beL': 204,
        'epss': 10.06, 'epsInf': 8.16,
        'hwLO': 49.8e-3
    },

    'AlSb':{ # from Vurgaftman[0] unless specified
        "Crystal": "ZincBlende", 
        "alc": 6.1355, "alc_T": 2.60e-5,
        'c11': 876.9, 'c12': 434.1,
        'EgG': 2.386, 'EgL': 2.329, 'EgX': 1.696, 'VBO': -0.41, 'DSO': 0.676,
        'acG': -4.5, 'acL': 0,     # NextNano DB
        'acX': 2.54,  # NextNano DB
        'av': -1.4, 'b': -1.35, 'XiG': 0, 'XiL': 11.35, 'XiX': 6.11,
        'me0': 0.14, 'Ep': 18.7, 'F': -0.56,
        'alG': 0.42e-3, 'beG': 140,
        'alX': 0.39e-3, 'beX': 140,
        'alL': 0.58e-3, 'beL': 140,
        'epss': 12.04,    # ISBN 0849389127
        'epsInf': 10.24,  # ISBN 0849389127
        'hwLO': 42.7      # http://prb.aps.org/pdf/PRB/v43/i9/p7231_1
    },

    'GaSb':{ # from Vurgaftman[0] unless specified
        "Crystal": "ZincBlende", 
        "alc": 6.0959, "alc_T": 4.72e-5,
        'c11': 884.2, 'c12': 402.6,
        'EgG': 0.812, 'EgL': 0.875, 'EgX': 1.141, 'VBO': -0.03, 'DSO': 0.76,
        'acG': -7.5, 
        'av': -0.8, 'b': -2.0,
        'me0': 0.039, 'Ep': 27.0, 'F': -1.63,
        'alG': 0.417e-3, 'beG': 140,
        'alX': 0.475e-3, 'beX': 94,
        'alL': 0.597e-3, 'beL': 140,
    },

    'InSb':{ # from Vurgaftman[0] unless specified
        "Crystal": "ZincBlende", 
        "alc": 6.4794, "alc_T": 3.48e-5,
        'c11': 684.7, 'c12': 373.5,
        'EgG': 0.235, 'EgL': 0.93, 'EgX': 0.63, 'VBO': 0, 'DSO': 0.81,
        'acG': -6.94,
        'av': -0.36, 'b': -2.0,
        'me0': 0.0135, 'Ep': 23.3, 'F': -0.23,
        'alG': 0.32e-3, 'beG': 170,
        # alX, beX, alL, beL are unknown... using gamma point data instead
        'alX': 0.32e-3, 'beX': 170,
        'alL': 0.32e-3, 'beL': 170,
    },

    'InP':{
        "Crystal": "ZincBlende", 
        "alc": 5.869, "alc_T": 2.79e-5,
        'c11': 1011, 'c12': 561, 
        'EgG': 1.4236, 'EgL': 2.014, 'EgX': 2.384, #-3.7e-4 * T
        'VBO': -0.94, 'DSO': 0.108, 
        'acG': -6.0, 'av': -0.6, 'b': -2.0, 
        'alG': 0.363, 'beG': 162, 
        'alL': 0.363, 'beL': 162, 
        # alX, beX are unknown... using gamma point data instead
        'alX': 0.363, 'beX': 162, 
        'me0': 0.0795, 'Ep': 20.7, 'F':-1.31,
    },
    
    "alpha-GaN":{
        "Crystal": "Wurtzite"
    },
    "beta-GaN":{
        "Crystal": "ZincBlende",
        # alc_T not given, use 0
        "alc": 4.50, "alc_T": 0,
        'c11': 293, 'c12': 159,
        'EgG': 3.299, 'EgL': 5.59, 'EgX': 4.52,
        'VBO': -2.64, 'DSO': 0.017,
        'acG': -2.2, 'av': -5.2, 'b': -2.2,
        'alG': 0.593e-3, 'beG': 600,
        'alL': 0.593e-3, 'beL': 600,
        'alX': 0.593e-3, 'beX': 600,
        'me0': 0.15, 'Ep': 25.0, 'F': -0.92,
    },

    "alpha-AlN":{
        "Crystal": "Wurtzite"
    },
    "beta-AlN":{
        "Crystal": "ZincBlende",
        # alc_T not given, use 0                              
        "alc": 4.38, "alc_T": 0,
        'c11': 304, 'c12': 160,
        'EgG': 4.9, 'EgL': 9.3, 'EgX': 6.0,
        'VBO': -3.44, 'DSO': 0.019,
        'acG': -6.0, 'av': -3.4, 'b': -1.9,
        'alG': 0.593e-3, 'beG': 600,
        'alL': 0.593e-3, 'beL': 600,
        'alX': 0.593e-3, 'beX': 600,
        'me0': 0.25, 'Ep': 27.1, 'F': 0.76,
    },

    "alpha-InN":{
        "Crystal": "Wurtzite"
    },
    "beta-InN":{
        "Crystal": "ZincBlende",
        # alc_T not given, use 0
        "alc": 4.98, "alc_T": 0,
        'c11': 187, 'c12': 125,
        'EgG': 1.94, 'EgL': 5.82, 'EgX': 2.51,
        'VBO': -2.38, 'DSO': 0.006,
        'acG': -1.85, 'av': -1.5, 'b': -1.2,
        'alG': 0.245e-3, 'beG': 624,
        'alL': 0.245e-3, 'beL': 624,
        'alX': 0.245e-3, 'beX': 624,
        'me0': 0.12, 'Ep': 25.0, 'F': -0.92,
    }


}

AParm = {
    'InGaAs':{
        'EgG': 0.477, 'EgL': 0.33, 'EgX': 1.4, 'VBO': -0.38, 'DSO': 0.15,
        'acG': 2.61, 'acL': 2.61,  # NextNano DB
        'acX': 2.61,  # NextNano DB
        'me0': 0.0091, 'Ep': -1.48, 'F': 1.77,
        'name': 'InxGa1-xAs',
        'composition': ('InAs', 'GaAs')
    }, 

    'AlInAs':{
        'EgG': 0.70, 'EgL': 0, 'EgX': 0, 'VBO': -0.64, 'DSO': 0.15,
        'acG': -1.4, 'acL': -1.4,  # NextNano DB
        'acX': -1.4,  # NextNano DB
        'me0': 0.049, 'Ep': -4.81, 'F': -4.44,
        'name': 'Al1-xInxAs',
        'composition': ('InAs', 'AlAs')
    }, 

    'AlGaAs':{
        'EgG': -0.127,  # + 1.310*x(Al)
        # To describe the band gap bowing at the Gamma point in AlxGa1-xAs,
        # a linear or even a quadratic interpolation is not sufficient.
        # Here, the bowing parameter is not constant,
        # it depends on the alloy composition x.
        'EgL': 0.055, 'EgX': 0, 'VBO': 0, 'DSO': 0,
        'acG': 0, 'acL': 0, 'acX': 0,
        'me0': 0, 'Ep': 0, 'F': 0,
        'name': 'AlxGa1-xAs',
        'composition': ('AlAs', 'GaAs')
    },

    'AlAsSb':{
        'EgG': 0.8, 'EgL': 0.28, 'EgX': 0.28, 'DSO': 0.15, 'VBO': -1.71,
        'name': 'AlAsxSb1-x',
        'composition': ('AlAs', 'AlSb')
    }, 

    'AlGaSb':{
        'EgG': -0.044,  # + 1.22*x(Al), same as AlGaAs
        'EgL': 0, 'EgX': 0, 'VBO': 0, 'DSO': 0.2,
        'acG': 0, 'acL': 0,  # NextNano DB
        'acX': 0,  # NextNano DB
        'me0': 0, 'Ep': 0, 'F': 0,
        'name': 'AlxGa1-xSb',
        'composition': ('AlSb', 'GaSb')
    },

    'InAsSb':{
        'EgG': 0.67, 'EgL': 0.6, 'EgX': 0.6, 'VBO': 0, 'DSO': 1.2,
        'acG': 0, 'acL': 0,  # NextNano DB
        'acX': 0,  # NextNano DB
        'me0': 0.035, 'Ep': 0, 'F': 0,
        'name': 'InAsxSb1-x',
        'composition': ('InAs', 'InSb')
    }
}

def main(material):
    print("Looking for parameters of %s:" % material)
    if material in AParm:
        A, B = AParm[material]['composition']
        print("Alloy with (%s)x(%s)1-x"%(A,B))
        print("%s: "%A, MParm[A])
        print("%s: "%B, MParm[B])
        print("Bowing parameters:", AParm[material])
        return 0
    elif material in MParm: 
        print(MParm[material])
        return 0
    else:
        print("Not found.")
        return 1

if __name__ == "__main__":
    import sys
    for material in sys.argv[1:]: 
        main(material)
        print("")

#m vim: ts=4 sw=4 sts=4 expandtab
