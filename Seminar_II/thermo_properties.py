import numpy as np
import databank_OP3 as db


def pure_vapor_pressure(compounds, T):
    # T = is the absolute temperature of the system in K
    # compound = is an array with the "names" of the compounds in the system
    # the "names" must match with the names stored in the data_bank
    # vapor pressure ,p_vap, is given in bar
    # vapor pressure is calculated using Equation #1 from Appendix A
    cts, critics, wace = db.thermodynamic_constants(compounds)
    p_vap = np.zeros(len(compounds))
    for i in range(len(compounds)):
        A = cts[i][0]
        B = cts[i][1]
        C = cts[i][2]
        exp_vap = A - (B / (T + C - 273.15))
        p_vap[i] = pow(10, exp_vap)
    return p_vap

