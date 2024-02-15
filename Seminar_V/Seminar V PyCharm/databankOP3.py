import numpy as np
# ---------------------------------------------------------------------------------------------------------
# NumPy is the fundamental package for scientific computing in Python.
# It is a Python library that provides a multidimensional array object, (arrays, vectors, and matrices)
# and  routines for fast operations on arrays, including mathematical, logical,
# basic linear algebra, basic statistical operations, random simulation and much more.
# ----------------------------------------------------------------------------------------------------------
# As this database is general and works for any number of compounds the variables are stored in MATRICES
# (i.e.,"cts_antoine"). For exmple for the antoine constants , the OUTPUT is a matrix where the number
# of rows correspond to the number of compounds and the columns to the antoine constants.
# ---------------------------------------------------------------------------------------------------------
# This script is a data bank with all the relevant constants to perform calculations of thermodynamic
# properties. 
# This is the script with the oficial database of the course UOP3
# All the values for the constants are taken from the appendix A
# of "Properties of Gases and Liquids". The outputs are defined as follows:
# 1) "antoine" correspond to the constant values of A, B, and C for the calculation of the vapor pressure,
# using the expression # 1 (antoine equation).
# 2) "critical" correspond to the values of critical temperature (K), pressure (bar),
# and volume (cm^3/mol) respectively
# 3) "Mw" molecular weight
# 4) "wac" acentric factor
#-----------------------------------------------------------------------------------------------------------
# the function thermodynamic_constants gives as output the thermodynamic constants and properties when
# a compound or MIXTURE of compounds are specified
#-------------------------------------------------------------------------------------------------------------
def thermodynamic_constants(compounds):
    properties = {
        
        "pentane": {"antoine": [3.97786, 1064.840, 232.014], "critical": [469.70, 33.70, 311.00], "Mw": 72.150, "wac": 0.252},
        "heptane": {"antoine": [4.02023, 1263.909, 216.432], "critical": [540.20, 27.40, 428.00], "Mw":100.204, "wac": 0.350},
        "hexane": {"antoine": [4.00139, 1170.875, 224.317], "critical": [507.60, 30.25, 368.00], "Mw":86.177, "wac":0.300},
        "i-butane": {"antoine": [4.00272, 947.5400, 248.870], "critical": [407.85, 36.40, 262.70], "Mw": 58.123,
                     "wac": 0.186},
        "n-butane": {"antoine": [3.93266, 935.7730,238.789], "critical": [407.85, 36.40, 262.70], "Mw": 58.123,
                     "wac": 0.186},
        "i-pentane": {"antoine": [3.92023, 1022.880, 233.460], "critical": [460.39, 33.81, 308.30], "Mw": 72.150,
                      "wac": 0.229},
        "n-pentane": {"antoine": [3.97786,1064.840,232.014], "critical": [460.39, 33.81, 308.30], "Mw": 72.150,
                      "wac": 0.229},
        "butane": {"antoine": [3.93266, 935.7730, 238.789], "critical": [425.12, 37.96, 255.00], "Mw": 58.123,
                   "wac": 0.200},
        "octane": {"antoine": [4.05075, 1356.360, 209.635], "critical": [425.12, 37.96, 255.00], "Mw": 58.123,
                   "wac": 0.200},
        "nonane": {"antoine": [ 4.07356, 1438.030, 202.694], "critical": [425.12, 37.96, 255.00], "Mw": 58.123,
                   "wac": 0.200}

    }
    cts_antoine  = np.zeros((len(compounds),3)).tolist()
    # cts_antoine is a matrix of Zeros (0), containing the required compounds in the ROWS and the antoine constants,
    # for each compound, in COLUMNS
    
    for i in range(len(compounds)):
        
        for k in range(3):
            cts_antoine[i][k] = properties[compounds[i]]["antoine"][k]
            

    return(cts_antoine)
