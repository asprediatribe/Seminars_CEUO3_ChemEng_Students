# This script (databankUOP3) is a data bank with all the relevant constants to perform calculations of thermodynamic
# properties. All the values for the constants are taken from the appendix A
# of "Properties of Gases and Liquids". The outputs are defined as follows:
# 1) "antoine" correspond to the constant values of A, B, and C for the calculation of the vapor pressure,
# using the expression # 1 (antoine equation).
# 2) "critical" correspond to the values of critical temperature (K), pressure (bar),
# and volume (cm^3/mol) respectively
# 3) "Mw" molecular weight
# 4) "wac" acentric factor
# properties is the name of the dictionary
properties = {
    "benzene": {"antoine": [3.98523, 1184.240, 217.572], "critical": [562.05, 48.95, 256.00], "Mw": 78.114,"wac": 0.210},
    "water"  : {"antoine": [5.11564, 1687.537, 230.17],  "critical": [647.14, 220.64, 55.95], "Mw": 18.015, "wac": 0.344},
    "pentane": {"antoine": [3.97786, 1064.840, 232.014], "critical": [469.70, 33.70, 311.00], "Mw": 72.150, "wac": 0.252},
    "heptane": {"antoine": [4.02023, 1263.909, 216.432], "critical": [540.20, 27.40, 428.00], "Mw":100.204, "wac": 0.350}
}

# create a python list "cts_antoine" to store the values of the  (3) Antoine's constants
cts_antoine = list(range(3))
for i in range(3):
    cts_antoine[i] = properties["Compound"]["antoine"][i]
    # print(cts_antoine)
print(cts_antoine)
A= cts_antoine[0]
B= cts_antoine[1]
C= cts_antoine[2]

T = # Temparature in K
logP = A - (B/(T+C-273.15))
Psat = pow(10, logP)
print(" sat pressure in bar", Psat)

