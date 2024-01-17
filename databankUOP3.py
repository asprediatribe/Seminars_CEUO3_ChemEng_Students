# This script (databankUOP3) is a data bank with all the relevant constants to perform calculations of thermodynamic
# properties. All the values for the constants are taken from the appendix A
# of "Properties of Gases and Liquids". The outputs are defined as follows:
# 1) "antoine" correspond to the constant values of A, B, and C for the calculation of the vapor pressure,
# using the expression # 1 (antoine equation).
# 2) "critical" correspond to the values of critical temperature (K), pressure (bar),
# and volume (cm^3/mol) respectively
# 3) "Mw" molecular weight
# 4) "wac" acentric factor
#def thermodynamic_constants(compounds):
properties = {
        "compound1"  : {"antoine" : [, , ], "critical": [, , ], "Mw": , "wac": },
        "compound2"  : {"antoine" : [, , ], "critical": [, , ], "Mw": , "wac": }
            }


wace        = properties["compound"]["wac"]
cts_antoine = properties["compound"]["antoine"][1]
critics     = properties["compound"]["critical"][1]

print(wace , cts_antoine, critics)
#return(cts_antoine,critics,wace)
