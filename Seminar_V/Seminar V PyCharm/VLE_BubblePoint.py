# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 16:40:40 2024



@author: s2216474
"""

import thermo_properties as TP

tol = 0.001
delta=0.0001

def VLE_Tbubble(x,P,compounds):
    # input: system (compounds) of interest
    # input: Pressure of the system, P.  Remember: in a Dew Temperature calculation
    # The pressure  is known and constant.
    # input: vapor composition, y
    # ----------------------------------------------------------#
    # Requesting access to the 'pure_vapor_temperature' function, which is available within
    # the 'thermo_properties' script."
    print(x,P,compounds)
    Tvap_mix = TP.pure_vapour_temperature(compounds, P)
    print(Tvap_mix)
    # Tvap_mix : returns the values of Tsat for the individual compounds
    # at the given temperature P
    # ----------------------------------------------------------#
    # first guess of dew temperature in K : ------------------#
    # ------------------------------------------------------------
    # 2 - Po is apx. the pressure given by Partial Pressure Law in ideal mixtures
    To = sum(x * Tvap_mix)
    # ----------------------------------------------------------------
    # Set a value (on purpose) for the equilibrium function.
    # This is necessary to allow the code to proceed with the instructions within the While loop.
    # notice that if  fequil = 1.0 , then  1.0 is > tol, therefore the code executes the While loop.
    fequil = 1.0
    # -----------------------------------------------------------#
    # Once inside the While loop, the instructions are executed over and over
    # until the condition is not longer satisfied, it is to say until : abs(fequil) <= tol, which means
    # until the system coexist at VLE.
    # -----------------------------------------------------------#
    while abs(fequil) > tol:
        Pvap_mix = TP.pure_vapour_pressure(compounds, To)
        Ki = Pvap_mix / P  # K-values ( Raoult's model only !!)
        fequil =  1.0 - sum(Ki * x)  # phase-equilibrium criteria
        y = Ki * x
        # ------------------------------------------------------------#
        # If the phase-equilibrium criteria is not satisfied a new pressure
        # is estimated using the Newton-Raphson (NR) method.
        # To apply the NR method, the derivative of the function (phase-equilibrium criteria)
        # Here the fequil is an implicit function of the temperature, therefore
        # The derivative is calculated numerically:
        Tup = To + delta
        Pvap_mixup = TP.pure_vapour_pressure(compounds, Tup)
        Kiup = Pvap_mixup / P
        fequil_up = 1.0 -  sum(Kiup * x)
        def_fob = (fequil_up - fequil) / delta
        # ---------------------------------------------------------------#
        # NR formula for a new estimate of the dew temperature, T1
        T1 = To - (fequil / def_fob)
        To = T1

    Tb = To
    return (Tb, y)




# --------------------------------------------------------#
tol = 1.0e-10
# input: system (compounds) of interest
# input: Temperature of the system.  Remember: in a Bubble Pressure calculation
# the temperatute is known and constant.
# ----------------------------------------------------------#
# Requesting access to the 'pure_vapor_pressure' function, which is available within
# the 'thermo_properties' script."
def VLE_Pbubble(x,T,compounds):

    Pvap_mix = TP.pure_vapour_pressure(compounds, T)
    # Pvap_mix : returns the values of Psat for the individual compounds
    # at the given temperature T
    # ----------------------------------------------------------#
    # input: liquid composition of the mixture
    # ----------------------------------------------------------#
    # first guess of bubble pressure in bar : ------------------#
    # This is just two different ways to specify a more
    # realistic initial guess :
    # 1 - Po is apx. the mean of the Psat of the pure compounds
    Po = sum(Pvap_mix) / 2
    # ------------------------------------------------------------
    # 2 - Po is appx. to the pressure  given by Partial Pressure Law in ideal mixtures
    # Po = sum(x*Pvap_mix)
    # ----------------------------------------------------------------
    # Set a value (on porpose) for the equilibrium function.
    # This is necessary to allow the code to proceed with the instructions within the While loop.
    # notice that if  fequil = 1.0 , then  1.0 is > tol, therefore the code executes the While loop.
    fequil = 1.0
    # -----------------------------------------------------------#
    # Once inside the While loop, the instructions are executed over and over
    # until the condition is not longer satisfied, it is to say until : abs(fequil) <= tol, which means
    # until the system coexist at VLE.
    # -----------------------------------------------------------#
    while abs(fequil) > tol:
        Ki = Pvap_mix / Po  # K-values ( Raoult's model only !!)
        y = Ki * x  # estimation of the composition of vapour #
        fequil = 1.0 - sum(Ki * x)  # phase-equilibrium criteria
        # ------------------------------------------------------------#
        # If the phase-equilibrium criteria is not satisfied a new pressure
        # is estimated using the Newton-Raphson (NR) method.
        # To apply the NR method, the derivative of the function (phase-equilibrium criteria)
        defun = 10*sum(Ki * x) / Po  # ( Raoult's model only !!)
        # ---------------------------------------------------------------#
        # NR formula for a new estimate of the bubble pressure, P1
        P1 = Po - (fequil / defun)
        Po = P1
        
        return Po,y

