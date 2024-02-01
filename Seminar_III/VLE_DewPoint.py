# requesting for the thermodynamic properties
# available in script "thermo_properties.py "
import thermo_properties as TP

tol = 1.0e-10
delta = 1.0e-5

def Dew_point_ideal(compounds, P, y):
    # input: system (compounds) of interest
    # input: Pressure of the system, P.  Remember: in a Dew Temperature calculation
    # the Pressure  is known and constant.
    # input: vapor composition , y
    # ----------------------------------------------------------#
    # Requesting access to the 'pure_vapor_temperature' function, which is available within
    # the 'thermo_properties' script."
    Tvap_mix = TP.pure_vapor_temperature(compounds, P)
    # Tvap_mix : returns the values of Tsat for the individual compounds
    # at the given temperature P
    # ----------------------------------------------------------#
    # first guess of dew temperature in K : ------------------#
    # ------------------------------------------------------------
    To = sum(y * Tvap_mix)
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
        Pvap_mix = TP.pure_vapor_pressure(compounds, To)
        Ki = Pvap_mix / P  # K-values ( Raoult's model only !!)
        fequil = 1.0 - sum(y / Ki)  # phase-equilibrium criteria
        xo = y / Ki
        # ------------------------------------------------------------#
        # If the phase-equilibrium criteria is not satisfied a new pressure
        # is estimated using the Newton-Raphson (NR) method.
        # To apply the NR method, the derivative of the function (phase-equilibrium criteria)
        # Here the fequil is an implict function of the temperaure , therefore
        # the derivative is calculated numerically:
        Tup = To + delta
        Pvap_mixup = TP.pure_vapor_pressure(compounds, Tup)
        Kiup = Pvap_mixup / P
        fequil_up = 1.0 - sum(y / Kiup)
        def_fob = (fequil_up - fequil) / delta
        # ---------------------------------------------------------------#
        # NR formula for a new estimate of the dew temperature, T1
        T1 = To - (fequil / def_fob)
        To = T1

    Td = To
    return (Td, xo)
