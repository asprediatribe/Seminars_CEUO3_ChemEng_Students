def pure_vapor_temperature(compounds, P):
    # P = is the pressure of the system in bar
    # compound = is an array with the "names" of the compounds in the system
    # "names" must match with the names stored in the data_bank
    # saturation temperature ,t_sat, is given in K
    # t_sat is calculated using Equation #1 from Appendix A, solving for T
    cts = db.thermodynamic_constants(compounds)
    t_sat = np.zeros(len(compounds))
    for i in range(len(compounds)):
        A = cts[i][0]
        B = cts[i][1]
        C = cts[i][2]
        t_sat[i] = (B / (A - np.log10(P))) - (C - 273.15)
    return t_sat
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

def Bubble_temperature_ideal(compounds, P, x):
    # input: system (compounds) of interest
    # input: Pressure of the system, P.  Remember: in a Dew Temperature calculation
    # The pressure  is known and constant.
    # input: vapor composition, y
    # ----------------------------------------------------------#
    # Requesting access to the 'pure_vapor_temperature' function, which is available within
    # the 'thermo_properties' script."
    Tvap_mix = TP.pure_vapor_temperature(compounds, P)
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
        Pvap_mix = TP.pure_vapor_pressure(compounds, To)
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
        Pvap_mixup = TP.pure_vapor_pressure(compounds, Tup)
        Kiup = Pvap_mixup / P
        fequil_up = 1.0 -  sum(Kiup * x)
        def_fob = (fequil_up - fequil) / delta
        # ---------------------------------------------------------------#
        # NR formula for a new estimate of the dew temperature, T1
        T1 = To - (fequil / def_fob)
        To = T1

    Tb = To
    return (Tb, y)
