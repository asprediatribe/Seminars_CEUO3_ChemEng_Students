# requesting for the thermodynamic properties
# available in script "thermo_properties.py "
import thermo_properties as TP

# Remember: 'TP' is simply a shorthand way to reference the script "thermo_properties"
# --------------------------------------------------------#
tol = 1.0e-10
# input: system (compounds) of interest
compounds = ["pentane", "heptane"]
# input: Temperature of the system.  Remember: in a Bubble Pressure calculation
# the temperatute is known and constant.
T = 400  # K
# ----------------------------------------------------------#
# Requesting access to the 'pure_vapor_pressure' function, which is available within
# the 'thermo_properties' script."
Pvap_mix = TP.pure_vapor_pressure(compounds, T)
# Pvap_mix : returns the values of Psat for the individual compounds
# at the given temperature T
# ----------------------------------------------------------#
# input: liquid composition of the mixture
x = [0.4, 0.6]
# ----------------------------------------------------------#
# first guess of bubble pressure in bar : ------------------#
# This is just two different ways to specify a more
# realistic initial guess :
# 1 - Po is apx. the mean of the Psat of the pure compounds
# Po = sum(Pvap_mix) / 2
# ------------------------------------------------------------
# 2 - Po is appx. to the pressure  given by Partial Pressure Law in ideal mixtures
 Po = sum(x*Pvap_mix)
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
    defun = sum(Ki * x) / Po  # ( Raoult's model only !!)
    # ---------------------------------------------------------------#
    # NR formula for a new estimate of the bubble pressure, P1
    P1 = Po - (fequil / defun)
    Po = P1
    # print(P1)
    # print(fequil)
print("Bubble pressure (bar) of system", compounds, "at T =", T, "K is:", P1)
print("the composition of vapor coexisting at equilibria is: ", y)
