# requesting for the thermodynamic properties
# available in script "thermo_properties.py "
import thermo_properties as TP
# Remember: 'TP' is simply a shorthand way to reference the script "thermo_properties"
# --------------------------------------------------------#
# input: system (compounds) of interest
compounds = ["pentane", "heptane"]
# input: Temperature of the system.  Remember: in a Bubble Pressure calculation
# the temperatute is known and constant.
T = 400 # K
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
# first guess of bubble pressure in bar
Po = 5.4067 # bar
# -----------------------------------------------------------#
Ki = Pvap_mix / Po  # K-values
y  = Ki * x         # estimation of the composition of vapour #
fequil = 1.0 - sum(Ki * x)  # phase-equilibrium criteria
# ------------------------------------------------------------#
# If the phase-equilibrium criteria is not satisfied a new pressure
# is estimated using the Newton-Raphson (NR) method.
# To apply the NR method, the derivative of the function (phase-equilibrium criteria)
defun = sum(Ki * x) / Po
# ---------------------------------------------------------------#
# NR formula for a new estimate of the bubble pressure, P1
P1 = Po - (fequil / defun)
print(P1)
print(fequil)
