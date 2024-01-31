import thermo_properties as prop
# Python script to solve the Isothermal Flash problem (PTFlash)
# using the Rachford & Rice Algorithm. Please review class material LUO2a
# and self-study material , Chapter 3, (Week 3).
tol = 1.0e-10
# ----- specify the system -----------------------------------------------
compounds = []
# ------------------------------------------------------------------------
# -- Specify feed and flash conditions -----------------------------------
T_op =   # K                    Flash drum operating  temperature, T*
P_op =   # bar                  Flash drum operating  pressure, P*
z = []   # Feed compositions, Zf
F =      #   [kmol/h] Feed stream molar flow, F
# If F  is not given set a  calculation base
# -------------------------------------------------------------------------
# Initial guesses ---------------------------------------------------------
psi_o =  # initial guess for the phase fraction (psi)
xo = []
# Initial values for K-values --------------------------------------------------
# The calculation of K-values, in ideal systems, depends on
# the saturation pressure (Raoult's Law).
Pi_sat = prop.pure_vapor_pressure(compounds, T_op)
# ------------------------------------------------------------------------------
# As we need to executed the algorithm ( Algorithm # 2 in Self_study
# material) over and over until the R & R function be satisfied, we set a value
# of the function to start the sequence of calculations.
fob_rr = 1.0
while abs(fob_rr) > tol:
    # 1 - K-values implementing Raoult's law
    Ki = Pi_sat / P_op
    # ---------------------------------------------------------------------------------
    # Evalute the R&R objective function-----------------------------------------------
    # If you check the R&R objective function (expression 3.1 in Self_study material)
    # it comprises a sum of products in numerator and denominator.
    num_rr = z * (1 - Ki)
    denom_rr = 1 + psi_o * (Ki - 1)
    fob_rr = sum(num_rr / denom_rr)
    # ----------------------------------------------------------------------------------
    # Variables to determine on a PTFlash
    # Liquid- molar composition
    xo = z / (1.0 + (psi_o * (Ki - 1.0)))
    # ----------------------------------------------------------------------------------
    # Vapor molar flow
    V = psi_o * F
    # ----------------------------------------------------------------------------------
    # Liquid molar flow
    L = F - V
    # vapor molar composition of vapour stream at equilibria with xo
    yo = Ki * xo
    # -----------------------------------------------------------------------------------
    # Newton - Raphson method to estimate new value of PSI-------------------------------
    # To estimate a new value of the phase fraction PSI, we implement Newton-Raphson method
    # therefore the derivative of the R&R function is required,hence:
    defo_rr = sum((z * ((Ki - 1.0) ** 2)) / ((1.0 + (psi_o * (Ki - 1.0))) ** 2))
    # Newton-Raphson formula applied to R&R function
    psi_1 = psi_o - (fob_rr / defo_rr)
    # ------------------------------------------------------------------------------------
    # return a new value of PSI, and recalculate all the variables
    psi_o = psi_1
    print(psi_o)

PSI = psi_o
print("PTFlash phase-fraction (PSI):", PSI)
print("Vapour stream  molar flow rate :", V)
print("Liquid stream  molar flow rate :", L)
print("Liquid  molar composition :", xo)
print("Vapour  molar composition :", yo)
