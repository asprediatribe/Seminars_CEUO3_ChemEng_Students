import thermo_properties as prop
import numpy as np
import VLE_BubblePoint as BP
import math as mt
import equations_fug_method as eq
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

P =  # conversion Psi to bar
compounds = []
Ti_sat = prop.pure_vapour_temperature(compounds, P)
print("Tsat:", Ti_sat)
HK =
LK =
# ----------------------------------------------------------------------------------------------
# Feed specifications, xF/ f
# ----------------------------------------------------------------------------------------------
compounds_in_F = []
f_iC4 =
f_nC4 =
f_iC5 =
f_nC5 =
f_C6 =
f_C7 =
f_C8 =
f_C9 =
# molar feed fluxes , fi
f = np.array([f_iC4, f_nC4, f_iC5, f_nC5, f_C6, f_C7, f_C8, f_C9])
# distillate molar fraction xD,i
zF = f / sum(f)
# ---------------------------------------------------------------------------
print("Feed molar fractions, zF:", zF, sum(zF))
# ---------------------------------------------------------------------------

# --------------------------------------------------------------
# Distillate specifications, XD/ d
# --------------------------------------------------------------
compounds_in_D = []
d_iC4 =  # lbmol/h
d_nC4 =
d_iC5 =
d_nC5 =
# molar distillate fluxes , di
d = np.array([d_iC4, d_nC4, d_iC5, d_nC5])
# ditillate molar fraction xD,i
xD = d / sum(d)
# --------------------------------------------------------------
print("Distillate molar fractions, xD:", xD, sum(xD))
# --------------------------------------------------------------
# Bottoms  specifications, XB/ b
# --------------------------------------------------------------
compounds_in_B = []
b_nC4 =
b_iC5 =
b_nC5 =
b_C6 =
b_C7 =
b_C8 =
b_C9 =
# molar bottoms fluxes , bi
b = np.array([b_nC4, b_iC5, b_nC5, b_C6, b_C7, b_C8, b_C9])
# distillate molar fraction xD,i
xB = b / sum(b)
# --------------------------------------------------------------
print("Bottoms molar fractions, xB:", xB, sum(xB))
# --------------------------------------------------------------
# --------------------------------------------------------------
# Calculation of  temperatures at stage 1 and N:
# --------------------------------------------------------------
# in this code the nomenclature is as follows:
# stage N (top of the column)
# stage 1 (Bottom of the column)-partial reboiler
# Notice that, this nomenclature is not rigid, some books named
# the stages in the other way around (i.e N-bottom ; 1-top)
# Here we follow the nomenclature used in the self_study material
# Chapter 4.
# --------------------------------------------------------------
# -- Establish Temperature in the top --------------------------
Td_refo, YD = BP.VLE_Tbubble(xD, P, compounds_in_D)
# --------------------------------------------------------------
P_bubble_D, YD = BP.VLE_Pbubble(xD, Td_refo, compounds_in_D)
print("Tbubble reflux_drum, TD (K):", Td_refo, YD)
# ---------------------------------------------------------------
# Establish bottom temperature TB:
# ---------------------------------------------------------------
T_bubble_B, YB = BP.VLE_Tbubble(xB, P, compounds_in_B)
print("Tbubble of bottoms, TB (K)", T_bubble_B)
# ---------------------------------------------------------------
# relative volatility of key compounds
# ---------------------------------------------------------------
Psat_keys_N = prop.pure_vapour_pressure([LK, HK], Td_refo)
print("Psat_Keys_N", Psat_keys_N)
Psat_keys_1 = prop.pure_vapour_pressure([LK, HK], T_bubble_B)
print("Psat_Keys_1", Psat_keys_1)
alpha_ij_N = Psat_keys_N[0] / Psat_keys_N[1]
alpha_ij_1 = Psat_keys_1[0] / Psat_keys_1[1]
print("alphaN", alpha_ij_N)
print("alpha1", alpha_ij_1)
alpha_m_ij = np.sqrt(alpha_ij_N * alpha_ij_1)
print("alphaM", alpha_m_ij)
# --------------------------------------------------------------------------------
# ----Fenske equation for minimum equilibrium stages ----Nmin---------------------
keyss = [] # put the corresponding inputs
ratio_keys_db = np.zeros(len(keyss))
for i in range(len(keyss)):
    index_kd = compounds_in_D.index(keyss[i])
    index_kb = compounds_in_B.index(keyss[i])
    ratio_keys_db[i] = xD[index_kd] / xB[index_kb]
Nmin = np.log10(ratio_keys_db[0] * (1.0 / ratio_keys_db[1])) / np.log10(alpha_m_ij)
print("Minimum number of stages", Nmin)
# --------------------------------------------------------------------------------
# Distribution of No-key components ----------------------------------------------
# --------------------------------------------------------------------------------
b_i_profile = np.zeros(len(compounds_in_F))
alpha_m_iHK = np.zeros(len(compounds_in_F))
index_HK_d = compounds_in_D.index(HK)
index_HK_b = compounds_in_B.index(HK)
dr = d[index_HK_d]
br = b[index_HK_b]
j = 0
for i in range(len(compounds_in_F)):
    Psat_NK_N = prop.pure_vapour_pressure([compounds_in_F[i], HK], Td_refo)
    Psat_NK_1 = prop.pure_vapour_pressure([compounds_in_F[i], HK], T_bubble_B)
    # ---------------------------------------------------------------------------
    # print("compo",compounds_in_F[i])
    alpha_ir_N = Psat_NK_N[0] / Psat_NK_N[1]
    alpha_ir_1 = Psat_NK_1[0] / Psat_NK_1[1]
    alpha_m_iHK[i] = np.sqrt(alpha_ir_N * alpha_ir_1)
    if compounds_in_F[i] != LK and compounds_in_F[i] != HK:
        dr = d[index_HK_d]
        br = b[index_HK_b]
        denom_fenske = (dr / br) * mt.pow(alpha_m_iHK[i], Nmin)
        b_i_profile[i] = f[i] / (1.0 + denom_fenske)
    else:
        index_kb = compounds_in_B.index(keyss[j])
        b_i_profile[i] = b[index_kb]
        j = j + 1
print("profile b", b_i_profile)
print("sum", sum(b_i_profile))
d_i_profile = f - b_i_profile
print("profile d", d_i_profile)
print("sum", sum(d_i_profile))
xd_i_profile = d_i_profile / sum(d_i_profile)
xb_i_profile = b_i_profile / sum(b_i_profile)
print("recalculated compositions in D:", xd_i_profile)
print("recalculated compositions in B:", xb_i_profile)
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# Solution of the First Underwood Equation
# ---------------------------------------------------------------------------------
# liquid fraction in the feed:
q =
# Feed temperature
TF =
# calculation of relative volativity at the pinch point (feed)
# r= HK
# delta_o = fisrt initial guess (value) to build first underwood equation
delta_oo =
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# The following lines plot the first underwood equation as a function of "phi",
# the "absorption factor", to visualize the shape and roots of the function
# ---------------------------------------------------------------------------------
fobo_delta, alpha_inf_ir, alpha_in_LK_HK = eq.underwood_1_fab(TF, HK, LK, compounds_in_F, zF, delta_oo, q)
print("FOBO", fobo_delta, "alpha_inf", alpha_inf_ir)

delta_feas = np.linspace(alpha_inf_ir[len(compounds_in_F) - 2], alpha_inf_ir[1], 1000)
fun_graph = np.zeros(len(delta_feas))
for i in range(len(delta_feas)):
    fun_graph[i], _, _ = eq.underwood_1_fab(TF, HK, LK, compounds_in_F, zF, delta_feas[i], q)

font1 = {'family': 'serif', 'color': 'k', 'size': 15}
plt.figure(1)
plt.plot(delta_feas, fun_graph, 'b')
plt.xlabel('$\delta$', fontdict=font1)
plt.ylabel('$U1(\delta)$', fontdict=font1)
# plt.xlim(0.0, 1.0)
plt.ylim(-2.5, 2.5)
plt.grid()
# ------------------------------------------------------------------------------------------------
# ---------------------------------- end of the plot Underwood 1 vs Phi---------------------------
# Here the code proceed with the root calculation
# Root finder procedure, calculation of "phi" ----------------------------------------------------
# m_root = potential roots os the absorption factor "delta"
# NCd = number of distributing components
NCd = 2
m_root = NCd - 1
#roots_solv_u1 = np.zeros(m_root)
voc = 0.015
# print("delta_limit", alpha_in_LK_HK)
# guess_delta_o = [1.0 + voc, alpha_in_LK_HK - voc]
guess_delta_o = [1.0 + voc]
delta_Root = fsolve(eq.solve_underwood_delta, guess_delta_o, args=(TF, HK, LK, compounds_in_F, zF, q))
roots_solv_u1 = delta_Root

# ------------------------------------------------------------------------------------------------
print("roots", roots_solv_u1)
# , "Fobo", fob_U1)
# ------------------------------------------------------------------------------------------------
# ----------- Minimum reflux calculation using the root(s) "phi"----------------------------------
R_min_Underwood, _ = eq.underwood_2_Rmin(TF, HK, compounds_in_D, xD, roots_solv_u1[0])
print("Minimum Reflux:", R_min_Underwood)
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# Guilliland Correlation for Actual Reflux Ratio
# and Equilibrium Stages
# ------------------------------------------------------------------------------------------------------
R =
X_gil = (R - R_min_Underwood) / (R + 1.0)
# molokanov's equation
molok_o = (1.0 + (54.4 * X_gil)) / (11.0 + (117.2 * X_gil))
molok_1 = (X_gil - 1.0) / np.sqrt(X_gil)
Y = 1.0 - np.exp(molok_o * molok_1)
N_stages_gil = (Y + Nmin) / (1.0 - Y)
print(" Guilliland X:", X_gil)
print(" Guilliland Y:", Y)
print("Ideal Number of Stages by Guilliland:", N_stages_gil)
# -----------------------------------------------------------------------------------------------------
plt.show()
