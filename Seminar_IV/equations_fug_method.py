import numpy as np
import thermo_properties as prop


# -----------------------------------------------------------------------------------------------------
# function to generate the objective function 1  of underwood approach
# -----------------------------------------------------------------------------------------------------
# inputs :
# tf = Feed temperature, TF
# r = heavy key compound, HK
# l = light key compound, LK
# compounds_feed = name of the compounds (All)  present in the feed stream.
# zFeed = composition in the feed stream , zF
# delta_o = initial guess for delta
# q_feed  = feed state- fraction of liquid in the feed (q)P.
def underwood_1_fab(tf, r, l, compounds_feed, zFeed, delta_o, q_feed):
    alpha_inf_ir = np.zeros(len(compounds_feed))
    Under_sum1_delta = 0
    for u in range(len(compounds_feed)):
        Psat_ir_F = prop.pure_vapor_pressure([compounds_feed[u], r], tf)
        alpha_inf_ir[u] = Psat_ir_F[0] / Psat_ir_F[1]
        if compounds_feed[u] == l:
            alpha_inf_LKHK = alpha_inf_ir[u]
        Under_sum1_delta = Under_sum1_delta + ((alpha_inf_ir[u] * zFeed[u]) / (alpha_inf_ir[u] - delta_o))
    fob_U_delta = Under_sum1_delta - (1.0 - q_feed)
    return fob_U_delta, alpha_inf_ir, alpha_inf_LKHK


# -----------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------
# function to generate the objective function 2  of underwood approach
# -----------------------------------------------------------------------------------------------------
def underwood_2_Rmin(tf, r, compounds_distillate, x_required_D, delta_roots):
    alpha_infd_ir = np.zeros(len(compounds_distillate))
    Under_sum2_rmin = 0
    # for k in range(len(delta_roots)):
    for u in range(len(compounds_distillate)):
        Psat_ird_F = prop.pure_vapor_pressure([compounds_distillate[u], r], tf)
        alpha_infd_ir[u] = Psat_ird_F[0] / Psat_ird_F[1]
        term_R = (alpha_infd_ir[u] * x_required_D[u]) / (alpha_infd_ir[u] - delta_roots)
        Under_sum2_rmin = Under_sum2_rmin + term_R
    R_min = Under_sum2_rmin - 1.0
    return R_min, alpha_infd_ir


# -----------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Bisection method to solve function 1  of underwood approach
# ------------------------------------------------------------------------------------------------------
def solve_underwood_1(tf, r, lo, compounds_feed, zFeed, delta_o, q_feed):
    tol = 1.0e-6
    inc = 1.0e-3
    fob_Uo, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_o, q_feed)
    if fob_Uo > 0.0:
        delta_up = delta_o  # set the upper positive value of delta
        delta_inc = delta_o + inc
        fob_Uinc, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_inc, q_feed)
        if fob_Uinc > fob_Uo:
            while fob_Uo > 0.0:
                delta_o = delta_o - inc
                # % loop to decrease the value of f(x) when decreases with lower values than xo
                fob_Uo, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_o, q_feed)
                delta_low = delta_o
        else:
            while fob_Uo > 0.0:
                delta_o = delta_o + inc
                # loop to decrease the value of f(x) when decreases with higher values than xo
                fob_Uo, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_o, q_feed)
                delta_low = delta_o
    else:  # the function f(x) is negative in xo ; f(xo) < 0
        delta_low = delta_o  # set the lower  (negative)  value of delta
        delta_inc = delta_o + inc
        fob_Uinc, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_inc, q_feed)
        if fob_Uinc > fob_Uo:
            while fob_Uo < 0.0:
                delta_o = delta_o + inc
                # loop to increase the value of f(x) when increases for higher values than xo
                fob_Uo, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_o, q_feed)
                delta_up = delta_o
        else:
            while fob_Uo < 0.0:
                delta_o = delta_o - inc
                fob_Uo, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_o, q_feed)
                delta_up = delta_o
    # ---------------------------------------------------------------------------------------
    # Bisection (secant ) algorithm
    delta_mop = (delta_up + delta_low) / 2.0
    fob_Uo, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_mop, q_feed)
    print(delta_up, delta_low, delta_mop, fob_Uo)
    while abs(fob_Uo) > tol:
        if fob_Uo > 0.0:
            delta_up = delta_mop
            delta_mean = (delta_up + delta_low) / 2.0
            fob_Uo, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_mean, q_feed)
            delta_mop = delta_mean
        else:
            delta_low = delta_mop
            delta_mean = (delta_up + delta_low) / 2.0
            fob_Uo, _, _ = underwood_1_fab(tf, r, lo, compounds_feed, zFeed, delta_mean, q_feed)
            delta_mop = delta_mean

    return delta_mean, fob_Uo
# ------------------------------------------------------------------------------------------------------
