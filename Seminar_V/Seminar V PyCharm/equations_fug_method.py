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
        Psat_ir_F = prop.pure_vapour_pressure([compounds_feed[u], r], tf)
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
        Psat_ird_F = prop.pure_vapour_pressure([compounds_distillate[u], r], tf)
        alpha_infd_ir[u] = Psat_ird_F[0] / Psat_ird_F[1]
        term_R = (alpha_infd_ir[u] * x_required_D[u]) / (alpha_infd_ir[u] - delta_roots)
        Under_sum2_rmin = Under_sum2_rmin + term_R
    R_min = Under_sum2_rmin - 1.0
    return R_min, alpha_infd_ir


def solve_underwood_delta(delta_o, tf, r, l, compounds_feed, zFeed, q_feed):
    delta, _, _ = underwood_1_fab(tf, r, l, compounds_feed, zFeed, delta_o, q_feed)
    return delta
