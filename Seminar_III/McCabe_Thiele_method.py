import numpy as np
import matplotlib.pyplot as plt
import VLE_BubblePoint as BP
import VLE_DewPoint as DP

# python script for the graphical construction
# of the McCabe-Thiele method for the design of binary
# distillation columns.

# -------------Specifications: ----------------------------
# ---------------------------------------------------------
# Binary system
compounds = ["pentane", "heptane"]
# constant pressure across the column
Po = 10.0
# ---------------------------------------------------------
# Composition (more volatile compound MVC)
# Feed composition:
zf = 0.5
# Bottoms composition:
xB = 0.05
# Distillate composition:
xD = 0.95
# Reflux ratio:
R = 3.35
# liquid fraction in the feed:
q = 1.95
# ----------------------------------------------------------
font1 = {'family': 'serif', 'color': 'k', 'size': 15}
# Plot 45o - Line (Y=X)-------------------------------------
xof = np.linspace(0.0, 1.0, 100)
yof = xof
plt.figure(1)
plt.plot(xof, yof, 'k')
plt.xlabel('$x$', fontdict=font1)
plt.ylabel('$y$', fontdict=font1)
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.0)
plt.grid()
# --------------------------------------------------------
# Plot equilibrium curve
# --------------------------------------------------------
# Here, the script invoke the function "VLE_BubblePoint"
# to calculate Y = VLE(X)
Yo = np.zeros((len(compounds), len(xof)))
for i in range(len(xof)):
    xo = [xof[i], 1 - xof[i]]
    T, Y = BP.Bubble_temperature_ideal(compounds, Po, xo)
    Yo[0][i] = Y[0]
    Yo[1][i] = Y[1]
plt.plot(xof, Yo[0])
# ---------------------------------------------------------
# Mark the given composition for feed, distillate & bottoms
# ----------------------------------------------------------
plt.plot(zf, zf, 'ro', xB, xB, 'mo', xD, xD, 'go')
plt.text(zf, 0.4, "z$_F$", fontsize=10, color="red")
plt.text(0.1, xB, "x$_B$", fontsize=10, color="magenta")
plt.text(xD, 0.89, "x$_D$", fontsize=10, color="green")
#plt.text(0.2, 0.8, "Equilibrium Curve", fontsize=10, color="blue")
#plt.text(0.2, 0.75, "$y=VLE(x, P)$", fontsize=10, color="blue", fontweight='bold')
#plt.show()
# ----------------------------------------------------------
# Rectifying operation line equation -----------------------
# It's a linear equation of type:y = m*x + b
m_rec = (R / (R + 1))  # slope (m)
b_rec = (xD / (R + 1))  # Y-intercept (b)
x_rec = np.linspace(xD, 0.0, 100)
y_rec = (m_rec * x_rec) + b_rec
#plt.plot(x_rec, y_rec, '--g')
# ----------------------------------------------------------
# Feed q - line equation -----------------------------------
# It's a linear equation of type: y = m*x + b
m_feed = (q / (q - 1))  # slope (m)
b_feed = (zf / (1 - q))  # Y-intercept (b)
#x_feed = np.linspace(zf, 0.3, 100)
#y_feed = (m_feed * x_feed) + b_feed
#plt.plot(x_feed, y_feed, '--r')
#plt.show()
# ----------------------------------------------------------
# Intersection of q-line and rectifying line
# ----------------------------------------------------------
xP = (b_feed - b_rec) / (m_rec - m_feed)  # intersection in x
# ----------------------------------------------------------
# Plot q-line from zf to intersection point ----------------
x_feed = np.linspace(zf, xP, 100)
y_feed = (m_feed * x_feed) + b_feed
plt.plot(x_feed, y_feed, 'r')
# ----------------------------------------------------------
# Plot rectifying line from xD to intersection point -------
x_rec_int = np.linspace(xD, xP, 100)
y_rec_int = (m_rec * x_rec_int) + b_rec
plt.plot(x_rec_int, y_rec_int, 'g')
# ----------------------------------------------------------
# Plot the intersection point (xp, yp)
yP = (m_rec * xP) + b_rec
plt.plot(xP, yP, 'sk')
# ----------------------------------------------------------
# Stripping operating line equation ------------------------
# It's a linear equation of type: y-y' = m*(x-x')
# ----------------------------------------------------------
m_strip = (yP - xB) / (xP - xB)
x_strip = np.linspace(xB, xP, 100)
y_strip = m_strip * (x_strip - xB) + xB
# plot stripping line:
plt.plot(x_strip, y_strip, 'm')
# plt.legend(
# ["45$^{o}$ - line", "Equilibrium line", "z$_F$", "x$_B$", "x$_D$", "q-line", "Rectifying line", "intercept",
# "Stripping line"])
# ---------------------------------------------------------------------------------------
# Graphical construction of the Staircase to determine th number of equilibrium - stages
# ---------------------------------------------------------------------------------------
# Set a given number of stages ----------------------------------------------------------
N_stages = 5
dg = 0.01
# Total condenser assumption: The staircase starts at xD = yD
x_Mc = [xD, 1 - xD]
yo = x_Mc
xoo = x_Mc[0]
# ----------------------------------------------------------------------------------------
for i in range(N_stages):
    # First stage Staircase: x1 at VLE with YD = XD :
    # the calculation involve the DewPont algorith X= VLE(P,Y)
    T_dew, X = DP.Dew_temperature_ideal(compounds, Po, yo)
    x_eq = X[0]
    plt.plot(x_eq, yo[0], 'yo')
    # ------------- indicates the number of equilibrium stage-----------------------------
    plt.text(x_eq - dg, yo[0] + (1.5 * dg), i + 1, fontsize=10)
    # ------------------------------------------------------------------------------------
    # horizontal lines in the staircase --------------------------------------------------
    plt.plot([xoo, x_eq], [yo[0], yo[0]], 'y')
    xoo = x_eq
    # ------------------------------------------------------------------------------------
    # If the liquid passing stream composition (x) is greater than the intersection point xP;
    # Then the Rectifying operating line must be used, otherwise implement Stripping
    # Operating line
    if xoo > xP:
        y_Mc = (m_rec * x_eq) + b_rec
    else:
        y_Mc = m_strip * (x_eq - xB) + xB
    # -------------------------------------------------------------------------------------
    # vertical lines in the staircase -----------------------------------------------------
    plt.plot([x_eq, x_eq], [yo[0], y_Mc], 'y')
    yo = [y_Mc, 1 - y_Mc]
    x_Mc = yo
    # --------------------------------------------------------------------------------------

plt.title('n-pentane + n-heptane')
plt.text(0.17, 0.05, "Partial Reboiler", fontsize=10)
plt.show()
