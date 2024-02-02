import numpy as np
# ------------------------------------------------------------------------------------------------------
# ---- Matplotlib is a library for creating static, animated, and interactive visualizations in Python.
# ---- Matplotlib Create  quality plots and figures
import matplotlib.pyplot as plt
# ------------------------------------------------------------------------------------------------------
import VLE_BubblePoint as BPid

x1 = np.linspace(0.0, 1.0, 100)
T = 300  # K
compounds = ["pentane", "heptane"]
Pbo = np.zeros(len(x1))
Yo = np.zeros((len(compounds), len(x1)))
Ybo = np.zeros((len(compounds), len(x1)))

for i in range(len(x1)):
    xo = [x1[i], 1.0 - x1[i]]
    P, Yb = BPid.Bubble_pressure_ideal(compounds, T, xo)
    Pbo[i] = P
    Ybo[0][i] = Yb[0]
    Ybo[1][i] = Yb[1]

# ---------------------------------------------------------------------------
#print(Pbo)
plt.figure(1)
#plt.plot(x1,  Ybo[0],  'b', x1, x1, 'k')
plt.plot(x1, Pbo, 'b', Ybo[0], Pbo, 'r')
#plt.legend(["bubble line", "dew line"])
#plt.xlabel('$x_{1}$, $y_{1}$')
plt.xlabel('$x_{1}$')
plt.ylabel('y')
#plt.ylabel('$y_{1}$')
plt.title('n-pentane(1) / n-heptane(2)')
# plt.ylim(1.0, 11.0)
#plt.text(0.4, 10, '$T = 400 K$', fontsize=12)
plt.grid()

plt.xlim(0.0, 1.0)

# fig.tight_layout()
plt.show()
