import numpy as np
#------------------------------------------------------------------------------------------------------
# ---- Matplotlib is a library for creating static, animated, and interactive visualizations in Python. 
# ---- Matplotlib Create  quality plots and figures
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------------------------------
import "name of your SCRIPT where Pbubble is performed" as BPid

x1 = np.linspace(0.0, 1.0, 100)
T = "Temperature of the system"  # K
compounds = ["compound1", "compound2"]
Pbo = np.zeros(len(x1))
Yo = np.zeros((len(compounds), len(x1)))
Ybo = np.zeros((len(compounds), len(x1)))

for i in range(len(x1)):
    xo = [x1[i], 1.0 - x1[i]]
    P, Yb = BPid."function's name where Pbubble is performed"("required inputs")
    Pbo[i] = P
    Ybo[0][i] = Yb[0]
    Ybo[1][i] = Yb[1]


# ---------------------------------------------------------------------------
plt.figure(1)
plt.plot(x1, Pbo, 'b', Ybo[0], Pbo, 'r')
plt.legend(["bubble line", "dew line"])
plt.xlabel('$x_{1}$, $y_{1}$')
plt.ylabel('$P$ / bar')
plt.title('n-pentane(1) / n-heptane(2)')
plt.ylim(1.0, 11.0)
plt.text(0.4, 10, '$T = 400 K$', fontsize=12)
plt.grid()

plt.xlim(0.0, 1.0)

#fig.tight_layout()
plt.show()
