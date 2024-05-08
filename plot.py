import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

# First 8 columns are for each individual in a generation, the last is
# the global best so far.
data = np.loadtxt("OUT.methane")

size = int(360//2.5)

A = np.linspace(0, 360, size)
D = np.linspace(0, 5, size)

KA = 0.36
KD = 4.6
A0 = 110.0
D0 = 1.113

AA, DD = np.array(np.meshgrid(A, D))
U = 6 * 0.5 * KA * ((AA - A0) * 2 * np.pi / 360)**2 \
    + 4 * 0.5 * KD * (DD - D0)**2

fig, ax = plt.subplots()

ax.contourf(AA, DD, U, levels=20)
ax.scatter(data[-1][-2], data[-1][-1], c="orange")

def pts(gen, Xi):
    """XI == 0 ⇒ angle, XI == 1 ⇒ distance."""
    d = data[gen]
    return d[Xi:8:2]

inds = ax.scatter(pts(0,0), pts(0,1), c="white")
glob = ax.scatter(data[0][8], data[0][9], c="red")

def update(gen):
    ax.set_title(f"Generation {gen}")
    # zip doesn't work for whatever reason.
    inds.set_offsets(np.stack([pts(gen,0),pts(gen,1)]).T)
    glob.set_offsets([data[gen][8], data[gen][9]])
    return ax, inds, glob

ani = anim.FuncAnimation(fig=fig, func=update, frames=len(data), interval=300)

ani.save(filename="./test.gif", writer="ffmpeg")

# plt.show()

# Local Variables:
# python-shell-interpreter: "ipython3"
# python-shell-interpreter-args: "--simple-prompt"
# End:
