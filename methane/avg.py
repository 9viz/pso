import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("glob-methane")
fig, axs = plt.subplots(2, 1)

axs[0].set_title("angle")
axs[0].hist(data[:,0])

axs[1].set_title("distance")
axs[1].hist(data[:,1])
