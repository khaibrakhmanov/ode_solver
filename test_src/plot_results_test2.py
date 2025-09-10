# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/usr/bin/env python

#import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.style.use('classic')

fsize=9
lw=1.0
ms=3.5

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Times New Roman'

plt.rcParams.update({'axes.titlesize':fsize})
plt.rcParams.update({'axes.labelsize':fsize})
plt.rcParams.update({'font.size':fsize})
plt.rcParams.update({'xtick.labelsize':fsize})
plt.rcParams.update({'ytick.labelsize':fsize})
plt.rcParams.update({'legend.fontsize':fsize})

data = np.loadtxt("../msvs_project/data.dat")
t = data[:, 0]
q = data[:, 1]
q_an = data[:, 2]

inch2cm = 2.54
fig_width_cm = 11
fig_height_cm = 11#fig_width_cm
fig, axs = plt.subplots(1, 1, figsize=(fig_width_cm/inch2cm, fig_height_cm/inch2cm), dpi=200)

#------------------------------------------------------------------------------

axs.set_xscale("log")
axs.set_yscale("log")
axs.set_xlabel(r'$r$'+u' [au]', size=fsize)
axs.set_ylabel(r'$T$ ' + u'[K]', size=fsize)
axs.plot(t, q, 'o-', linewidth=lw, markersize=ms, label=u'num')

# def qa(r):
#     return np.exp(-1000*r)
# ts=np.linspace(0, 0.1, 2000)#10**logts
# axs.plot(ts, qa(ts), '--', linewidth=lw, markersize=ms, label=u'an')

axs.plot(t, q_an, '--', linewidth=lw, markersize=ms, label=u'an')

# axs.set_ylim(1e-45, 1.1)
axs.legend(loc='best')

fig.show()