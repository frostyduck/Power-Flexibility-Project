import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
sns.set()


df=pd.DataFrame({u'costs': {0: 19276, 1: 17905},
 u'cr': {0: 0.9449,1: 0.9714},
 u'Optimization models': {0: u'CED',1: u'IHPD'}})

df = df.set_index('Optimization models')

fig = plt.figure(figsize=(7,7)) # Create matplotlib figure

ax = fig.add_subplot(111) # Create matplotlib axes
ax2 = ax.twinx() # Create another axes that shares the same x-axis as a
width = .3

df.costs.plot(kind='bar',color='green', ax=ax, width=width, position=1)
df.cr.plot(kind='bar',color='blue', ax=ax2, width = width, position=0)

ax.grid(None, axis=1)
ax2.grid(None)
ax.set_ylabel('Production сosts')
ax2.set_ylabel('Wind Utilization')

ax.set_ylim(10000, 20000)
ax2.set_ylim(0.75, 1.0)
ax.set_xlim(-1,2)

plt.text(5, .30, r'gamma: $\gamma$', {'color': 'r', 'fontsize': 20})
plt.text(5, .18, r'Omega: $\Omega$', {'color': 'b', 'fontsize': 20})

plt.text(1.42, 0.995, 'Costs', {'color': 'g', 'fontsize': 14}, va="top", ha="right")
plt.text(1.95, 0.985, 'Wind Utilization', {'color': 'b', 'fontsize': 14}, va="top", ha="right")
plt.show()

