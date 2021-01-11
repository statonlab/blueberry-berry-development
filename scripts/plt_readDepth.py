import os
import matplotlib as matpl
if os.environ.get('DISPLAY','') == '':
    print('Currently no display found. Using the non-interactive Agg backend')
    matpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys
file1 = sys.argv[1]
file2 = sys.argv[2]

data1 = pd.read_table(file1, sep='\t', names=['chr','pos','depth'])
data2 = pd.read_table(file2, sep='\t', names=['chr','pos','depth'])
x1 = data1['pos'].values
y1 = data1['depth'].values
x2 = data2['pos'].values
y2 = data2['depth'].values

#fig,axs = plt.subplots(2)
#axs[0].sns.distplot(y1, hist=False)
#axs[0].set(ylabel='read depth')
#axs[1].sns.distplot(y2, hist=False)
#axs[1].set(xlabel='position',ylabel='read depth')


fig, ax = plt.subplots(1,1)
sns.distplot(y1, hist=False)
sns.distplot(y2, hist=False)
plt.xlabel('Read depth')
plt.ylabel('Density')
plt.savefig("histogram.png")
