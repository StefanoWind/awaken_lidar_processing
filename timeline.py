# -*- coding: utf-8 -*-
'''
Plot timeline of available files
'''
import os
cd=os.path.dirname(__file__)
import sys
import warnings
from datetime import datetime
import numpy as np
import matplotlib.dates as mdates
import glob
import re
import matplotlib
from matplotlib import pyplot as plt
warnings.filterwarnings('ignore')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['savefig.dpi']=300
plt.close("all")

warnings.filterwarnings('ignore')

#%% Inputs

#users inputs
if len(sys.argv)==1:
    path=os.path.join(cd,'data/awaken/rt1.lidar.z02.a0/*nc')
else:
    path=sys.argv[1]

#%% Initialization
files=glob.glob(path)
dates=[]
print(f"{len(files)} files found in {path}")
os.makedirs(os.path.join(cd,'figures'),exist_ok=True)

#%% Main
for f in files:
    match = re.search(r"(\d{8})\.(\d{6})", f)
    if match:
        date_str, time_str = match.groups()
        dt = datetime.strptime(date_str + time_str, "%Y%m%d%H%M%S")
        
    dates=np.append(dates,dt)
    
    
#%% Plots
plt.figure(figsize=(18,3))
plt.plot(dates,np.zeros(len(dates)),'.g',markersize=1)
plt.xlabel('Time (UTC)')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))
plt.title(path)
plt.grid()
plt.tight_layout()
plt.yticks([])
plt.savefig(os.path.join(cd,'figures','timeline.png'))