# -*- coding: utf-8 -*-
'''
Apply LiSBOA to statistics scans

Inputs (both hard-coded and available as command line inputs in this order):
    sdate [%Y-%m-%d]: start date in UTC
    edate [%Y-%m-%d]: end date in UTC
    delete [bool]: whether to delete raw data
    path_config: path to general config file
    mode [str]: serial or parallel
'''
import os
cd=os.path.dirname(__file__)
import sys
import traceback
import warnings
from lisboa import statistics as stats
from datetime import datetime
import yaml
from multiprocessing import Pool
import logging
import re
import xarray as xr
import ast
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib

warnings.filterwarnings('ignore')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 18

warnings.filterwarnings('ignore')

#%% Inputs

#users inputs
if len(sys.argv)==1:
    sdate='2023-02-20' #start date
    edate='2025-03-01' #end date
    delete=False #delete input files?
    replace=True #replace existing files?
    path_config=os.path.join(cd,'configs/config_awaken.yaml') #config path
    mode='serial'#serial or parallel
else:
    sdate=sys.argv[1]
    edate=sys.argv[2] 
    delete=sys.argv[3]=="True"
    replace=sys.argv[4]=="True"
    path_config=sys.argv[5]
    mode=sys.argv[6]#
    
#%% Initalization

#configs
with open(path_config, 'r') as fid:
    config = yaml.safe_load(fid)

#initialize main logger
logfile_main=os.path.join(cd,'log',datetime.strftime(datetime.now(), '%Y%m%d.%H%M%S'))+'_errors.log'
os.makedirs('log',exist_ok=True)

#%% Functions
def _load_config_from_file(config_file: str, source: str):
    """
    Load configuration from an Excel file.

    Args:
        config_file (str): Path to Excel configuration file

    Returns:
        LidarConfig or None: Configuration parameters or None if loading fails
    """
    configs = pd.read_excel(config_file,header=None).set_index(0)
    date_source = np.int64(re.search(r"\d{8}", source).group(0))

    matches = []
    for c in configs.columns:
        regex=configs[c]['regex']
        if "start_date" not in  configs[c]:
            sdate=19700101
        else:
            sdate = configs[c]["start_date"]
        if "end_date" not in  configs[c]:
            edate=30000101
        else:
            edate = configs[c]["end_date"]
        
        match = re.findall(regex, source)
        if len(match) > 0 and sdate <= date_source <= edate:
            matches.append(c)

    if not matches:
        return None
        
    elif len(matches) > 1:
        return None
    
    config=configs[matches[0]].to_dict()
    for s in ['mins','maxs','Dn0']:
        config[s]=list(np.array(ast.literal_eval(config[s]))*config['diameter'])
    del config['regex']
    del config['start_date']
    del config['end_date']
    del config['diameter']

    return config
    
def lisboa_file(file,save_path_stats,config_path,logfile_main,sdate,edate,delete,replace):
    date=re.search(r'\d{8}.\d{6}',file).group(0)[:8]
    if datetime.strptime(date,'%Y%m%d')>=datetime.strptime(sdate,'%Y-%m-%d') and datetime.strptime(date,'%Y%m%d')<=datetime.strptime(edate,'%Y-%m-%d'):
        try:
            logfile=os.path.join(cd,'log',os.path.basename(file).replace('nc','log'))
         
            #load config
            config=_load_config_from_file(config_path,file)
            if config is not None:
                #load data
                Data=xr.open_dataset(file)
                Data=Data.where(Data.qc_wind_speed==0)
                x_exp=[Data.x.values.ravel(),Data.y.values.ravel(),Data.z.values.ravel()]
                f=Data.wind_speed.values.ravel()
                
                lproc=stats.statistics(config)
                grid,Dd,excl,avg,hom=lproc.calculate_statistics(x_exp,f)
                
                if delete:
                    os.remove(file)
        except:
            with open(logfile_main, 'a') as lf:
                lf.write(f"{datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')} - ERROR - Error standardizing file {os.path.basename(file)}: \n")
                traceback.print_exc(file=lf)
                lf.write('\n --------------------------------- \n')

#%% Main
for s in config['channels_stand']:
        
    #standardize all files within date range
    channel=config['channels_stand'][s]
    files=glob.glob(os.path.join(config['path_data'],channel,'*nc'))
    if mode=='serial':
        for f in files:
              lisboa_file(f,None,config['path_config_lisboa'],logfile_main,sdate,edate,delete,replace)
    elif mode=='parallel':
        args = [(files[i],None, config['path_config_lisboa'],logfile_main,sdate,edate,delete,replace) for i in range(len(files))]
        with Pool() as pool:
            pool.starmap(lisboa_file, args)
    else:
        raise BaseException(f"{mode} is not a valid processing mode (must be serial or parallel)")
          
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
        
