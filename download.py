# -*- coding: utf-8 -*-
'''
Download lidar data

Inputs (both hard-coded and available as command line inputs in this order):
    t_start [%Y-%m-%d]: start date in UTC
    t_end [%Y-%m-%d]: end date in UTC
    download [bool]: whether to download new data
    path_config: path to general config file
'''
import os
cd=os.path.dirname(__file__)
import sys
import warnings
from datetime import datetime
from datetime import timedelta
import yaml
from doe_dap_dl import DAP
from matplotlib import pyplot as plt

warnings.filterwarnings('ignore')

#%% Inputs

#users inputs
if len(sys.argv)==1:
    t_start='2023-11-21' #start date
    t_end='2023-11-22' #end date
    mfa=False
    path_config=os.path.join(cd,'configs/config_awaken.yaml') #config path
else:
    t_start=sys.argv[1] #start date
    t_end=sys.argv[2]  #end date
    mfa=sys.argv[3]=="True" #use MFA on WDH
    path_config=os.path.join(cd,'config',sys.argv[4])#config path
    
#%% Initalization
print(f'Downloading lidar data from {t_start} to {t_end}: MFA={mfa}, config={path_config}.')

#configs
with open(path_config, 'r') as fid:
    config = yaml.safe_load(fid)

#DAP setup
a2e = DAP('a2e.energy.gov',confirm_downloads=False)
if mfa:
    a2e.setup_two_factor_auth(username=config['username'], password=config['password'])
else:
    a2e.setup_cert_auth(username=config['username'], password=config['password'])

N_periods=(datetime.strptime(t_end, '%Y-%m-%d')-datetime.strptime(t_start, '%Y-%m-%d'))/timedelta(hours=config['time_increment'])
time_bin=[datetime.strptime(t_start, '%Y-%m-%d') + timedelta(hours=config['time_increment']*x) for x in range(int(N_periods)+1)]

#%% Main
for t1,t2 in zip(time_bin[:-1],time_bin[1:]):
    for c in config['channels']:
        channel=config['channels'][c]
        save_path=os.path.join(config['path_data'],channel.replace('a0','b0'))
        
        if config['ext1'][c]=='':
            _filter = {
                'Dataset': channel,
                'date_time': {
                    'between':  [datetime.strftime(t1, '%Y%m%d%H%M%S'),
                                 datetime.strftime(t2-timedelta(seconds=1), '%Y%m%d%H%M%S')]
                },
                'file_type': config['format'][c]}
        else:
            _filter = {
                'Dataset': channel,
                'date_time': {
                    'between':  [datetime.strftime(t1, '%Y%m%d%H%M%S'),
                                 datetime.strftime(t2-timedelta(seconds=1), '%Y%m%d%H%M%S')]
                },
                'file_type': config['format'][c],
                'ext1':config['ext1'][c], 
            }
        
        a2e.download_with_order(_filter, path=os.path.join(cd,'data',channel), replace=False)
                
        
        
