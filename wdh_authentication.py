# -*- coding: utf-8 -*-
'''
Authenticate in MFA
'''
import os
cd=os.path.dirname(__file__)
from doe_dap_dl import DAP
import yaml

#%% Inputs
path_config=os.path.join(cd,'configs/config_awaken.yaml') #config path

#%% Initialization
a2e = DAP('wdh.energy.gov',confirm_downloads=False)
#configs
with open(path_config, 'r') as fid:
    config = yaml.safe_load(fid)

#%% Main
a2e.setup_two_factor_auth(username=config['username'],password=config['password'])
