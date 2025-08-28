# -*- coding: utf-8 -*-
'''
Authenticate in MFA
'''
from doe_dap_dl import DAP
import sys

#%% Inputs
if len(sys.argv)==1:
    login='mfa'
else:
    login=sys.argv[1]

#%% Main
a2e = DAP('wdh.energy.gov',confirm_downloads=False)

if login=='mfa':
    a2e.setup_two_factor_auth()
elif login=='basic':
    a2e.setup_cert_auth()
else:
    raise ValueError(r'login type should be mfa or basic, not {login}')