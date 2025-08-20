# -*- coding: utf-8 -*-
'''
Authenticate in MFA
'''
from doe_dap_dl import DAP


a2e = DAP('a2e.energy.gov',confirm_downloads=False)
password=input('Password: ')
a2e.setup_two_factor_auth(username='sletizia', password=password)
