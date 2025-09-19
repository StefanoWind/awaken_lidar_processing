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
import socket
import getpass
import matplotlib.pyplot as plt
import matplotlib
plt.close('all')

warnings.filterwarnings('ignore')
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['savefig.dpi'] = 300

#%% Inputs

#users inputs
if len(sys.argv)==1:
    sdate='2023-02-20' #start date
    edate='2025-12-31' #end date
    delete=False #delete input files?
    replace=False #replace existing files?
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
    
    for s in ['mins','maxs','Dn0','plot_locations','u_limits']:
        config[s]=list(np.array(ast.literal_eval(config[s]))*config['diameter'])
    
    config_lisboa=config.copy()
    del config_lisboa['regex']
    del config_lisboa['start_date']
    del config_lisboa['end_date']
    del config_lisboa['diameter']
    del config_lisboa['hub_height']
    del config_lisboa['plot_locations']
    del config_lisboa['data_level_in']
    del config_lisboa['data_level_out']
    del config_lisboa['u_limits']
    
    return config,config_lisboa
    
def lisboa_file(file,save_path_stats,config_path,logfile_main,sdate,edate,delete,replace):
    date=re.search(r'\d{8}.\d{6}',file).group(0)[:8]
    if datetime.strptime(date,'%Y%m%d')>=datetime.strptime(sdate,'%Y-%m-%d') and datetime.strptime(date,'%Y%m%d')<=datetime.strptime(edate,'%Y-%m-%d'):
        try:
            logfile=os.path.join(cd,'log',os.path.basename(file).replace('nc','log'))
         
            #load config
            config,config_lisboa=_load_config_from_file(config_path,file)
            save_path=file.replace(config['data_level_in'],config['data_level_out'])
            if not os.path.isfile(save_path) or replace:
                
                if config is not None:
                    #load data
                    Data=xr.open_dataset(file)
                    time=Data.time
                    Data=Data.where(Data.qc_wind_speed==0)
                    
                    if len(config['Dn0'])==3:
                        x_exp=[Data.x.values.ravel(),Data.y.values.ravel(),Data.z.values.ravel()]
                    else:
                        x_exp=[Data.x.values.ravel(),Data.y.values.ravel()]
                    
                    #de-projection
                    proj=Data.x/Data.range
                    f=(Data.wind_speed/proj).values.ravel()
                    
                    #thresholding
                    f[f<config['u_limits'][0]]=np.nan
                    f[f>config['u_limits'][1]]=np.nan
                    
                    #run LiSBOA
                    lproc=stats.statistics(config_lisboa,logfile=logfile)
                    grid,Dd,excl,avg,hom=lproc.calculate_statistics(x_exp,f,2)
                    avg[avg<config['u_limits'][0]]=np.nan
                    avg[avg>config['u_limits'][1]]=np.nan
                    
                    #%% Output
                    Output=xr.Dataset()
                    if len(grid)==3:
                        coords={'x':grid[0],'y':grid[1],'z':grid[2]}
                    else:
                        coords={'x':grid[0],'y':grid[1]}
                    Output['u_avg']=xr.DataArray(avg,coords=coords,
                                                 attrs={'units':'m/s','description':'mean streamwise velocity'})
                    Output['u_std']=xr.DataArray(hom**0.5,coords=coords,
                                                 attrs={'units':'m/s','description':'std of streamwise velocity'})

                    Output.attrs['start_time']=str(time.isel(beamID=0,scanID=0).values)
                    Output.attrs['end_time']=  str(time.isel(beamID=-1,scanID=-1).values)
                    
                    for c in config:
                        Output.attrs[f'config_{c}']=config[c]
                    
                    Output.attrs["data_level"]=config['data_level_out']
                    Output.attrs['input_source']=os.path.basename(file)
                    Output.attrs["contact"]= "stefano.letizia@nrel.gov"
                    Output.attrs["institution"]= "NREL"
                    Output.attrs["description"]= "Statistics of de-projected wind speed calculated through LiSBOA"
                    Output.attrs["reference"]= "Letizia et al. LiSBOA (LiDAR Statistical Barnes Objective Analysis) for optimal design of lidar scans and retrieval of wind statistics – Part 1: Theoretical framework. AMT, 14, 2065–2093, 2021, 10.5194/amt-14-2065-2021"
                    Output.attrs["history"]= (
                        f"Generated by {getpass.getuser()} on {socket.gethostname()} on "
                        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} using {os.path.basename(sys.argv[0])}"
                    )
                    Output.attrs["code1"]="https://github.com/NREL/FIEXTA/tree/main/lisboa"
                    Output.attrs["code2"]="https://github.com/StefanoWind/awaken_lidar_processing" 
                    Output.attrs["location_meaning"]=Data.attrs['location_meaning']
                    Output.attrs["location_id"]=Data.attrs['location_id']
                    
                    os.makedirs(os.path.dirname(save_path),exist_ok=True)
                    Output.to_netcdf(save_path)
                    
                    visualization(Output,config,save_path)
                    
                    if delete:
                        os.remove(file)
        except:
            with open(logfile_main, 'a') as lf:
                lf.write(f"{datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')} - ERROR - Error processing file {os.path.basename(file)}: \n")
                traceback.print_exc(file=lf)
                lf.write('\n --------------------------------- \n')
                
                
def visualization(Data,config,save_path):
    
    
    #Extrac data at hub-height plane
    x=Data.x.values
    y=Data.y.values
    
    if 'z' in Data.coords:
        z=Data.z.values
        u_avg_int=Data['u_avg'].sel(z=0,method='nearest').values
        u_std_int=Data['u_std'].sel(z=0,method='nearest').values
    else:
        u_avg_int=Data['u_avg'].values
        u_std_int=Data['u_std'].values
    TI_int=u_std_int.T/np.abs(u_avg_int.T)*100

    #xlims
    if config['mins'][0]<0:
        xlim=[config['mins'][0]/config['diameter'],0.5]
    else:
        xlim=[-0.5,config['maxs'][0]/config['diameter']]
        
    #contour levels
    u_avg=Data['u_avg'].values
    TI=Data['u_std'].values/u_avg*100
    levels_u=np.unique(np.round(np.linspace(np.nanpercentile(u_avg,5)-0.5, np.nanpercentile(u_avg,95)+0.5, 20),1))
    levels_TI=np.unique(np.round(np.linspace(np.nanpercentile(TI,5)-0.5, np.nanpercentile(TI,95)+0.5, 20),1))
    
    #Plot mean velocity at hub_height
    fig=plt.figure(figsize=(18,6.6))
    ax = plt.subplot(1,2,1)
    ax.set_facecolor((0,0,0,0.2))
    
    cf=plt.contourf(x/config['diameter'],y/config['diameter'],u_avg_int.T,levels_u, cmap='coolwarm',extend='both')
    plt.contour(x/config['diameter'],y/config['diameter'],u_avg_int.T,levels_u, colors='k',linewidths=1,alpha=0.25,extend='both')
    plt.xlim(xlim)
    plt.ylim([config['mins'][1]/config['diameter'],config['maxs'][1]/config['diameter']])
    plt.grid(alpha=0.5)
    ax.set_aspect('equal')
    plt.xlabel(r'$x/D$')
    plt.ylabel(r'$y/D$')
    plt.title('Mean streamwise velocity on '+Data.attrs['start_time'][:10]+'\n File: '+os.path.basename(Data.attrs['input_source'])\
              +'\n Time (UTC): '+Data.attrs['start_time'][11:19]+' - '+Data.attrs['end_time'][11:19])
    fig.subplots_adjust(left=0.1, right=0.9, wspace=0.4)
    draw_turbine_side(0, 0, 1, 270)
    
    cax=fig.add_axes([ax.get_position().x0+ax.get_position().width+0.01,ax.get_position().y0,0.015,ax.get_position().height])
    plt.colorbar(cf, cax=cax,label=r'Mean streamwise velocity [m s$^{-1}$]')

    #plot TI at hub_height
    ax = plt.subplot(1,2,2)
    ax.set_facecolor((0,0,0,0.2))
    cf=plt.contourf(x/config['diameter'],y/config['diameter'],TI_int,levels_TI,cmap='hot',extend='both')
    plt.contour(x/config['diameter'],y/config['diameter'],TI_int,levels_TI, colors='k',linewidths=1,alpha=0.25,extend='both')
    plt.xlim(xlim)
    plt.ylim([config['mins'][1]/config['diameter'],config['maxs'][1]/config['diameter']])
    plt.grid(alpha=0.5)
    ax.set_aspect('equal')
    plt.xlabel(r'$x/D$')
    plt.ylabel(r'$y/D$')
    plt.title('Turbulence intensity on '+Data.attrs['start_time'][:10]+'\n File: '+os.path.basename(Data.attrs['input_source'])\
              +'\n Time (UTC): '+Data.attrs['start_time'][11:19]+' - '+Data.attrs['end_time'][11:19])
    draw_turbine_side(0, 0, 1, 270)
    cax=fig.add_axes([ax.get_position().x0+ax.get_position().width+0.01,ax.get_position().y0,0.015,ax.get_position().height])
    plt.colorbar(cf, cax=cax,label=r'Turbulence intensity [%]')
    
    fig.savefig(save_path.replace('.nc','_hub_height.png'))
    plt.close()
    
    #y-z planes (for 3D scans only)
    if 'z' in Data.coords:
        
        fig=plt.figure(figsize=(18,10))
        ctr=1
        for x_plot in config['plot_locations']:
            u_avg_int=Data['u_avg'].interp(x=[x_plot],method='linear').values.squeeze()
            u_std_int=Data['u_std'].interp(x=[x_plot],method='linear').values.squeeze()
            TI_int=u_std_int/np.abs(u_avg_int)*100
            
            #plot mean velocity
            ax = plt.subplot(len(config['plot_locations']),2,(ctr-1)*2+1)
            if ctr==1:
                plt.title('Mean streamwise velocity on '+Data.attrs['start_time'][:10]+'\n File: '+os.path.basename(Data.attrs['input_source'])\
                          +'\n Time (UTC): '+Data.attrs['start_time'][11:19]+' - '+Data.attrs['end_time'][11:19]\
                          +'\n' +r'$x/D='+str(x_plot/config['diameter'])+'$')
            else:
                plt.title(r'$x/D='+str(x_plot/config['diameter'])+'$')
                
            if ctr<len(config['plot_locations']):
                ax.set_xticklabels([])
            else:
                plt.xlabel(r'$y/D$') 
                    
            ax.set_facecolor((0,0,0,0.2))
            cf1=plt.contourf(y/config['diameter'],z/config['diameter'],u_avg_int.T,levels_u, cmap='coolwarm',extend='both')
            plt.contour(y/config['diameter'],z/config['diameter'],u_avg_int.T,levels_u, colors='k',linewidths=1,alpha=0.25,extend='both')
            plt.grid(alpha=0.5)

            draw_turbine_front(0,0,1,0)
            ax.fill_between(y,config['mins'][2]/config['diameter'], 
                            y*0-config['hub_height']/config['diameter'],color='k')
            
            plt.ylabel(r'$z/D$')
            
            plt.xlim([config['mins'][1]/config['diameter'],config['maxs'][1]/config['diameter']])
            plt.ylim([config['mins'][2]/config['diameter'],config['maxs'][2]/config['diameter']])
            ax.set_aspect('equal')
            
            #plot TI
            ax = plt.subplot(len(config['plot_locations']),2,(ctr-1)*2+2)
            if ctr==1:
                plt.title('Turbulence intensity on '+Data.attrs['start_time'][:10]+'\n File: '+os.path.basename(Data.attrs['input_source'])\
                          +'\n Time (UTC): '+Data.attrs['start_time'][11:19]+' - '+Data.attrs['end_time'][11:19]\
                          +'\n' +r'$x/D='+str(x_plot/config['diameter'])+'$')
            else:
                plt.title(r'$x/D='+str(x_plot/config['diameter'])+'$')
            
            if ctr<len(config['plot_locations']):
                ax.set_xticklabels([])
            else:
                plt.xlabel(r'$y/D$') 
                
            ax.set_facecolor((0,0,0,0.2))
            cf2=plt.contourf(y/config['diameter'],z/config['diameter'],TI_int.T,levels_TI, cmap='hot',extend='both')
            plt.contour(y/config['diameter'],z/config['diameter'],TI_int.T,levels_TI, colors='k',linewidths=1,alpha=0.25,extend='both')
            plt.grid(alpha=0.5)
            
            draw_turbine_front(0,0,1,0)
            ax.fill_between(y,config['mins'][2]/config['diameter'], 
                            y*0-config['hub_height']/config['diameter'],color='k')
            
            plt.ylabel(r'$z/D$')
            
            plt.xlim([config['mins'][1]/config['diameter'],config['maxs'][1]/config['diameter']])
            plt.ylim([config['mins'][2]/config['diameter'],config['maxs'][2]/config['diameter']])
            ax.set_aspect('equal')
            
            ctr+=1

        axs=fig.axes
        cax = fig.add_axes([axs[-2].get_position().x0+axs[-2].get_position().width+0.0075,axs[-2].get_position().y0,
                            0.01,axs[0].get_position().y0+axs[0].get_position().height-axs[-2].get_position().y0])
        plt.colorbar(cf1,cax=cax,label=r'Mean streamwise velocity [m s$^{-1}$]')
        cax = fig.add_axes([axs[-1].get_position().x0+axs[-1].get_position().width+0.0075,axs[-1].get_position().y0,
                            0.01,axs[1].get_position().y0+axs[1].get_position().height-axs[-1].get_position().y0])
        plt.colorbar(cf2,cax=cax,label=r'Turbulence intensity [%]')
        fig.subplots_adjust(hspace=0.4,wspace=0.4)
        
        fig.savefig(save_path.replace('.nc','_cross-planes.png'))
        plt.close()
        
def draw_turbine_side(x,y,D,yaw):
    import matplotlib.image as mpimg
    from matplotlib import transforms
    from matplotlib import pyplot as plt
    img = mpimg.imread('Turbine_side.png')
    ax=plt.gca()
    xlim=ax.get_xlim()
    ylim=ax.get_ylim()
    tr = transforms.Affine2D().scale(D/1800,D/1800).translate(-30*D/700,-350*D/700).rotate_deg(90-yaw).translate(x,y)
    ax.imshow(img, transform=tr + ax.transData)
    plt.xlim(xlim)
    plt.ylim(ylim)
    
def draw_turbine_front(x,y,D,yaw=0):
    import matplotlib.image as mpimg
    from matplotlib import transforms
    from matplotlib import pyplot as plt
    
    img = mpimg.imread('Turbine_front.png')
    ax=plt.gca()
    xlim=ax.get_xlim()
    ylim=ax.get_ylim()
    tr = transforms.Affine2D().translate(-500,-570).scale(1/530/2,1/530/2).scale(D*np.cos(np.radians(yaw)),D).translate(0,0).rotate_deg(180)
    ax.imshow(img, transform=tr + ax.transData, zorder=2)
    plt.xlim(xlim)
    plt.ylim(ylim)
    

#%% Main
for s in config['channels_stand']:
        
    #standardize all files within date range
    channel=config['channels_stand'][s]
    files=sorted(glob.glob(os.path.join(config['path_data'],channel,'*nc')))
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
        
