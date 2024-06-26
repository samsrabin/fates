# =======================================================================================
#
# For usage: $python HydroUTestDriver.py --help
#
# This script runs unit tests on the hydraulics functions.
#
#
# =======================================================================================

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
#from matplotlib.backends.backend_pdf import PdfPages
import platform
import numpy as np
import os
import sys
import getopt
import code  # For development: code.interact(local=dict(globals(), **locals()))
import time
import imp
import ctypes
from ctypes import *
from operator import add


CDLParse = imp.load_source('CDLParse','../shared/py_src/CDLParse.py')
F90ParamParse = imp.load_source('F90ParamParse','../shared/py_src/F90ParamParse.py')
PyF90Utils = imp.load_source('PyF90Utils','../shared/py_src/PyF90Utils.py')


from CDLParse import CDLParseDims, CDLParseParam, cdl_param_type
from F90ParamParse import f90_param_type, GetSymbolUsage, GetPFTParmFileSymbols, MakeListUnique
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr

# Load the fortran objects via CTYPES

f90_unitwrap_obj = ctypes.CDLL('bld/UnitWrapMod.o',mode=ctypes.RTLD_GLOBAL)
f90_constants_obj = ctypes.CDLL('bld/FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_wftfuncs_obj = ctypes.CDLL('bld/FatesHydroWTFMod.o',mode=ctypes.RTLD_GLOBAL)
f90_hydrounitwrap_obj = ctypes.CDLL('bld/HydroUnitWrapMod.o',mode=ctypes.RTLD_GLOBAL)

# Alias the F90 functions, specify the return type
# -----------------------------------------------------------------------------------

initalloc_wtfs = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_initallocwtfs
setwrf = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_setwrf
setwkf = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_setwkf
th_from_psi = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapthfrompsi
th_from_psi.restype = c_double
psi_from_th = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrappsifromth
psi_from_th.restype = c_double
dpsidth_from_th = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapdpsidth
dpsidth_from_th.restype = c_double
ftc_from_psi = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapftcfrompsi
ftc_from_psi.restype = c_double
dftcdpsi_from_psi = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapdftcdpsi
dftcdpsi_from_psi.restype = c_double
ftcminwt_from_psi = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapftcminweightfrompsi
ftcminwt_from_psi.restype = c_double

# Some constants
rwccap = [1.0,0.947,0.947,0.947]
pm_leaf = 1
pm_stem = 2
pm_troot = 3
pm_aroot = 4
pm_rhiz = 5

# These parameters are matched with the indices in FATES-HYDRO
vg_type = 1
cch_type = 2
tfs_type = 3

isoil1 = 0  # Top soil layer parameters (@BCI)
isoil2 = 1  # Bottom soil layer parameters

# Constants for rhizosphere
watsat = [0.567, 0.444]
sucsat = [159.659, 256.094]
bsw    = [6.408, 9.27]

unconstrained = True


# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================


class vg_wrf:
    def __init__(self,index,alpha, psd, th_sat, th_res):
        self.alpha = alpha
        self.psd   = psd
        self.th_sat = th_sat
        self.th_res = th_res
        init_wrf_args = [self.alpha, self.psd, self.th_sat, self.th_res]
        iret = setwrf(ci(index),ci(vg_type),ci(len(init_wrf_args)),c8_arr(init_wrf_args))

class cch_wrf:
    def __init__(self,index,th_sat,psi_sat,beta):
        self.th_sat  = th_sat
        self.psi_sat = psi_sat
        self.beta    = beta
        init_wrf_args = [self.th_sat,self.psi_sat,self.beta]
        iret = setwrf(ci(index),ci(cch_type),ci(len(init_wrf_args)),c8_arr(init_wrf_args))

class vg_wkf:
    def __init__(self,index,alpha, psd, th_sat, th_res, tort):
        self.alpha  = alpha
        self.psd    = psd
        self.th_sat = th_sat
        self.th_res = th_res
        self.tort   = tort
        init_wkf_args = [self.alpha, self.psd,self.th_sat,self.th_res,self.tort]
        iret = setwkf(ci(index),ci(vg_type),ci(len(init_wkf_args)),c8_arr(init_wkf_args))

class cch_wkf:
    def __init__(self,index,th_sat,psi_sat,beta):
        self.th_sat  = th_sat
        self.psi_sat = psi_sat
        self.beta    = beta
        init_wkf_args = [self.th_sat,self.psi_sat,self.beta]
        iret = setwkf(ci(index),ci(cch_type),ci(len(init_wkf_args)),c8_arr(init_wkf_args))


class tfs_wrf:
    def __init__(self,index,th_sat,th_res,pinot,epsil,rwc_fd,cap_corr,cap_int,cap_slp,pmedia):
        self.th_sat = th_sat
        self.th_res = th_res
        self.pinot  = pinot
        self.epsil  = epsil
        self.rwc_fd = rwc_fd
        self.cap_corr = cap_corr
        self.cap_int  = cap_int
        self.cap_slp  = cap_slp
        self.pmedia   = pmedia
        init_wrf_args = [self.th_sat,self.th_res,self.pinot,self.epsil,self.rwc_fd,self.cap_corr,self.cap_int,self.cap_slp,self.pmedia]
        iret = setwrf(ci(index),ci(tfs_type),ci(len(init_wrf_args)),c8_arr(init_wrf_args))

class tfs_wkf:
    def __init__(self,index,p50,avuln):
        self.avuln = avuln
        self.p50   = p50
        init_wkf_args = [self.p50,self.avuln]
        iret = setwkf(ci(index),ci(tfs_type),ci(len(init_wkf_args)),c8_arr(init_wkf_args))


def OMParams(zsoi):

    zsapric= 0.5
    om_watsat         = max(0.93 - 0.1   *(zsoi/zsapric), 0.83)
    om_bsw            = min(2.7  + 9.3   *(zsoi/zsapric), 12.0)
    om_sucsat         = min(10.3 - 0.2   *(zsoi/zsapric), 10.1)

    #om_sucsat = SuctionMMtoMPa(om_sucsat)
    
    return(om_watsat,om_sucsat,om_bsw)
    
def CCHParmsCosby84T5(zsoi,om_frac,sand_frac,clay_frac):

    # Cosby, B.J., Hornberger, G.M., Clapp, R.B., and Ginn, T.R. 1984.
    # A statistical exploration of the relationships of soil moisture
    # characteristics to the physical properties of soils. Water Resour.
    # Res. 20:682-690.

    # cosby_1984_table5

    # Get pedotransfer for soil matrix
    watsat = 0.489 - 0.00126*sand_frac
    bsw    = 2.91 + 0.159*clay_frac
    sucsat = 10. * ( 10.**(1.88-0.0131*sand_frac) ) 

    [om_watsat,om_sucsat,om_bsw] = OMParams(zsoi)
    
    # Update pedotransfer to include organic material
    watsat  = (1.0 - om_frac) * watsat  + om_watsat * om_frac
    bsw     = (1.0 - om_frac) * bsw     + om_bsw * om_frac
    sucsat  = (1.0 - om_frac) * sucsat  + om_sucsat * om_frac

    # Convert from mm to MPa
    sucsat = SuctionMMtoMPa(sucsat)

    
    # psi = psi_sat*(th/th_sat)**(-beta)
    # th  = th_sat*(psi / psi_sat)^(-1/beta)
    
    return(watsat,sucsat,bsw)


def SuctionMMtoMPa(suction_mm):

    denh2o     = 1.0e3    # kg/m3
    grav_earth = 9.8      # m/s2
    mpa_per_pa = 1.0e-6   # MPa per Pa
    m_per_mm   = 1.0e-3   # meters per millimeter
    
    suction_mpa = (-1.0)*suction_mm*denh2o*grav_earth*mpa_per_pa*m_per_mm

    return(suction_mpa)
                
def main(argv):


    version = platform.python_version()
    verlist = version.split('.')


    # Read in the arguments
    # =======================================================================================

#    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
#    parser.add_argument('--cdl-file', dest='cdlfile', type=str, \
#                        help="Input CDL filename.  Required.", required=True)

#    args = parser.parse_args()


    # Set number of analysis points
    npts = 1000


    #    min_theta = np.full(shape=(2),dtype=np.float64,fill_value=np.nan)
    #    wrf_type = [vg_type, vg_type, cch_type, cch_type]
    #    wkf_type = [vg_type, tfs_type, cch_type, tfs_type]
    
    #    th_ress = [0.01, 0.10, -9, -9]
    #    th_sats = [0.55, 0.55, 0.65, 0.65]
    #    alphas  = [1.0, 1.0, 1.0, 1.0]
    #    psds    = [2.7, 2.7, 2.7, 2.7]
    #    tort    = [0.5, 0.5, 0.5, 0.5]
    #    beta    = [-9, -9, 6, 9]
    #    avuln   = [2.0, 2.0, 2.5, 2.5]
    #    p50     = [-1.5, -1.5, -2.25, -2.25]


    ncomp = 4
    ncomp_tot = 20
    
    rwc_fd  = [1.0,0.958,0.958,0.958]
    rwccap  = [1.0,0.947,0.947,0.947]
    cap_slp = []
    cap_int = []
    cap_corr= []
    hydr_psi0 = 0.0
    hydr_psicap = -0.6

    for pm in range(4):
        if (pm == 0):
            cap_slp.append(0.0)
            cap_int.append(0.0)
            cap_corr.append(1.0)
        else:
            cap_slp.append((hydr_psi0 - hydr_psicap )/(1.0 - rwccap[pm]))
            cap_int.append(-cap_slp[pm] + hydr_psi0)
            cap_corr.append(-cap_int[pm]/cap_slp[pm])


    # Allocate memory to our objective classes
    iret = initalloc_wtfs(ci(ncomp_tot),ci(ncomp_tot))
    print('Allocated')


    #min_psi = -10.
    #min_psi_falloff = 1
    #test_psi = np.linspace(min_psi, 0, num=1000)
    #weight = y = e^(2*(-10-x)) 
    
    

    # Define the funcions and their parameters
#    vg_wrf(1,alpha=1.0,psd=2.7,th_sat=0.55,th_res=0.1)
#    vg_wkf(1,alpha=1.0,psd=2.7,th_sat=0.55,th_res=0.1,tort=0.5)

    cch_wrf(1,th_sat=0.55, psi_sat=-1.56e-3, beta=3)
    cch_wkf(1,th_sat=0.55, psi_sat=-1.56e-3, beta=3)

#    cch_wrf(3,th_sat=0.55, psi_sat=-1.56e-3, beta=6)
#    tfs_wkf(3,p50=-2.25, avuln=2.0)

    names=['Soil','ARoot','Stem','Leaf']

    theta_sat = [0.55,0.65,0.65,0.75]
    theta_res = [0.15,0.16,0.21,0.11]

    # Absorbing root
    tfs_wrf(2,th_sat=theta_sat[1],th_res=theta_res[1],pinot=-1.043478, \
            epsil=8,rwc_fd=rwc_fd[3],cap_corr=cap_corr[3], \
            cap_int=cap_int[3],cap_slp=cap_slp[3],pmedia=4)
    tfs_wkf(2,p50=-2.25, avuln=2.0)

    # Stem
    tfs_wrf(3,th_sat=theta_sat[2],th_res=theta_res[2],pinot=-1.22807, \
            epsil=10,rwc_fd=rwc_fd[2],cap_corr=cap_corr[2], \
            cap_int=cap_int[2],cap_slp=cap_slp[2],pmedia=2)
    tfs_wkf(3,p50=-2.25, avuln=4.0)

    # Leaf
    tfs_wrf(4,th_sat=theta_sat[3],th_res=theta_res[3],pinot=-1.465984, \
            epsil=12,rwc_fd=rwc_fd[0],cap_corr=cap_corr[0], \
            cap_int=cap_int[0],cap_slp=cap_slp[0],pmedia=1)
    tfs_wkf(4,p50=-2.25, avuln=2.0)

    print('initialized WRF')

#    theta = np.linspace(0.10, 0.7, num=npts)

    theta = np.full(shape=(ncomp,npts),dtype=np.float64,fill_value=np.nan)
    psi   = np.full(shape=(ncomp,npts),dtype=np.float64,fill_value=np.nan)
    dpsidth = np.full(shape=(ncomp,npts),dtype=np.float64,fill_value=np.nan)
    cdpsidth = np.full(shape=(ncomp,npts),dtype=np.float64,fill_value=np.nan)



    for ic in range(ncomp):
        theta[ic,:] = np.linspace(theta_res[ic], 1.2*theta_sat[ic], num=npts)
        for i in range(npts):
            psi[ic,i] = psi_from_th(ci(ic+1),c8(theta[ic,i]))


    # Theta vs psi plots

    fig0, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
        ax1.plot(theta[ic,:],psi[ic,:],label='{}'.format(names[ic]))

    ax1.set_ylim((-10,1))
    ax1.set_ylabel('Psi [MPa]')
    ax1.set_xlabel('VWC [m3/m3]')
    ax1.legend(loc='lower right')

    for ic in range(ncomp):
        for i in range(npts):
            dpsidth[ic,i]  = dpsidth_from_th(ci(ic+1),c8(theta[ic,i]))
        for i in range(1,npts-1):
            cdpsidth[ic,i] = (psi[ic,i+1]-psi[ic,i-1])/(theta[ic,i+1]-theta[ic,i-1])

    # Theta vs dpsi_dth (also checks deriv versus explicit)

    fig1, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
        ax1.plot(theta[ic,],dpsidth[ic,],label='func')
        ax1.plot(theta[ic,],cdpsidth[ic,],label='check')
    #ax1.set_ylim((0,1000))

    ax1.set_ylabel('dPSI/dTh [MPa m3 m-3]')
    ax1.set_xlabel('VWC [m3/m3]')
    ax1.legend(loc='upper right')

    fig11, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
        ax1.plot(theta[ic,],1.0/dpsidth[ic,],label='{}'.format(names[ic]))

    ax1.set_ylabel('dTh/dPSI/ [m3 m-3 MPa-1]')
    ax1.set_xlabel('VWC [m3/m3]')
    ax1.legend(loc='upper right')


    # Push parameters to WKF classes
    # -------------------------------------------------------------------------
    # Generic VGs

    ftc   = np.full(shape=(ncomp,npts),dtype=np.float64,fill_value=np.nan)
    dftcdpsi = np.full(shape=(ncomp,npts),dtype=np.float64,fill_value=np.nan)
    cdftcdpsi = np.full(shape=(ncomp,npts),dtype=np.float64,fill_value=np.nan)

    for ic in range(ncomp):
        for i in range(npts):
            ftc[ic,i] = ftc_from_psi(ci(ic+1),c8(psi[ic,i]))

            if( (ftc[ic,i]>0.9) and (theta[ic,i]<0.4) ):
                print('tpf: ',theta[ic,i],psi[ic,i],ftc[ic,i])

    for ic in range(ncomp):
        for i in range(npts):
            dftcdpsi[ic,i]  = dftcdpsi_from_psi(ci(ic+1),c8(psi[ic,i]))
        for i in range(1,npts-1):
            cdftcdpsi[ic,i] = (ftc[ic,i+1]-ftc[ic,i-1])/(psi[ic,i+1]-psi[ic,i-1])


    fig2, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
        ax1.plot(psi[ic,:],ftc[ic,:],label='{}'.format(names[ic]))

    ax1.set_ylabel('FTC')
    ax1.set_xlabel('Psi [MPa]')
    ax1.set_xlim([-10,3])

    ax1.legend(loc='upper left')


    # FTC versus theta

    fig4, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
        ax1.plot(theta[ic,:],ftc[ic,:],label='{}'.format(names[ic]))

    ax1.set_ylabel('FTC')
    ax1.set_xlabel('Theta [m3/m3]')
    ax1.legend(loc='upper left')

    # dFTC/dPSI

    fig3,ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
#        ax1.plot(psi[ic,:],abs(dftcdpsi[ic,:]-cdftcdpsi[ic,:])/abs(cdftcdpsi[ic,:]),label='{}'.format(ic))
        ax1.plot(psi[ic,:],dftcdpsi[ic,:],label='{}'.format(names[ic]))

    ax1.set_ylabel('dFTC/dPSI')
    ax1.set_xlabel('Psi [MPa]')
    ax1.set_xlim([-10,3])
    ax1.set_ylim([0,2])
#    ax1.set_ylim([0,10])
    ax1.legend(loc='upper right')


    ndepth = 1000
    zdepth = np.linspace(0.01,10.0,num = ndepth)

    om_watsat_z = np.zeros(shape=(ndepth,1))
    om_sucsat_z = np.zeros(shape=(ndepth,1))
    om_bsw_z    = np.zeros(shape=(ndepth,1))

    for icl in range(ndepth):
        [om_watsat_z[icl],om_sucsat_z[icl],om_bsw_z[icl]] = OMParams(zdepth[icl])

    fig4,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(9,6))
    
    ax1.plot(om_watsat_z,-zdepth)
    ax2.plot(om_sucsat_z,-zdepth)
    ax3.plot(om_bsw_z,-zdepth)
    ax4.axis('off')

    
    ax1.set_ylabel('soil depth')
    ax1.set_xlabel('Saturated WC [m3/m3]')
    ax2.set_xlabel('Saturated Suction [mm]')
    ax3.set_ylabel('soil depth')
    ax3.set_xlabel('beta ')
    #    ax1.set_xlim([-10,3])
    #    ax1.set_ylim([0,2])
    #    ax1.set_ylim([0,10])
    #    ax1.legend(loc='upper right')
    plt.tight_layout()

        
    # Soil texture distributions

    #clay_frac = np.zeros(shape=(100,1))
    #sand_frac = np.zeros(shape=(100,1))
    #om_frac   = np.zeros(shape=(100,1))
    #clay_frac = np.linspace(0.0, 0.7, num=npts)

    watsat_v = []
    sucsat_v = []
    bsw_v    = []
    
    ntex = 5
    for icl in range(ntex):
        clay_frac = float(icl)/float(ntex)
        for isa in range(ntex):
            sand_frac =  float(isa)/float(ntex)
            
            if( (clay_frac+sand_frac)<1.0 ):
                for om_frac in [0.0,1.0]:
                    for zsoi in [0.01,10.0]:
                        [watsat,sucsat,bsw] = CCHParmsCosby84T5(zsoi,om_frac,sand_frac,clay_frac)
                        watsat_v.append(watsat)
                        sucsat_v.append(sucsat)
                        bsw_v.append(bsw)
                

    fig5,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(9,6))
    
    ax1.hist(watsat_v,bins=50)
    ax2.hist(sucsat_v,bins=50)
    ax3.hist(bsw_v,bins=50)
    ax4.axis('off')

    
    ax1.set_xlabel('Sat. WC [m3/m3]')
    ax2.set_xlabel('Sat. Sucation [MPa]')
    ax3.set_xlabel('Beta')

    plt.tight_layout()

    ind = []
    ind.append(np.argmin(bsw_v))
    ind.append(np.argmax(bsw_v))
    ind.append(np.argmin(sucsat_v))
    ind.append(np.argmax(sucsat_v))
    ind.append(np.argmin(watsat_v))
    ind.append(np.argmax(watsat_v))
    
    # Now lets test the lowest possible water contents using these 

    #psi_v = -np.linspace(0.01,24,24)
    psi_v = -np.logspace(-3,1.5,60)
    fig10,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(6,9))
    
    
    ii=4
    th_v = []
    ftc_v = []
    ftc_min_wt_v = []
    for i in ind:
        ii = ii+1
        cch_wrf(ii,th_sat=watsat_v[i], psi_sat=sucsat_v[i], beta=bsw_v[i])
        cch_wkf(ii,th_sat=watsat_v[i], psi_sat=sucsat_v[i], beta=bsw_v[i])
        print('---th sat: {}, psi sat: {}, beat: {}'.format(watsat_v[i],sucsat_v[i],bsw_v[i]))
        for psi in psi_v:
            th_v.append(th_from_psi(ci(ii),c8(psi)))
            ftc_v.append(ftc_from_psi(ci(ii),c8(psi)))
            ftc_min_wt_v.append(ftcminwt_from_psi(ci(ii),c8(psi)))
                
        ax1.plot(psi_v,th_v)
        ax2.plot(psi_v,ftc_v)
        ax3.plot(psi_v,ftc_min_wt_v)
        th_v = []
        ftc_v = []
        ftc_min_wt_v = []
        
    ax2.set_xlabel('Suction MPa')
    ax1.set_ylabel('Theta')
    ax2.set_ylabel('FTC')
    ax3.set_ylabel('FTC Min Weight')
    
    plt.show()

#    code.interact(local=dict(globals(), **locals()))

# Helper code to plot negative logs

def semilogneg(x):

    y = np.sign(x)*np.log(abs(x))
    return(y)

def semilog10net(x):

    y = np.sign(x)*np.log10(abs(x))
    return(y)


# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
