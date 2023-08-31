#!/usr/bin/env python
# coding: utf-8

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys

def ketc(slicer,grating,grat_wave,f_lam_index,seeing,exposure_time,ccd_bin='1x1',mag_AB=False,flux=False,Nframes=1,nas=False,spatial_bin=[1,1],spectral_bin=1,sb=False,plot=True,plotout=False,emline_width=False,printout=False,wrange=False):
    bin_factor=1.
    if ccd_bin=='2x2': bin_factor=0.25
    if ccd_bin=='2x2' and slicer=='S': print('******** WARNING: DO NOT USE 2x2 BINNING WITH SMALL SLICER')
    read_noise={'B':2.7,'R':3.0}[grating[0]] #electrons
    nas_overhead = 10. # seconds per half cycle
    seeing1 = seeing
    seeing2 = seeing
    pixels_per_arcsec = 1./0.147
    if slicer=='L': 
        seeing2 = 1.38
        snr_spatial_bin = seeing1*seeing2
        pixels_spectral = 8
        arcsec_per_slice = 1.35
    elif slicer=='M': 
        seeing2 = max([0.69,seeing])
        snr_spatial_bin = seeing1*seeing2
        pixels_spectral = 4
        arcsec_per_slice = 0.69
    elif slicer=='S': 
        seeing2 = seeing
        snr_spatial_bin = seeing1*seeing2
        pixels_spectral = 2
        arcsec_per_slice = 0.35
    N_slices = seeing/arcsec_per_slice
    if spatial_bin[0]>1 or spatial_bin[1]>1:
        N_slices = spatial_bin[1]/arcsec_per_slice
        snr_spatial_bin = spatial_bin[0]*spatial_bin[1]
    pixels_spatial_bin = pixels_per_arcsec * N_slices
    A_per_pixel={'BL':0.625,'BM':0.28,'BH1':0.125,'BH2':0.125,'BH3':0.125,'RL':1.0,'RM1':0.24,'RM2':0.25,'RH1':0.125,'RH2':0.125,'RH3':0.125,'RH4':0.125}[grating]
    
    print(f'f_lam ~ lam^{f_lam_index}')
    print(f'reference wavelength = {grat_wave}')
    print(f'SLICER: {slicer}')
    print(f'GRATING: {grating}')
    print(f'SEEING: {seeing} arcsec')
    print(f'Å/pixel: {A_per_pixel}')
    print(f'spectral pixels in 1 spectral resolution element = {pixels_spectral}')
    A_per_spectral_bin = pixels_spectral*A_per_pixel
    print(f'Å/resolution element: {A_per_spectral_bin}')
    if spectral_bin>1: snr_spectral_bin = spectral_bin
    else: snr_spectral_bin = A_per_spectral_bin
    print(f'Å/SNR bin: {snr_spectral_bin}')
    pixels_per_snr_spec_bin = snr_spectral_bin/A_per_pixel
    print(f'Pixels/Spectral SNR bin: {pixels_per_snr_spec_bin}')
    print(f'SNR Spatial Bin [arcsec^2]: {snr_spatial_bin}')
    print(f'SNR Spatial Bin [pixels^2]: {pixels_spatial_bin}')
    
    if emline_width: flux = flux/emline_width
    if not flux and emline_width: 
        print('Do not use mag_AB for emission line')
        return
    if mag_AB and ~sb: flux = (10**(-0.4*(mag_AB+48.6)))*(3e18/grat_wave)/grat_wave/snr_spatial_bin
    if mag_AB and sb: flux = (10**(-0.4*(mag_AB+48.6)))*(3e18/grat_wave)/grat_wave
        
    w, p_A = make_obj(flux,grat_wave,f_lam_index,grating)
        
    if mag_AB and ~sb:
        flux_input = 'mag_AB'
        print(f'OBJECT mag: {mag_AB} {flux_input}')
    if mag_AB and sb:
        flux_input = 'mag_AB / arcsec^2'
        print(f'OBJECT mag: {mag_AB} {flux_input}')
    if flux and ~sb and ~bool(emline_width): flux_input = 'erg cm^-2 s^-1 Å^-1'
    if flux and ~sb and emline_width: flux_input = f'erg cm^-2 s^-1 Å^-1 over {emline_width:4.1f} Å'
    if flux and sb and ~bool(emline_width): flux_input = 'erg cm^-2 s^-1 Å^-1 arcsec^-2'
    if flux and sb and emline_width: flux_input = f'erg cm^-2 s^-1 Å^-1 arcsec^-2 over {emline_width:4.1f} Å'
    if flux: print(f'OBJECT Flux {flux} {flux_input}')
    if emline_width: print('EMISSION LINE OBJECT --> flux is not per unit Å')
    
    if not nas:
        c_o = obj_cts(w,p_A,grating,exposure_time)*snr_spatial_bin*snr_spectral_bin
        c_s = sky_cts(w,grating,exposure_time)*snr_spatial_bin*snr_spectral_bin
        c_r = Nframes*read_noise**2*pixels_per_snr_spec_bin*pixels_spatial_bin*bin_factor
        snr = c_o/np.sqrt(c_s+c_o+c_r)
    else:
        n_cyc = int((exposure_time-nas_overhead)/2./(nas+nas_overhead)+0.5)
        total_exposure = (2*n_cyc*(nas+nas_overhead))+nas_overhead
        print(f'NAS: Rounding up to {n_cyc} Cycles of NAS for total exposure of {total_exposure} s')
        exposure_time = n_cyc*nas
        c_o = obj_cts(w,p_A,grating,exposure_time)*snr_spatial_bin*snr_spectral_bin
        c_s = sky_cts(w,grating,exposure_time)*snr_spatial_bin*snr_spectral_bin
        c_r = 2.*Nframes*read_noise**2*pixels_per_snr_spec_bin*pixels_spatial_bin*bin_factor
        snr = c_o/np.sqrt(2.*c_s+c_o+c_r)
    
    snr[np.isnan(snr)]=0.
    rat = c_o/c_s
    rat[np.isinf(rat)]=0.
    rat[np.isnan(rat)]=0.
    rat[rat<0]=0.
    
    if plot or plotout:
        fig,ax=plt.subplots(6,1,figsize=(10,12))
        fig.suptitle('SNR')
        xrange = [np.min(w),np.max(w)]
        if wrange: xrange=wrange
        lw=1.0
        ax[0].plot(w,snr,lw=lw)
        ax[0].set_ylabel(f'SNR / {snr_spectral_bin:4.1f} Å')
        ax[1].plot(w,c_o,lw=lw)
        ax[1].set_ylabel(f'Obj cts / {snr_spectral_bin:4.1f} Å')
        ax[2].plot(w,c_s,lw=lw)
        ax[2].set_ylabel(f'Sky cts / {snr_spectral_bin:4.1f} Å')
        ax[3].plot(w,c_r*np.ones_like(w))
        ax[3].set_ylabel(f'Read Noise cts / {snr_spectral_bin:4.1f} Å')
        ax[4].plot(w,rat,lw=lw)
        ax[4].set_ylabel('Obj cts / Sky cts')
        ax[5].plot(w,p_A,lw=lw)
        ax[5].set_ylabel(r'Flux [ph cm$^{-2}$ s$^{-1}$ Å${^-1}$]')
        ax[5].set_xlabel('Wavelength (Å)')
        for i in ax: i.set_xlim(*xrange)
        plt.tight_layout()
        if plotout: plt.savefig(f'{plotout}.pdf',bbox_inches='tight')
        plt.show()
    
    if printout:
        with open(printout+'.txt', 'w') as f:
            f.write("Wave(Å)      SNR        cts_obj     cts_sky     cts_rn      obj/sky     p_lam\n")
            for i in range(len(w)):
                f.write(f"{w[i]:<12.1f}{snr[i]:<12.1f}{c_o[i]:<12.3f}{c_s[i]:<12.3f}{c_r:<12.3f}{rat[i]:<12.3f}{p_A[i]:<12.3e}\n")
    return w,snr

def make_obj(flux, grat_wave, f_lam_index, grat):
    w={'B':np.arange(3000,6000),'R':np.arange(5000,11000)}[grat[0]]
    p_A = flux/(2e-8/w)*(w/grat_wave)**f_lam_index
    return w, p_A
    
def obj_cts(w, fo, grat, exposure_time):
    A_geo = np.pi/4.*(1e3)**2
    eff = inst_throughput(w, grat)
    cts = eff*A_geo*exposure_time*fo
    return cts

def inst_throughput(wave, grat):
    eff_bl = np.array([0.1825,0.38,0.40,0.46,0.47,0.44])
    eff_bm = np.array([0.1575, 0.33, 0.36, 0.42, 0.48, 0.45])
    eff_bh1 = np.array([0., 0.0, 0.0, 0.0, 0.0, 0.])
    eff_bh2 = np.array([0.,  0.18, 0.3, 0.4, 0.28, 0.])
    eff_bh3 = np.array([0., 0., 0., 0.2, 0.29, 0.31])
    wave_bl = np.array([355., 530.])*10.
    wave_bm = np.array([355., 530.])*10.
    wave_bh1 = np.array([350., 450.])*10.
    wave_bh2 = np.array([405., 486.])*10.
    wave_bh3 = np.array([405., 530.])*10.
    wave_rl = np.array([540.,1100.])*10.
    wave_rm1 = np.array([540.,1100.])*10.
    wave_rm2 = np.array([540.,1100.])*10.
    wave_rh1 = np.array([540.,800.])*10.
    wave_rh2 = np.array([540.,900.])*10.
    wave_rh3 = np.array([540.,1100.])*10.
    wave_rh4 = np.array([540.,1100.])*10.
    wave_r,e_rl,erm1,e_rm1,erm2,e_rm2,erh1,e_rh1,erh2,e_rh2,erh3,e_rh3,erh4,e_rh4=np.loadtxt('KCRM_eff.csv',skiprows=1,delimiter=',').T
    if grat[0]=='B': 
        trans_atmtel = np.array([0.54, 0.55, 0.56, 0.56, 0.56, 0.55])
        wave_0 = np.array([355.,380.,405.,450.,486.,530.])*10.
    elif grat[0]=='R':
        wave_0 = wave_r*10.
        trans_atmtel = 0.56
    else: 
        print(f'Grating {grat} not recognized.')
        return 0.*wave
    eff=trans_atmtel*{'BL':eff_bl,'BM':eff_bm,'BH1':eff_bh1,'BH2':eff_bh2,'BH3':eff_bh3,'RL':e_rl,'RM1':e_rm1,'RM2':e_rm2,'RH1':e_rh1,'RH2':e_rh2,'RH3':e_rh3,'RH4':e_rh4}[grat]
    wave_range={'BL':wave_bl,'BM':wave_bm,'BH1':wave_bh1,'BH2':wave_bh2,'BH3':wave_bh3,'RL':wave_rl,'RM1':wave_rm1,'RM2':wave_rm2,'RH1':wave_rh1,'RH2':wave_rh2,'RH3':wave_rh3,'RH4':wave_rh4}[grat]
    eff_int = np.interp(wave,wave_0,eff)
    eff_int[(wave<wave_range[0]) | (wave>wave_range[1])] = 0.
    return eff_int

def sky_cts(w, grat, exposure_time, airmass=1.2, area=1.):
    A_geo = np.pi/4.*(1e3)**2
    eff = inst_throughput(w, grat)
    cts = eff*A_geo*exposure_time*sky_mk(w, grat)*airmass*area 
    return cts

def sky_mk(wave, grat):
    with fits.open('lris_esi_skyspec_fnu_uJy.fits') as temp:
        f_nu=temp[0].data
        dw=float(temp[0].header['CDELT1'])
        w0=float(temp[0].header['CRVAL1'])
        ws=np.arange(len(f_nu))*dw+w0
        f_lam = f_nu*1e-29*3e18/ws/ws
        p_lam = f_lam/(2e-8/ws)
#         p_lam[p_lam<0]=0
        ps_int = np.interp(wave,ws,p_lam)
    return  ps_int

if __name__=='__main__':

	args=sys.argv[1:]
	
	if len(args)==6: 
		print('You must specify a flux or magnitude')
		quit()
	else:
		reqargs=args[:6]
		for i in [2,3,4,5]:
			reqargs[i]=float(reqargs[i])
		optional=args[6:]
		argsdict={'ccd_bin':'1x1','mag_AB':False,'flux':False,'Nframes':1,'nas':False,'spatial_bin':[1,1],'spectral_bin':1,'sb':False,'plot':True,'plotout':False,'emline_width':False,'printout':False,'wrange':False}
		for arg in optional:
			arglist=arg.strip().split('=')
			if arglist[0] in ['mag_AB','flux','Nframes','nas','spectral_bin','emline_width']:
				argsdict[arglist[0]]=float(arglist[1])
			elif arglist[0]=='spatial_bin':
				argsdict[arglist[0]]=[float(i) for i in arglist[1].split(',')]
			elif arglist[0] in ['sb','plot']:
				if arglist[1]=='False': argsdict[arglist[0]]=False
				if arglist[1]=='True': argsdict[arglist[0]]=True
			else: argsdict[arglist[0]]=arglist[1]
		ketc(*reqargs,ccd_bin=argsdict['ccd_bin'],mag_AB=argsdict['mag_AB'],flux=argsdict['flux'],Nframes=argsdict['Nframes'],nas=argsdict['nas'],spatial_bin=argsdict['spatial_bin'],spectral_bin=argsdict['spectral_bin'],sb=argsdict['sb'],plot=argsdict['plot'],plotout=argsdict['plotout'],emline_width=argsdict['emline_width'],printout=argsdict['printout'])

quit()