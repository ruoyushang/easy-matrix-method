
import os
import sys,ROOT
import array
import math
from array import *
from ROOT import *
from astropy import units as my_unit
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.time import Time
from scipy import special
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import cycle
from scipy import fftpack
from astropy import wcs
from astropy.io import fits

import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

import CommonPlotFunctions

fig, ax = plt.subplots()
figsize_x = 8.6
figsize_y = 6.4
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

new_nbins_x = 100
new_nbins_y = 100
#new_nbins_x = 80
#new_nbins_y = 80
#new_nbins_x = 50
#new_nbins_y = 50

effective_area_cut = 10000.
energy_bin = CommonPlotFunctions.energy_bin
energy_bin_cut_low = int(sys.argv[4])
energy_bin_break = int(sys.argv[5])
energy_bin_cut_up = int(sys.argv[6])
source_name = sys.argv[1]
input_epoch = sys.argv[2] # 'V5' or 'V6' or 'V5V6'
isON = sys.argv[3]  # 'ON' or 'OFF'

doImposter = False
if isON:
    doImposter = True

doBiasCorrect = False
#doBiasCorrect = True
#if not doImposter:
#    doBiasCorrect = False

analysis_method = CommonPlotFunctions.analysis_method
n_xoff_bins = CommonPlotFunctions.n_xoff_bins
n_yoff_bins = CommonPlotFunctions.n_yoff_bins
smooth_size_spectroscopy = CommonPlotFunctions.smooth_size_spectroscopy

n_imposters = 5
if not doImposter:
    n_imposters = 0


plot_tag = source_name
plot_tag += '_'+analysis_method
plot_tag += '_E'+sys.argv[4]+'_'+sys.argv[6]

list_epoch = []
if 'V5' in input_epoch:
    list_epoch += ['V5']
if 'V6' in input_epoch:
    list_epoch += ['V6']

total_data_expo = 0.

MSCW_lower_blind = -0.5
MSCL_lower_blind = -0.6
MSCW_upper_blind = 0.6
MSCL_upper_blind = 0.4
n_extra_lower_bins = 0
n_extra_upper_bins = 6
mtx_dim_w_fine = 6
mtx_dim_l_fine = 6
MSCW_bin_size = (MSCW_upper_blind-MSCW_lower_blind)/(mtx_dim_w_fine)
MSCL_bin_size = (MSCL_upper_blind-MSCL_lower_blind)/(mtx_dim_l_fine)
MSCW_plot_upper_fine = MSCW_upper_blind+n_extra_upper_bins*MSCW_bin_size
MSCL_plot_upper_fine = MSCL_upper_blind+n_extra_upper_bins*MSCL_bin_size
MSCW_plot_lower_fine = MSCW_lower_blind-n_extra_lower_bins*MSCW_bin_size
MSCL_plot_lower_fine = MSCL_lower_blind-n_extra_lower_bins*MSCL_bin_size

def SaveFITS(hist_input):

    image_length_x = hist_input.GetNbinsX()
    image_length_y = hist_input.GetNbinsY()
    array_shape = (image_length_y,image_length_x)
    new_image_data = np.zeros(array_shape)

    ref_pix = 1
    central_ra = hist_input.GetXaxis().GetBinCenter(ref_pix)
    central_dec = hist_input.GetYaxis().GetBinCenter(ref_pix)
    pixel_size = hist_input.GetXaxis().GetBinCenter(ref_pix+1)-hist_input.GetXaxis().GetBinCenter(ref_pix)
    reduced_wcs = wcs.WCS(naxis=2)
    reduced_wcs.wcs.ctype = ['RA---TAN','DEC--TAN']
    reduced_wcs.wcs.crpix = [ref_pix,ref_pix] # Set the reference pixel coordinates
    reduced_wcs.wcs.crval = [central_ra, central_dec] # Set the reference values for the physical coordinates
    reduced_wcs.wcs.cdelt = [pixel_size,pixel_size]

    for pix_x in range(0,new_image_data.shape[0]):
        for pix_y in range(0,new_image_data.shape[1]):
            new_image_data[pix_x,pix_y] += hist_input.GetBinContent(pix_y+1,pix_x+1)

    filename = hist_input.GetName()
    fits.writeto('output_plots/%s.fits'%(filename), new_image_data, reduced_wcs.to_header(), overwrite=True)

def FillSkyMapHistogram(hist_input,hist_output,scale=1.):

    temp_nbins_x = hist_output.GetNbinsX()
    temp_nbins_y = hist_output.GetNbinsY()
    temp_map_left = hist_output.GetXaxis().GetBinLowEdge(1)
    temp_map_right = hist_output.GetXaxis().GetBinLowEdge(hist_output.GetNbinsX()+1)
    temp_map_lower = hist_output.GetYaxis().GetBinLowEdge(1)
    temp_map_upper = hist_output.GetYaxis().GetBinLowEdge(hist_output.GetNbinsY()+1)
    hist_temp = ROOT.TH2D("hist_temp","",temp_nbins_x,temp_map_left,temp_map_right,temp_nbins_y,temp_map_lower,temp_map_upper)

    for binx in range(0,hist_input.GetNbinsX()):
        for biny in range(0,hist_input.GetNbinsY()):
            bin_center_x = hist_input.GetXaxis().GetBinCenter(binx+1)
            bin_center_y = hist_input.GetYaxis().GetBinCenter(biny+1)
            bin_content = hist_input.GetBinContent(binx+1,biny+1)
            bin_error = hist_input.GetBinError(binx+1,biny+1)
            bx2 = hist_temp.GetXaxis().FindBin(bin_center_x)
            by2 = hist_temp.GetYaxis().FindBin(bin_center_y)
            hist_temp.SetBinContent(bx2,by2,bin_content)
            hist_temp.SetBinError(bx2,by2,bin_error)

    hist_output.Add(hist_temp,scale)

def GetMapChi2Distribution(hist_z,hist_cnt):

    count_max = hist_cnt.GetMaximum()
    zscores = []
    for binx in range(0,hist_z.GetNbinsX()):
        for biny in range(0,hist_z.GetNbinsY()):
            content = hist_z.GetBinContent(binx+1,biny+1)
            count = hist_cnt.GetBinContent(binx+1,biny+1)
            #if count/count_max<0.5: continue
            if count==0.: continue
            zscores += [content]
    return zscores

def GetMapNormalDistribution(hist_mean):

    zscores = []
    hist_mean_skymap = ROOT.TH2D("hist_mean_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    hist_noise_skymap = ROOT.TH2D("hist_noise_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    hist_mean_skymap.Add(hist_mean)

    n_trials = 10
    for trial in range(0,n_trials):
        hist_noise_skymap.Reset()
        for binx in range(0,hist_noise_skymap.GetNbinsX()):
            for biny in range(0,hist_noise_skymap.GetNbinsY()):
                if hist_mean_skymap.GetBinContent(binx+1,biny+1)<1.: continue
                error = pow(hist_mean_skymap.GetBinContent(binx+1,biny+1),0.5)
                random_number = np.random.normal(loc=0, scale=error)
                hist_noise_skymap.SetBinContent(binx+1,biny+1,random_number)
        #hist_mean_skymap = CommonPlotFunctions.Smooth2DMap(hist_mean_skymap,smooth_size_spectroscopy,False)
        #hist_noise_skymap = CommonPlotFunctions.Smooth2DMap(hist_noise_skymap,smooth_size_spectroscopy,False)
        for binx in range(0,hist_noise_skymap.GetNbinsX()):
            for biny in range(0,hist_noise_skymap.GetNbinsY()):
                if hist_mean_skymap.GetBinContent(binx+1,biny+1)<1.: continue
                error = pow(hist_mean_skymap.GetBinContent(binx+1,biny+1),0.5)
                content = hist_noise_skymap.GetBinContent(binx+1,biny+1)/error
                zscores += [content]

    return zscores

def GetFluxCalibration(energy,elev):

    if energy_bin_cut_low==0:
        return 1.

    # The energy threshold needs to be as low as 100 GeV for this method to work.

    # 3-tel cut
    #str_flux_calibration = ['1.55e+02', '3.94e+00', '1.66e+00', '1.22e+00', '1.03e+00', '1.26e+00', '1.20e+00', '1.13e+00', '1.30e+00', '1.36e+00', '1.37e+00']
    # v490
    str_flux_calibration = ['2.10e+01', '1.59e+00', '1.26e+00', '1.17e+00', '1.07e+00', '1.39e+00', '1.33e+00', '1.27e+00', '1.28e+00', '1.23e+00', '1.37e+00']

    flux_calibration = []
    for string in str_flux_calibration:
        flux_calibration.append(float(string))

    return flux_calibration[energy]

def flux_crab_func(x):
    # TeV^{-1}cm^{-2}s^{-1}
    # Crab https://arxiv.org/pdf/1508.06442.pdf
    return 37.5*pow(10,-12)*pow(x*1./1000.,-2.467-0.16*log(x/1000.))

def PrintSpectralDataForNaima(energy_axis,src_flux,src_flux_err,data_name):
    
    energy_mean_log = [] 
    energy_mean = [] 
    energy_edge_lo = [] 
    energy_edge_hi = [] 
    flux_mean = [] 
    flux_error = []
    ul = []
    for eb in range(0,len(energy_axis)):
        energy_mean_log += [math.log10(energy_axis[eb]/1000.)]
    for eb in range(0,len(energy_axis)):
        energy_log_delta = 0.
        if eb+1<len(energy_axis):
            energy_log_delta = energy_mean_log[eb+1]-energy_mean_log[eb]
        else:
            energy_log_delta = energy_mean_log[eb]-energy_mean_log[eb-1]
        energy_mean += [pow(10,energy_mean_log[eb])]
        energy_edge_lo += [pow(10,energy_mean_log[eb]-0.5*energy_log_delta)]
        energy_edge_hi += [pow(10,energy_mean_log[eb]+0.5*energy_log_delta)]
        flux_mean += [src_flux[eb]/((energy_axis[eb]/1000.)*(energy_axis[eb]/1000.))]
        flux_error += [src_flux_err[eb]/((energy_axis[eb]/1000.)*(energy_axis[eb]/1000.))]
    print ('=======================================================')
    print ('NAIMA flux points')
    print ('data_name = %s'%(data_name))
    for eb in range(0,len(energy_axis)):
        print ('%.2f %.2f %.2f %.2e %.2e %s'%(energy_mean[eb],energy_edge_lo[eb],energy_edge_hi[eb],flux_mean[eb],flux_error[eb],0))
    print ('=======================================================')

def GetHawcDiffusionFluxJ1908():

    energies = [1.19,1.82,3.12,5.52,9.96,18.65,34.17,59.71,103.07,176.38]
    fluxes = [1.95e-11,1.98e-11,2.00e-11,1.57e-11,1.18e-11,7.19e-12,4.70e-12,2.75e-12,2.13e-12,1.38e-12]
    flux_errs = [1.95e-11,1.98e-11,2.00e-11,1.57e-11,1.18e-11,7.19e-12,4.70e-12,2.75e-12,2.13e-12,1.38e-12]
    flux_errs_up = [+0.14e-11,+0.14e-11,+0.13e-11,+0.09e-11,+0.07e-11,+0.55e-12,+0.46e-12,+0.43e-12,+0.44e-12,+0.54e-12]
    flux_errs_low = [-0.15e-11,-0.13e-11,-0.13e-11,-0.09e-11,-0.07e-11,-0.53e-11,-0.45e-12,-0.42e-12,-0.47e-12,-0.54e-12]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,2)
        flux_errs_up[entry] = flux_errs[entry]+flux_errs_up[entry]
        flux_errs_low[entry] = flux_errs[entry]+flux_errs_low[entry]
        flux_errs[entry] = 0.25*fluxes[entry]

    return energies, fluxes, flux_errs

def flux_lhaaso_wcda_j1908_func(x):
    # TeV^{-1}cm^{-2}s^{-1}
    # https://arxiv.org/pdf/2305.17030.pdf
    Flux_N0 = 7.97 
    Gamma_index = 2.42
    return Flux_N0*pow(10,-13)*pow(x*1./3000.,-Gamma_index)

def flux_gamma_cygni_func(x):
    return 1.5*pow(10,-12)*pow(x*1./1000.,-2.37)

def flux_ic443_func(x):
    # IC 443 https://arxiv.org/pdf/0905.3291.pdf
    #return 0.838*pow(10,-12)*pow(x*1./1000.,-2.99)
    # IC 443 https://arxiv.org/pdf/1512.01911.pdf
    return 9.92*pow(10,-12)*pow(x*1./550.,-2.8)
    #return 9.92*pow(10,-13)*pow(x*1./1000.,-2.8)
def flux_ic443_hawc_func(x):
    # 3HWC J0617+224 https://arxiv.org/pdf/2007.08582.pdf
    return 4.5*pow(10,-15)*pow(x*1./7000.,-3.05)

def GetHawcSaraFluxJ1908():

    #energies = [1.53,2.78,4.75,7.14,11.15,18.68,36.15,61.99,108.69,187.51]
    #fluxes = [7.0009e-12,9.5097e-12,8.4629e-12,6.6242e-12,5.6764e-12,4.4924e-12,3.2932e-12,1.5250e-12,9.1235e-13,4.1833e-13]
    #flux_errs = [7.0009e-12,9.5097e-12,8.4629e-12,6.6242e-12,5.6764e-12,4.4924e-12,3.2932e-12,1.5250e-12,9.1235e-13,4.1833e-13]
    #flux_errs_up = [+7.2024e-13,+6.3288e-13,+5.4679e-13,+3.9318e-13,+2.6768e-13,+2.9978e-13,+2.2130e-13,+1.8650e-13,+1.8756e-13,+1.5458e-13]
    #flux_errs_low = [-7.1498e-13,-6.6198e-13,-5.2961e-13,-3.8152e-13,-2.8404e-13,-3.1157e-13,-2.0721e-13,-1.8818e-13,-1.7827e-13,-1.5612e-13]

    energies = [1.53,2.78,4.75,7.14,11.15,18.68,36.15,61.99,108.69,187.51]
    fluxes = [7.0009e-12,9.5097e-12,8.4629e-12,6.6242e-12,5.6764e-12,4.4924e-12,3.2932e-12,1.5250e-12,9.1235e-13,4.1833e-13]
    flux_errs = [7.0009e-12,9.5097e-12,8.4629e-12,6.6242e-12,5.6764e-12,4.4924e-12,3.2932e-12,1.5250e-12,9.1235e-13,4.1833e-13]
    flux_errs_low = [-1.293e-12,-1.425e-12,-1.623e-12,-1.292e-12,-8.837e-13,-4.499e-13,-4.926e-13,-2.418e-13,-1.947e-13,-1.549e-13]
    flux_errs_up = [+8.199e-13,+6.952e-13,+7.313e-13,+5.473e-13,+4.411e-13,+4.389e-13,+3.015e-13,+1.893e-13,+2.074e-13,+1.649e-13]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,2)
        flux_errs_up[entry] = flux_errs[entry]+flux_errs_up[entry]
        flux_errs_low[entry] = flux_errs[entry]+flux_errs_low[entry]
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,2)

    return energies, fluxes, flux_errs

def GetHessFluxJ1908():
    energies = [pow(10.,-0.332),pow(10.,0.022),pow(10.,0.396),pow(10.,0.769),pow(10.,1.124),pow(10.,1.478)]
    fluxes = [pow(10.,-10.981),pow(10.,-10.967),pow(10.,-11.057),pow(10.,-11.169),pow(10.,-11.188),pow(10.,-11.386)]
    flux_errs = [pow(10.,-0.332),pow(10.,0.022),pow(10.,0.396),pow(10.,0.769),pow(10.,1.124),pow(10.,1.478)]
    flux_errs_up = [pow(10.,-10.895),pow(10.,-10.916),pow(10.,-11.003),pow(10.,-11.101),pow(10.,-11.101),pow(10.,-11.264)]
    flux_errs_low = [pow(10.,-11.086),pow(10.,-11.010),pow(10.,-11.126),pow(10.,-11.264),pow(10.,-11.292),pow(10.,-11.556)]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,2)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,2)

    return energies, fluxes, flux_errs

def GetVeritasFluxJ1908():
    energies = [pow(10.,-0.270),pow(10.,-0.126),pow(10.,0.022),pow(10.,0.175),pow(10.,0.323),pow(10.,0.467),pow(10.,0.618),pow(10.,0.776),pow(10.,0.922),pow(10.,1.070),pow(10.,1.219)]
    fluxes = [pow(10.,-11.061),pow(10.,-11.028),pow(10.,-11.036),pow(10.,-11.097),pow(10.,-11.448),pow(10.,-11.166),pow(10.,-11.213),pow(10.,-11.068),pow(10.,-11.209),pow(10.,-11.231),pow(10.,-11.318)]
    flux_errs = [pow(10.,-0.270),pow(10.,-0.126),pow(10.,0.022),pow(10.,0.175),pow(10.,0.323),pow(10.,0.467),pow(10.,0.618),pow(10.,0.776),pow(10.,0.922),pow(10.,1.070),pow(10.,1.219)]
    flux_errs_up = [pow(10.,-10.952),pow(10.,-10.960),pow(10.,-10.974),pow(10.,-11.028),pow(10.,-11.303),pow(10.,-11.083),pow(10.,-11.112),pow(10.,-10.974),pow(10.,-11.083),pow(10.,-11.097),pow(10.,-11.141)]
    flux_errs_low = [pow(10.,-11.245),pow(10.,-11.155),pow(10.,-11.137),pow(10.,-11.195),pow(10.,-11.661),pow(10.,-11.282),pow(10.,-11.339),pow(10.,-11.162),pow(10.,-11.372),pow(10.,-11.408),pow(10.,-11.538)]

    scale_factor = 1./2.03
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]*scale_factor/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,2)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*scale_factor/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,2)

    return energies, fluxes, flux_errs

def GetFermiJordanFluxJ1908():

    energies = [42571.11253606245,85723.52082084052,172617.57055787765,347592.1821687443]
    fluxes = [2.856783157929038e-06,3.89109583469775e-06,5.0680678657082445e-06,9.271213817855382e-06]
    flux_stat_errs = [1.1604625384485099e-06,1.556189798998829e-06,2.2448723890895238e-06,3.4737117958614837e-06]
    flux_syst_errs = [1.135182978267407e-06,7.805371450450492e-07,1.6102184866176414e-06,1.5283362877401339e-06]
    flux_errs = []

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]/1e3
        fluxes[entry] = fluxes[entry]/1e6
        flux_stat_errs[entry] = flux_stat_errs[entry]/1e6
        flux_syst_errs[entry] = flux_syst_errs[entry]/1e6
        flux_errs += [pow(pow(flux_stat_errs[entry],2)+pow(flux_syst_errs[entry],2),0.5)]

    return energies, fluxes, flux_errs

def GetCTASensitivity():

    log_energy_tev = [-1.20,-1.00,-0.80,-0.60,-0.39,-0.19,0.00,0.20,0.40,0.60,0.80,0.99,1.21,1.39,1.59,1.80,2.00,2.20]
    log_flux_e2dnde_erg = [-11.44,-12.03,-12.23,-12.42,-12.57,-12.73,-12.87,-12.98,-13,-13.02,-12.98,-12.98,-12.79,-12.70,-12.47,-12.35,-12.15,-11.95]
    log_energy_tev = np.array(log_energy_tev)
    log_flux_e2dnde_erg = np.array(log_flux_e2dnde_erg)
    energies = pow(10,log_energy_tev)*1000.
    upper_limit = pow(10,log_flux_e2dnde_erg)*0.62

    return energies, upper_limit

def GetVeritasSensitivity():

    log_energy_tev = [-0.863,-0.753,-0.626,-0.448,-0.248,-0.012,0.255,0.507,0.718,0.912,1.117]
    log_flux_e2dnde_erg = [-11.55,-11.72,-11.88,-12.01,-12.15,-12.21,-12.18,-12.14,-12.04,-11.89,-11.78]
    log_energy_tev = np.array(log_energy_tev)
    log_flux_e2dnde_erg = np.array(log_flux_e2dnde_erg)
    energies = pow(10,log_energy_tev)*1000.
    upper_limit = pow(10,log_flux_e2dnde_erg)*0.62

    return energies, upper_limit

def GetTobiasFluxSS433():

    #log_energy_tev = [0.00,0.22,0.43,0.64]
    #log_flux_dnde = [-12.36,-12.56,-12.66,-13.57]
    #log_flux_dnde_hi = [-12.18,-12.46,-12.61,-13.44]
    #log_flux_dnde_lo = [-12.47,-12.70,-12.73,-13.77]
    #log_energy_tev = np.array(log_energy_tev)
    #log_flux_dnde = np.array(log_flux_dnde)
    #log_flux_dnde_hi = np.array(log_flux_dnde_hi)
    #log_flux_dnde_lo = np.array(log_flux_dnde_lo)
    #energies = pow(10,log_energy_tev)*1000.
    #fluxes = pow(10,log_flux_dnde)*pow(energies/1000.,2)
    #flux_hi = pow(10,log_flux_dnde_hi)*pow(energies/1000.,2)
    #flux_lo = pow(10,log_flux_dnde_lo)*pow(energies/1000.,2)
    #flux_errs = 0.5*(flux_hi-flux_lo)

    log_energy_tev = [0.00,0.22,0.43,0.64]
    log_energy_tev = np.array(log_energy_tev)
    energies = pow(10,log_energy_tev)*1000.
    fluxes = [3.57,2.94,1.52,0.36]
    flux_errs = [1.02,0.56,0.31,0.12]
    fluxes = np.array(fluxes)
    flux_errs = np.array(flux_errs)
    fluxes = fluxes*pow(10,-13)*pow(energies/1000.,2)
    flux_errs = flux_errs*pow(10,-13)*pow(energies/1000.,2)

    return energies, fluxes, flux_errs

def GetVeritasTobiasFluxJ1908():

    energy_edges = [794,1580,3160,6310,12600]
    energies = []
    for edge in range(0,len(energy_edges)-1):
        energies += [0.5*(energy_edges[edge]+energy_edges[edge+1])]
    fluxes = [8.96e-12, 8.99e-12, 7.49e-12, 4.41e-12]
    flux_errs = []
    fluxes_imp = []
    #fluxes_imp += [[4.79e-13, 1.23e-12, 2.00e-13, -9.25e-13]]
    #fluxes_imp += [[-5.14e-13, 8.54e-14, 4.01e-13, -3.21e-13]]
    #fluxes_imp += [[5.16e-13, 7.70e-13, 1.16e-12, -5.53e-13]]
    #fluxes_imp += [[-2.54e-13, 2.10e-12, 6.83e-13, -4.54e-13]]
    #fluxes_imp += [[-5.66e-13, -4.65e-13, -5.27e-13, -3.10e-13]]
    fluxes_imp += [[ 3.37038558e-13,  2.32611159e-13,  1.13754728e-14, -1.10503701e-14]]
    fluxes_imp += [[-4.54858909e-13,  1.21982883e-14,  1.70790747e-14, -3.47624119e-15]]
    fluxes_imp += [[ 4.03344751e-13,  1.46921705e-13,  6.16650492e-14, -6.37096588e-15]]
    fluxes_imp += [[-2.39568041e-13,  3.99530818e-13,  3.29603237e-14, -6.06316366e-15]]
    fluxes_imp += [[-4.13799877e-13, -1.02845122e-13, -2.65226278e-14, -3.41802916e-15]]

    for entry in range(0,len(energies)):
        syst_err = 0.
        for imp in range(0,len(fluxes_imp)):
            syst_err += pow(fluxes_imp[imp][entry],2)
        syst_err = pow(syst_err/float(len(fluxes_imp)-1),0.5)
        flux_errs += [pow(pow(syst_err,2)+pow(0.25*fluxes[entry],2),0.5)]

    print ('Tobias_energies = %s'%(energies))
    print ('Tobias_fluxes = %s'%(fluxes))
    return energies, fluxes, flux_errs

def GetLHAASOFluxJ1908():
    energies = [pow(10.,1.102),pow(10.,1.302),pow(10.,1.498),pow(10.,1.700),pow(10.,1.900),pow(10.,2.099),pow(10.,2.299),pow(10.,2.498),pow(10.,2.697)]
    fluxes = [pow(10.,-11.033),pow(10.,-10.988),pow(10.,-11.201),pow(10.,-11.324),pow(10.,-11.553),pow(10.,-11.860),pow(10.,-11.921),pow(10.,-12.346),pow(10.,-12.653)]
    flux_errs = [pow(10.,1.102),pow(10.,1.302),pow(10.,1.498),pow(10.,1.700),pow(10.,1.900),pow(10.,2.099),pow(10.,2.299),pow(10.,2.498),pow(10.,2.697)]
    flux_errs_up = [pow(10.,-10.966),pow(10.,-10.949),pow(10.,-11.167),pow(10.,-11.296),pow(10.,-11.513),pow(10.,-11.798),pow(10.,-11.854),pow(10.,-12.173),pow(10.,-12.391)]
    flux_errs_low = [pow(10.,-11.094),pow(10.,-11.027),pow(10.,-11.240),pow(10.,-11.368),pow(10.,-11.597),pow(10.,-11.944),pow(10.,-12.022),pow(10.,-12.536),pow(10.,-13.128)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,2)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,2)

    return energies, fluxes, flux_errs

def MakeSpectrum(roi_x,roi_y,roi_r,roi_name,excl_roi_x,excl_roi_y,excl_roi_r):

    imposter_flux_list = []
    imposter_flux_err_list = []
    for imposter in range(0,n_imposters):
        energy_axis, energy_error, imposter_flux, imposter_flux_err = CommonPlotFunctions.GetRegionSpectrum(hist_imposter_flux_skymap[imposter],energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
        imposter_flux_list += [imposter_flux]
        imposter_flux_err_list += [imposter_flux_err]
    energy_axis, energy_error, real_flux, real_flux_stat_err = CommonPlotFunctions.GetRegionSpectrum(hist_real_flux_skymap,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)

    energy_axis, energy_error, real_data, real_data_stat_err = CommonPlotFunctions.GetRegionSpectrum(hist_real_data_skymap,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
    energy_axis, energy_error, real_bkgd, real_bkgd_stat_err = CommonPlotFunctions.GetRegionSpectrum(hist_real_bkgd_skymap,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)

    real_flux_syst_err = []
    for ebin in range(0,len(energy_axis)):
        syst_err = 0.
        for imposter in range(0,n_imposters):
            syst_err += pow(imposter_flux_list[imposter][ebin],2)
        if n_imposters>0:
            syst_err = pow(syst_err/float(n_imposters),0.5)
        real_flux_syst_err += [syst_err]

    vectorize_f_crab = np.vectorize(flux_crab_func)
    xdata_array = []
    for binx in range(0,len(energy_axis)):
        xdata_array += [energy_axis[binx]]
    ydata_crab_ref = pow(np.array(xdata_array)/1e3,2)*vectorize_f_crab(xdata_array)

    if 'Crab' in source_name:
        log_energy = np.linspace(log10(2e2),log10(1.2e4),50)
        xdata = pow(10.,log_energy)
        ydata_crab = pow(xdata/1e3,2)*vectorize_f_crab(xdata)
        calibration_new = []
        for binx in range(0,len(energy_axis)):
            if real_flux[binx]>0.:
                calibration_new += [ydata_crab_ref[binx]/real_flux[binx]]
            else:
                calibration_new += [0.]
        print ('=======================================================================')
        formatted_numbers = ['%0.2e' % num for num in calibration_new]
        print ('new flux_calibration = %s'%(formatted_numbers))
        print ('=======================================================================')

    imposter_data_list = []
    imposter_bkgd_list = []
    for imposter in range(0,n_imposters):
        energy_axis, energy_error, imposter_data, imposter_data_err = CommonPlotFunctions.GetRegionSpectrum(hist_imposter_data_skymap[imposter],energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
        energy_axis, energy_error, imposter_bkgd, imposter_bkgd_err = CommonPlotFunctions.GetRegionSpectrum(hist_imposter_bkgd_skymap[imposter],energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
        imposter_data_list += [imposter_data]
        imposter_bkgd_list += [imposter_bkgd]

    real_bkgd_syst_err = []
    for ebin in range(0,len(energy_axis)):
        syst_err = 0.
        for imposter in range(0,n_imposters):
            syst_err += pow(imposter_data_list[imposter][ebin]-imposter_bkgd_list[imposter][ebin],2)
        if n_imposters>0:
            syst_err = pow(syst_err/float(n_imposters),0.5)
        real_bkgd_syst_err += [syst_err]


    imposter_data_avg = []
    for ebin in range(0,len(energy_axis)):
        avg = 0.
        for imposter in range(0,n_imposters):
            if n_imposters>0:
                avg += imposter_data_list[imposter][ebin]/n_imposters
        imposter_data_avg += [avg]
    imposter_s2b_list = []
    imposter_s2b_err_list = []
    for imposter in range(0,n_imposters):
        imposter_s2b = []
        imposter_s2b_err = []
        for ebin in range(0,len(energy_axis)):
            s2b_lower_bound = 0.
            s2b = 0.
            s2b_err = 0.
            if doImposter and imposter_bkgd_list[imposter][ebin]>0. and imposter_data_avg[ebin]>0.:
                s2b = (imposter_data_list[imposter][ebin]-imposter_bkgd_list[imposter][ebin])/imposter_bkgd_list[imposter][ebin]
                s2b_err = pow(imposter_data_list[imposter][ebin],0.5)/imposter_bkgd_list[imposter][ebin]
            imposter_s2b += [s2b]
            imposter_s2b_err += [s2b_err]
        imposter_s2b_list += [imposter_s2b]
        imposter_s2b_err_list += [imposter_s2b_err]

    real_rel_syst_bias = []
    real_rel_syst_err = []
    for ebin in range(0,len(energy_axis)):
        syst_bias = 0.
        syst_err = 0.
        for imposter in range(0,n_imposters):
            syst_bias += imposter_s2b_list[imposter][ebin]
            syst_err += pow(imposter_s2b_list[imposter][ebin],2)
        if n_imposters>0:
            syst_bias = syst_bias/float(n_imposters)
            syst_err = pow(syst_err/float(n_imposters),0.5)
        if syst_err>0:
            real_rel_syst_bias += [pow(n_imposters,0.5)*syst_bias/syst_err]
            real_rel_syst_err += [syst_err]
        else:
            real_rel_syst_bias += [0.]
            real_rel_syst_err += [0.]

    energy_error = np.array(energy_error)
    real_rel_syst_err = np.array(real_rel_syst_err)
    real_flux_syst_err = np.array(real_flux_syst_err)
    real_flux_stat_err = np.array(real_flux_stat_err)
    #real_flux_stat_err = np.array(real_bkgd_stat_err)/np.array(real_bkgd)*np.array(real_flux)
    real_flux_total_err = np.sqrt(np.square(real_flux_stat_err)+np.square(real_flux_syst_err))
    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    cycol = cycle('rgbcmy')
    for imposter in range(0,n_imposters):
        next_color = next(cycol)
        axbig.errorbar(energy_axis,imposter_s2b_list[imposter],imposter_s2b_err_list[imposter],xerr=energy_error,color=next_color,marker='_',ls='none')
    axbig.bar(energy_axis, 2.*real_rel_syst_err, bottom=0.-real_rel_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
    for ebin in range(0,len(energy_axis)):
        axbig.text(energy_axis[ebin]-0.5*energy_error[ebin], 2.2*real_rel_syst_err[len(real_rel_syst_err)-1], '%0.1f %%'%(real_rel_syst_err[ebin]*100.))
        axbig.text(energy_axis[ebin]-0.5*energy_error[ebin], 2.0*real_rel_syst_err[len(real_rel_syst_err)-1], '%0.1f $\sigma$'%(real_rel_syst_bias[ebin]))
    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('relative error')
    axbig.set_xscale('log')
    plotname = 'ImposterRelError_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    if 'Crab' in source_name:
        axbig = fig.add_subplot()
        axbig.plot(xdata, ydata_crab,'r-',label='1508.06442', zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',zorder=3)

        axbig.set_xlabel('Energy [GeV]')
        axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
        axbig.set_xscale('log')
        axbig.set_yscale('log')
        axbig.legend(loc='best')
        plotname = 'RealSpectrum_%s_r%s'%(roi_name,roi_r)
        fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
        axbig.remove()

    elif 'PSR_J1907_p0602' in source_name:
        log_energy = np.linspace(log10(1e2),log10(1e5),50)
        xdata_ref = pow(10.,log_energy)
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHawcDiffusionFluxJ1908()
        Sara_energies, Sara_fluxes, Sara_flux_errs = GetHawcSaraFluxJ1908()
        HESS_energies, HESS_fluxes, HESS_flux_errs = GetHessFluxJ1908()
        OldV_energies, OldV_fluxes, OldV_flux_errs = GetVeritasFluxJ1908()
        Jordan_energies, Jordan_fluxes, Jordan_flux_errs = GetFermiJordanFluxJ1908()
        Tobias_energies, Tobias_fluxes, Tobias_flux_errs = GetVeritasTobiasFluxJ1908()
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ1908()

        vectorize_f_wcda = np.vectorize(flux_lhaaso_wcda_j1908_func)
        log_energy = np.linspace(log10(1000.),log10(25000.),50)
        xdata = pow(10.,log_energy)
        ydata_wcda = pow(xdata/1e3,2)*vectorize_f_wcda(xdata)

        axbig = fig.add_subplot()
        axbig.errorbar(Jordan_energies,Jordan_fluxes,Jordan_flux_errs,color='g',marker='s',ls='none',label='Fermi-LAT',zorder=1)

        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='r', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='r',marker='_',ls='none',label='VERITAS',zorder=3)

        #axbig.errorbar(HESS_energies,HESS_fluxes,HESS_flux_errs,color='cyan',marker='s',ls='none',label='HESS',zorder=4)

        axbig.errorbar(Sara_energies,Sara_fluxes,Sara_flux_errs,color='purple',marker='o',ls='none',label='HAWC',zorder=5)

        # In the case of WCDA data, the overall systematic uncertainty can be as large as +8% −24% on the flux.
        # https://arxiv.org/pdf/2305.17030.pdf, sec 3.3
        axbig.fill_between(xdata, ydata_wcda-0.24*ydata_wcda, ydata_wcda+0.08*ydata_wcda,color='goldenrod', alpha=0.2, zorder=7)
        axbig.plot(xdata, ydata_wcda,color='goldenrod',label='LHAASO (WCDA)', zorder=8)
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='goldenrod',marker='^',ls='none',label='LHAASO (KM2A)',zorder=9)

        axbig.set_xlabel('Energy [GeV]')
        axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
        axbig.set_xscale('log')
        axbig.set_yscale('log')
        axbig.legend(loc='best')
        plotname = 'RealSpectrum_%s_r%s'%(roi_name,roi_r)
        fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
        axbig.remove()

        axbig = fig.add_subplot()
        axbig.errorbar(Jordan_energies,Jordan_fluxes,Jordan_flux_errs,color='gray',marker='s',ls='none',label='Fermi-LAT',zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='r', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='r',marker='_',ls='none',label='VERITAS (Matrix method)',zorder=3)
        axbig.errorbar(Sara_energies,Sara_fluxes,Sara_flux_errs,color='gray',marker='s',ls='none',label='HAWC',zorder=5)
        axbig.errorbar(Tobias_energies,Tobias_fluxes,Tobias_flux_errs,color='blue',marker='s',ls='none',label='VERITAS (Gammapy-3D)',zorder=1)
        axbig.fill_between(xdata, ydata_wcda-0.24*ydata_wcda, ydata_wcda+0.08*ydata_wcda,color='gray', alpha=0.2, zorder=7)
        axbig.plot(xdata, ydata_wcda,color='gray',label='LHAASO (WCDA)', zorder=8)
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='gray',marker='s',ls='none',label='LHAASO (KM2A)',zorder=9)

        axbig.set_xlabel('Energy [GeV]')
        axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
        axbig.set_xscale('log')
        axbig.set_yscale('log')
        axbig.legend(loc='best')
        plotname = 'TobiasSpectrum_%s_r%s'%(roi_name,roi_r)
        fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
        axbig.remove()


        PrintSpectralDataForNaima(Tobias_energies,Tobias_fluxes,Tobias_flux_errs,'Tobias')
        PrintSpectralDataForNaima(HESS_energies,HESS_fluxes,HESS_flux_errs,'HESS')
        PrintSpectralDataForNaima(Jordan_energies,Jordan_fluxes,Jordan_flux_errs,'Fermi')
        PrintSpectralDataForNaima(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,'LHAASO')
        PrintSpectralDataForNaima(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,'HAWC')
        PrintSpectralDataForNaima(Sara_energies,Sara_fluxes,Sara_flux_errs,'Sara')
        PrintSpectralDataForNaima(energy_axis,real_flux,real_flux_total_err,'VERITAS')

    elif 'SNR_G189_p03' in source_name:

        log_energy = np.linspace(log10(1e2),log10(1e4),50)
        xdata = pow(10.,log_energy)
        vectorize_f_veritas_paper = np.vectorize(flux_ic443_func)
        ydata_veritas_paper = pow(xdata/1e3,2)*vectorize_f_veritas_paper(xdata)
        vectorize_f_hawc = np.vectorize(flux_ic443_hawc_func)
        ydata_hawc = pow(xdata/1e3,2)*vectorize_f_hawc(xdata)

        vts_energies, vts_ul = GetVeritasSensitivity()

        axbig = fig.add_subplot()
        axbig.plot(vts_energies,vts_ul,color='orange',label='VERITAS sensitivity')
        axbig.plot(xdata, ydata_veritas_paper,'r-',label='VERITAS (0905.3291)',zorder=1)
        axbig.plot(xdata, ydata_hawc,'g-',label='HAWC (2007.08582)',zorder=2)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='k',marker='_',ls='none',label='Data')

        axbig.set_xlabel('Energy [GeV]')
        axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
        axbig.set_xscale('log')
        axbig.set_yscale('log')
        axbig.legend(loc='best')
        plotname = 'RealSpectrum_%s_r%s'%(roi_name,roi_r)
        fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
        axbig.remove()

    elif 'PSR_J2021_p4026' in source_name:

        log_energy = np.linspace(log10(1e2),log10(1e4),50)
        xdata = pow(10.,log_energy)
        vectorize_f_veritas_paper = np.vectorize(flux_gamma_cygni_func)
        ydata_veritas_paper = pow(xdata/1e3,2)*vectorize_f_veritas_paper(xdata)

        vts_energies, vts_ul = GetVeritasSensitivity()

        axbig = fig.add_subplot()
        axbig.plot(vts_energies,vts_ul,color='orange',label='VERITAS sensitivity')
        axbig.plot(xdata, ydata_veritas_paper,'r-',label='VERITAS (1305.6508)',zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='k',marker='_',ls='none',label='Data')

        axbig.set_xlabel('Energy [GeV]')
        axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
        axbig.set_xscale('log')
        axbig.set_yscale('log')
        axbig.legend(loc='best')
        plotname = 'RealSpectrum_%s_r%s'%(roi_name,roi_r)
        fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
        axbig.remove()

    elif 'SS433' in source_name:
        axbig = fig.add_subplot()

        vts_energies, vts_ul = GetVeritasSensitivity()
        cta_energies, cta_ul = GetCTASensitivity()
        Tobias_energies, Tobias_fluxes, Tobias_flux_errs = GetTobiasFluxSS433()

        axbig.plot(vts_energies,vts_ul,color='orange',label='50-h VERITAS sensitivity')
        #axbig.plot(cta_energies,cta_ul,color='green',label='50-h CTA-south sensitivity')
        axbig.errorbar(Tobias_energies,Tobias_fluxes,Tobias_flux_errs,color='red',marker='s',ls='none',label='VERITAS (Gammapy-3D)')
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (Matrix method)')
        PrintSpectralDataForNaima(energy_axis,real_flux,real_flux_total_err,'VERITAS')

        axbig.set_xlabel('Energy [GeV]')
        axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
        axbig.set_xscale('log')
        axbig.set_yscale('log')
        axbig.legend(loc='best')
        plotname = 'RealSpectrum_%s_r%s'%(roi_name,roi_r)
        fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
        axbig.remove()

    else:
        axbig = fig.add_subplot()

        vts_energies, vts_ul = GetVeritasSensitivity()

        axbig.plot(vts_energies,vts_ul,color='orange',label='50-h VERITAS sensitivity')
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='k',marker='_',ls='none',label='Data')
        PrintSpectralDataForNaima(energy_axis,real_flux,real_flux_total_err,'VERITAS')

        axbig.set_xlabel('Energy [GeV]')
        axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
        axbig.set_xscale('log')
        axbig.set_yscale('log')
        axbig.legend(loc='best')
        plotname = 'RealSpectrum_%s_r%s'%(roi_name,roi_r)
        fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
        axbig.remove()

    energy_mean_log = [] 
    energy_mean = [] 
    energy_edge_lo = [] 
    energy_edge_hi = [] 
    flux_mean = [] 
    flux_error = []
    ul = []
    for eb in range(0,len(energy_axis)):
        energy_mean_log += [math.log10(energy_axis[eb]/1000.)]
    for eb in range(0,len(energy_axis)):
        energy_log_delta = 0.
        if eb+1<len(energy_axis):
            energy_log_delta = energy_mean_log[eb+1]-energy_mean_log[eb]
        else:
            energy_log_delta = energy_mean_log[eb]-energy_mean_log[eb-1]
        energy_mean += [pow(10,energy_mean_log[eb])]
        energy_edge_lo += [pow(10,energy_mean_log[eb]-0.5*energy_log_delta)]
        energy_edge_hi += [pow(10,energy_mean_log[eb]+0.5*energy_log_delta)]

    zscore = []
    total_data = 0.
    total_bkgd = 0.
    total_bkgd_syst_err = 0.
    for eb in range(0,len(energy_axis)):
        real_bkgd_stat_err[eb] = pow(real_data[eb],0.5)
        if real_data[eb]==0.:
            zscore += [0.]
        else:
            zscore += [(real_data[eb]-real_bkgd[eb])/pow(pow(real_bkgd_stat_err[eb],2)+pow(real_bkgd_syst_err[eb],2),0.5)]
        total_data += real_data[eb]
        total_bkgd += real_bkgd[eb]
        total_bkgd_syst_err += pow(real_bkgd_syst_err[eb],2)
        print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print ('Energy                  = %0.2f-%0.2f TeV'%(energy_edge_lo[eb],energy_edge_hi[eb]))
        print ('significance in RoI     = %0.2f'%(zscore[eb]))
        print ('total count in RoI      = %s'%(real_data[eb]))
        print ('background count in RoI = %0.2f +/- %0.2f (stat error) +/- %0.2f (syst error)'%(real_bkgd[eb],real_bkgd_stat_err[eb],real_bkgd_syst_err[eb]))
        print ('flux in RoI             = %0.2e +/- %0.2e +/- %0.2e TeV/cm2/s'%(real_flux[eb],real_flux_stat_err[eb],real_flux_syst_err[eb]))
    total_bkgd_syst_err = pow(total_bkgd_syst_err,0.5)
    zscore_total = (total_data-total_bkgd)/pow(pow(total_data,1)+pow(total_bkgd_syst_err,2),0.5)
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('total significance in RoI     = %0.2f'%(zscore_total))
    zscore = np.array(zscore)
    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.step(energy_axis, zscore, where='mid', color='b')
    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('Significance')
    axbig.set_xscale('log')
    axbig.legend(loc='best')
    plotname = 'ZscoreSpectrum_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()



def MakeFluxMap(flux_map, data_map, bkgd_map, norm_map, elev_map):

    #norm_map_smooth = CommonPlotFunctions.Smooth2DMap(norm_map,0.3,True)
    skymap_bin_size_x = data_map[0].GetXaxis().GetBinCenter(2)-data_map[0].GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = data_map[0].GetYaxis().GetBinCenter(2)-data_map[0].GetYaxis().GetBinCenter(1)
    calibration_radius = CommonPlotFunctions.calibration_radius
    for ebin in range(0,len(energy_bin)-1):
        flux_map[ebin].Reset()
        #data_map_smooth = CommonPlotFunctions.Smooth2DMap(data_map[ebin],smooth_size_spectroscopy,True)
        #bkgd_map_smooth = CommonPlotFunctions.Smooth2DMap(bkgd_map[ebin],smooth_size_spectroscopy,True)
        #norm_map_smooth = CommonPlotFunctions.Smooth2DMap(bkgd_map[ebin],smooth_size_spectroscopy,True)
        norm_content_max = norm_map.GetMaximum()
        for binx in range(0,bkgd_map[ebin].GetNbinsX()):
            for biny in range(0,bkgd_map[ebin].GetNbinsY()):
                data_content = data_map[ebin].GetBinContent(binx+1,biny+1)
                data_error = data_map[ebin].GetBinError(binx+1,biny+1)
                bkgd_content = bkgd_map[ebin].GetBinContent(binx+1,biny+1)
                norm_content = norm_map.GetBinContent(binx+1,biny+1)
                binx_center = data_map[ebin].GetXaxis().GetBinCenter(binx+1)
                biny_center = data_map[ebin].GetYaxis().GetBinCenter(biny+1)
                new_binx = elev_map.GetXaxis().FindBin(binx_center)
                new_biny = elev_map.GetYaxis().FindBin(biny_center)
                elev_content = elev_map.GetBinContent(new_binx,new_biny)

                #if data_error*data_error<10.: continue
                if norm_content==0.: continue
                if norm_content/norm_content_max<0.1: continue

                correction = GetFluxCalibration(ebin,elev_content)

                norm_weight = 1.
                norm_ratio = norm_content/norm_content_max
                norm_weight = correction*1./(1.+np.exp(-(norm_ratio-0.3)/0.1))

                #correction = correction*norm_weight
                #stat_data_err = pow(max(data_content,0.),0.5)
                #flux_stat_err = max(stat_data_err,1.)/norm_content*correction*pow(energy_bin[ebin]/1e3,2)
                #flux_content = (data_content-bkgd_content)/norm_content*correction*pow(energy_bin[ebin]/1e3,2)
                #flux_map[ebin].SetBinContent(binx+1,biny+1,flux_content)
                #flux_map[ebin].SetBinError(binx+1,biny+1,flux_stat_err)

                delta_E = (energy_bin[ebin+1]-energy_bin[ebin])/(3.12e+00)
                local_areatime = hist_areatime_skymap.GetBinContent(binx+1,biny+1)
                if local_areatime<=0.: continue
                flux_content = (data_content-bkgd_content)/local_areatime*pow(energy_bin[ebin]/1e3,2)/(100.*100.*3600.)/delta_E
                stat_data_err = pow(max(bkgd_content,0.),0.5)
                flux_stat_err = max(stat_data_err,1.)/local_areatime*pow(energy_bin[ebin]/1e3,2)/(100.*100.*3600.)/delta_E
                flux_map[ebin].SetBinContent(binx+1,biny+1,flux_content*norm_weight)
                flux_map[ebin].SetBinError(binx+1,biny+1,flux_stat_err*norm_weight)


# Our function to fit is going to be a sum of two-dimensional Gaussians
def gaussian(x, y, x0, y0, sigma, A):
    #return A * np.exp( -((x-x0)/(2.*sigma))**2 -((y-y0)/(2.*sigma))**2)
    return A * np.exp(-((x-x0)**2+(y-y0)**2)/(2*sigma*sigma))/(2*np.pi*sigma*sigma)
# https://scipython.com/blog/non-linear-least-squares-fitting-of-a-two-dimensional-data/
# This is the callable that is passed to curve_fit. M is a (2,N) array
# where N is the total number of data points in Z, which will be ravelled
# to one dimension.
def _gaussian(M, *args):
    x, y = M
    arr = np.zeros(x.shape)
    for i in range(len(args)//4):
       arr += gaussian(x, y, *args[i*4:i*4+4])
    return arr

def fit_2d_model(hist_map_data, hist_map_bkgd, src_x, src_y):

    nbins_x = hist_map_data.GetNbinsX()
    nbins_y = hist_map_data.GetNbinsY()
    lon_min = MapEdge_left
    lon_max = MapEdge_right
    lat_min = MapEdge_lower
    lat_max = MapEdge_upper
    x_axis = np.linspace(lon_min,lon_max,nbins_x)
    y_axis = np.linspace(lat_min,lat_max,nbins_y)
    X_grid, Y_grid = np.meshgrid(x_axis, y_axis)
    # We need to ravel the meshgrids of X, Y points to a pair of 1-D arrays.
    XY_stack = np.vstack((X_grid.ravel(), Y_grid.ravel()))

    image_data = np.zeros((hist_map_data.GetNbinsX(),hist_map_data.GetNbinsY()))
    image_error = np.zeros((hist_map_data.GetNbinsX(),hist_map_data.GetNbinsY()))
    for binx in range (0,hist_map_data.GetNbinsX()):
        for biny in range (0,hist_map_data.GetNbinsY()):
            image_data[biny,binx] = hist_map_data.GetBinContent(binx+1,biny+1) - hist_map_bkgd.GetBinContent(binx+1,biny+1)
            error = pow(abs(hist_map_data.GetBinContent(binx+1,biny+1)),0.5)
            image_error[biny,binx] = max(1.,pow(error,0.5))

    #print ('set initial avlues and bounds')
    initial_prms = []
    bound_upper_prms = []
    bound_lower_prms = []
    lon = src_x
    lat = src_y
    sigma = 0.03807
    initial_prms += [(lon,lat,sigma,10.)]
    centroid_range = 0.5
    bound_lower_prms += [(lon-centroid_range,lat-centroid_range,sigma+0.0,0.)]
    bound_upper_prms += [(lon+centroid_range,lat+centroid_range,sigma+2.0,1e10)]
    # Flatten the initial guess parameter list.
    p0 = [p for prms in initial_prms for p in prms]
    p0_lower = [p for prms in bound_lower_prms for p in prms]
    p0_upper = [p for prms in bound_upper_prms for p in prms]
    print ('p0 = %s'%(p0))

    popt, pcov = curve_fit(_gaussian, XY_stack, image_data.ravel(), p0, sigma=image_error.ravel(), absolute_sigma=True, bounds=(p0_lower,p0_upper))
    fit_src_x = popt[0*4+0]
    fit_src_x_err = pow(pcov[0*4+0][0*4+0],0.5)
    print ('fit_src_x = %0.3f +/- %0.3f'%(fit_src_x,fit_src_x_err))
    fit_src_y = popt[0*4+1]
    fit_src_y_err = pow(pcov[0*4+1][0*4+1],0.5)
    print ('fit_src_y = %0.3f +/- %0.3f'%(fit_src_y,fit_src_y_err))
    fit_src_sigma = popt[0*4+2]
    fit_src_sigma_err = pow(pcov[0*4+2][0*4+2],0.5)
    print ('fit_src_sigma = %0.3f +/- %0.3f'%(fit_src_sigma,fit_src_sigma_err))
    true_src_sigma = pow(fit_src_sigma*fit_src_sigma-pow(CommonPlotFunctions.smooth_size_spectroscopy,2),0.5)
    print ('true_src_sigma = %0.3f +/- %0.3f'%(true_src_sigma,fit_src_sigma_err))
    fit_src_A = popt[0*4+3]
    print ('fit_src_A = %0.1e'%(fit_src_A))

    distance_to_psr = pow(pow(fit_src_x-src_x,2)+pow(fit_src_y-src_y,2),0.5)
    distance_to_psr_err = pow(pow(fit_src_x_err,2)+pow(fit_src_y_err,2),0.5)
    print ('distance_to_psr = %0.3f +/- %0.3f'%(distance_to_psr,fit_src_x_err))

    profile_fit = _gaussian(XY_stack, *popt)
    residual = image_data.ravel() - profile_fit
    chisq = np.sum((residual/image_error.ravel())**2)
    dof = len(image_data.ravel())-4
    print ('chisq/dof = %0.3f'%(chisq/dof))


def MakeExtensionProfile(roi_x,roi_y,roi_r,fit_profile,roi_name,real_map,imposter_maps,erange_tag):

    if real_map.Integral()==0.:
        return

    plot_radius = 0.5*(MapEdge_upper-MapEdge_lower)
    real_profile, real_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension(real_map,roi_x,roi_y,plot_radius)
    imposter_profile_list = []
    imposter_profile_err_list = []
    for imposter in range(0,n_imposters):
        imposter_profile, imposter_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension(imposter_maps[imposter],roi_x,roi_y,plot_radius)
        imposter_profile_list += [imposter_profile]
        imposter_profile_err_list += [imposter_profile_stat_err]

    real_profile_syst_err = []
    for ubin in range(0,len(theta2)):
        syst_err = 0.
        if n_imposters>0:
            for imposter in range(0,n_imposters):
                syst_err += pow(imposter_profile_list[imposter][ubin],2)
            syst_err = pow(syst_err/float(n_imposters),0.5)
        real_profile_syst_err += [syst_err]

    end_of_array = False
    while not end_of_array:
        for ubin in range(0,len(theta2)):
            delete_entry = False
            if real_profile[ubin]==0.:
                delete_entry = True
            if delete_entry:
                del theta2[ubin]
                del theta2_err[ubin]
                del real_profile[ubin]
                del real_profile_stat_err[ubin]
                del real_profile_syst_err[ubin]
                for imposter in range(0,n_imposters):
                    del imposter_profile_list[imposter][ubin]
                    del imposter_profile_err_list[imposter][ubin]
                break
            if ubin==len(theta2)-1:
                end_of_array = True

    theta2 = np.array(theta2)
    theta2_err = np.array(theta2_err)
    real_profile = np.array(real_profile)
    real_profile_stat_err = np.array(real_profile_stat_err)
    real_profile_syst_err = np.array(real_profile_syst_err)
    real_profile_total_err = np.sqrt(np.square(real_profile_stat_err)+np.square(real_profile_syst_err))
    #real_profile_total_err = []
    #for ubin in range(0,len(theta2)):
    #    stat_err = real_profile_stat_err[ubin]
    #    syst_err = real_profile_syst_err[ubin]
    #    real_profile_total_err += [pow(stat_err*stat_err+syst_err*syst_err,0.5)]
    #real_profile_total_err = np.array(real_profile_total_err)

    profile_sum = 0.
    for ubin in range(0,len(theta2)):
        profile_sum += real_profile[ubin]
    if fit_profile==1:
        start = (profile_sum, 0.5)
        popt, pcov = curve_fit(diffusion_func,theta2,real_profile,p0=start,sigma=real_profile_total_err,absolute_sigma=True,bounds=((0, 0.01), (np.inf, np.inf)))
        profile_fit = diffusion_func(theta2, *popt)
        residual = real_profile - profile_fit
        chisq = np.sum((residual/real_profile_stat_err)**2)
        dof = len(theta2)-2
        print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print (erange_tag)
        print ('diffusion flux = %0.2E +/- %0.2E'%(popt[0],pow(pcov[0][0],0.5)))
        print ('diffusion radius = %0.2f +/- %0.2f deg (chi2/dof = %0.2f)'%(popt[1],pow(pcov[1][1],0.5),chisq/dof))
    elif fit_profile==2:
        start = (profile_sum, 0.5)
        popt, pcov = curve_fit(gauss_func,theta2,real_profile,p0=start,sigma=real_profile_stat_err,absolute_sigma=True,bounds=((0, 0.01), (np.inf, np.inf)))
        profile_fit = gauss_func(theta2, *popt)
        residual = real_profile - profile_fit
        chisq = np.sum((residual/real_profile_stat_err)**2)
        dof = len(theta2)-2
        print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print (erange_tag)
        print ('gaussian flux = %0.2E +/- %0.2E'%(popt[0],pow(pcov[0][0],0.5)))
        print ('gaussian radius = %0.2f +/- %0.2f deg (chi2/dof = %0.2f)'%(popt[1],pow(pcov[1][1],0.5),chisq/dof))

    delta_theta = theta2[1]-theta2[0]
    baseline_xaxis = []
    for i in range(-1,len(theta2)+1):
        baseline_xaxis += [theta2[0]+delta_theta*i]
    baseline_yaxis = [0. for i in range(0,len(baseline_xaxis))]
    fig.clf()
    fig.set_figheight(figsize_y*0.8)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.plot(baseline_xaxis, baseline_yaxis, color='b', ls='dashed')
    axbig.errorbar(theta2,real_profile,real_profile_total_err,color='k',marker='s',ls='none',label='ON data')
    axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err+real_profile, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    #if fit_profile!=0:
    #    axbig.plot(theta2,diffusion_func(theta2,*popt),color='r')
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
    axbig.set_xlabel('angular distance from source [degree]')
    axbig.legend(loc='best')
    plotname = 'ProfileVsTheta2_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s_%s.png"%(plotname,erange_tag,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    fig.set_figheight(figsize_y*0.8)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    cycol = cycle('rgbcmy')
    axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    for imposter in range(0,n_imposters):
        next_color = next(cycol)
        axbig.errorbar(theta2,imposter_profile_list[imposter],imposter_profile_err_list[imposter],color=next_color,marker='s',ls='none',label='Mimic #%s'%(imposter+1))
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
    axbig.set_xlabel('angular distance from source [degree]')
    axbig.legend(loc='best')
    plotname = 'ProfileVsTheta2_Imposter_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s_%s.png"%(plotname,erange_tag,plot_tag),bbox_inches='tight')
    axbig.remove()


folder_path = CommonPlotFunctions.folder_path

epoch_idx = 0
SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_%s_G0_X0_Y0.root"%(folder_path,source_name,list_epoch[0],isON)
if os.path.exists(SourceFilePath):
    epoch_idx = 0
else:
    epoch_idx = 1
InputFile = ROOT.TFile("/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_%s_G0_X0_Y0.root"%(folder_path,source_name,list_epoch[epoch_idx],isON))
HistName = "Hist_OnData_SR_Skymap_Sum_ErecS100to200"

nbins_x = InputFile.Get(HistName).GetNbinsX()
nbins_y = InputFile.Get(HistName).GetNbinsY()
MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
MapCenter_x = (MapEdge_left+MapEdge_right)/2.
MapCenter_y = (MapEdge_lower+MapEdge_upper)/2.
binsize_x = (MapEdge_right-MapEdge_left)/nbins_x
binsize_y = (MapEdge_upper-MapEdge_lower)/nbins_y

nbins_x = new_nbins_x
nbins_y = new_nbins_y
MapEdge_right = MapCenter_x+0.5*nbins_x*binsize_x
MapEdge_left = MapCenter_x-0.5*nbins_x*binsize_x
MapEdge_upper = MapCenter_y+0.5*nbins_y*binsize_y
MapEdge_lower = MapCenter_y-0.5*nbins_y*binsize_y

excl_region_x = [MapCenter_x]
excl_region_y = [MapCenter_y]
excl_region_r = [0.0]
region_x = [MapCenter_x]
region_y = [MapCenter_y]
region_r = [0.2]
region_name = 'Center'
do_fit = 0
if 'Crab' in source_name:
    region_x = [MapCenter_x]
    region_y = [MapCenter_y]
    region_r = [CommonPlotFunctions.calibration_radius]
    #region_r = [1.0]
    region_name = 'Center'
elif 'Geminga' in source_name:
    region_x = [MapCenter_x]
    region_y = [MapCenter_y]
    region_r = [1.5]
    region_name = 'Center'
elif 'PSR_J1907_p0602' in source_name:

    #3HWC J1908+063, 287.05, 6.39 
    region_x = [287.05]
    region_y = [6.39]
    region_r = [1.2]
    region_name = '3HWC'

    #region_x = [286.56]
    #region_y = [7.20]
    #region_r = [0.3]
    #region_name = 'LHAASO'

elif 'SS433' in source_name:

    #SS 433 SNR
    #region_x = [288.0833333]
    #region_y = [4.9166667]
    #region_r = [0.3]
    #region_name = 'SS433'

    #SS 433 e1
    region_x = [288.404]
    region_y = [4.930]
    region_r = [0.3]
    region_name = 'SS433e1'
    #region_x = [288.35,288.50,288.65,288.8]
    #region_y = [4.93,4.92,4.93,4.94]
    #region_r = [0.1,0.1,0.1,0.1]
    #region_name = 'SS433e1'

    #SS 433 w1
    #region_x = [287.654]
    #region_y = [5.037]
    #region_r = [0.3]
    #region_name = 'SS433w1'

elif 'PSR_J1856_p0245' in source_name:

    region_x = [284.2958333]
    region_y = [2.6666667]
    region_r = [1.0]
    region_name = 'HESS'

elif 'PSR_J2021_p4026' in source_name:
    region_x = [305.0200000]
    region_y = [40.7572222]
    region_r = [0.3]
    region_name = 'Center'
elif '2HWC_J1953_p294' in source_name:
    # G067.6+00.9
    region_x = [299.44]
    region_y = [30.88]
    region_r = [1.0]
    region_name = 'Center'


InputFile.Close()

hist_areatime_skymap = ROOT.TH2D("hist_areatime_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_elev_skymap = ROOT.TH2D("hist_elev_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_azim_skymap = ROOT.TH2D("hist_azim_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_nsb_skymap = ROOT.TH2D("hist_nsb_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)

hist_real_flux_skymap_le = ROOT.TH2D("hist_real_flux_skymap_le","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_flux_skymap_he = ROOT.TH2D("hist_real_flux_skymap_he","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_flux_skymap_sum = ROOT.TH2D("hist_real_flux_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_data_skymap_le = ROOT.TH2D("hist_real_data_skymap_le","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_bkgd_skymap_le = ROOT.TH2D("hist_real_bkgd_skymap_le","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_data_skymap_he = ROOT.TH2D("hist_real_data_skymap_he","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_bkgd_skymap_he = ROOT.TH2D("hist_real_bkgd_skymap_he","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_data_skymap_sum = ROOT.TH2D("hist_real_data_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_bkgd_skymap_sum = ROOT.TH2D("hist_real_bkgd_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_norm_skymap_sum = ROOT.TH2D("hist_real_norm_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_diff_skymap_sum = ROOT.TH2D("hist_real_diff_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_diff_skymap_le = ROOT.TH2D("hist_real_diff_skymap_le","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_diff_skymap_he = ROOT.TH2D("hist_real_diff_skymap_he","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_data_msclw_sum = ROOT.TH2D("hist_real_data_msclw_sum","",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower_fine,MSCL_plot_upper_fine,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower_fine,MSCW_plot_upper_fine)
hist_real_bkgd_msclw_sum = ROOT.TH2D("hist_real_bkgd_msclw_sum","",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower_fine,MSCL_plot_upper_fine,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower_fine,MSCW_plot_upper_fine)
hist_real_diff_msclw_sum = ROOT.TH2D("hist_real_diff_msclw_sum","",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower_fine,MSCL_plot_upper_fine,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower_fine,MSCW_plot_upper_fine)
hist_real_sign_msclw_sum = ROOT.TH2D("hist_real_sign_msclw_sum","",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower_fine,MSCL_plot_upper_fine,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower_fine,MSCW_plot_upper_fine)
hist_real_significance_skymap_sum = ROOT.TH2D("hist_real_significance_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_significance_skymap_le = ROOT.TH2D("hist_real_significance_skymap_le","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_significance_skymap_he = ROOT.TH2D("hist_real_significance_skymap_he","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_flux_skymap = []
hist_real_data_skymap = []
hist_real_bkgd_skymap = []
hist_real_diff_skymap = []
hist_real_data_msclw = []
hist_real_bkgd_msclw = []
for ebin in range(0,len(energy_bin)-1):
    hist_real_flux_skymap += [ROOT.TH2D("hist_real_flux_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_data_skymap += [ROOT.TH2D("hist_real_data_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_bkgd_skymap += [ROOT.TH2D("hist_real_bkgd_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_diff_skymap += [ROOT.TH2D("hist_real_diff_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_data_msclw += [ROOT.TH2D("hist_real_data_msclw_E%s"%(ebin),"",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower_fine,MSCL_plot_upper_fine,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower_fine,MSCW_plot_upper_fine)]
    hist_real_bkgd_msclw += [ROOT.TH2D("hist_real_bkgd_msclw_E%s"%(ebin),"",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower_fine,MSCL_plot_upper_fine,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower_fine,MSCW_plot_upper_fine)]

hist_imposter_data_skymap_sum = []
hist_imposter_bkgd_skymap_sum = []
hist_imposter_diff_skymap_sum = []
hist_imposter_flux_skymap_sum = []
hist_imposter_flux_skymap_le = []
hist_imposter_flux_skymap_he = []
hist_imposter_significance_skymap_sum = []
hist_imposter_data_skymap = []
hist_imposter_bkgd_skymap = []
hist_imposter_diff_skymap = []
hist_imposter_flux_skymap = []
for imposter in range(0,n_imposters):
    hist_imposter_significance_skymap_sum += [ROOT.TH2D("hist_imposter_significance_skymap_sum_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap_sum += [ROOT.TH2D("hist_imposter_data_skymap_sum_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_bkgd_skymap_sum += [ROOT.TH2D("hist_imposter_bkgd_skymap_sum_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_diff_skymap_sum += [ROOT.TH2D("hist_imposter_diff_skymap_sum_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_flux_skymap_sum += [ROOT.TH2D("hist_imposter_flux_skymap_sum_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_flux_skymap_le += [ROOT.TH2D("hist_imposter_flux_skymap_le_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_flux_skymap_he += [ROOT.TH2D("hist_imposter_flux_skymap_he_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap_sublist = []
    hist_imposter_bkgd_skymap_sublist = []
    hist_imposter_diff_skymap_sublist = []
    hist_imposter_flux_skymap_sublist = []
    for ebin in range(0,len(energy_bin)-1):
        hist_imposter_data_skymap_sublist += [ROOT.TH2D("hist_imposter_data_skymap_I%s_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_bkgd_skymap_sublist += [ROOT.TH2D("hist_imposter_bkgd_skymap_I%s_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_diff_skymap_sublist += [ROOT.TH2D("hist_imposter_diff_skymap_I%s_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_flux_skymap_sublist += [ROOT.TH2D("hist_imposter_flux_skymap_I%s_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap += [hist_imposter_data_skymap_sublist]
    hist_imposter_bkgd_skymap += [hist_imposter_bkgd_skymap_sublist]
    hist_imposter_diff_skymap += [hist_imposter_diff_skymap_sublist]
    hist_imposter_flux_skymap += [hist_imposter_flux_skymap_sublist]


effective_area = ROOT.std.vector("double")(20)

n_samples = 0.
for xoff_idx in range(0,n_xoff_bins):
    for yoff_idx in range(0,n_yoff_bins):
        for epoch in list_epoch:
            n_groups = 0
            file_exists = True
            while file_exists:
                SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_%s_G%d_X%d_Y%d.root"%(folder_path,source_name,epoch,isON,n_groups,xoff_idx,yoff_idx)
                print ('Read file: %s'%(SourceFilePath))
                if os.path.exists(SourceFilePath):
                    n_groups += 1
                    print ('file exists.')
                else:
                    file_exists = False
                    print ('file does not exist.')
            
            for group in range(0,n_groups):
                n_samples += 1.
                InputFile = ROOT.TFile("/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_%s_G%d_X%d_Y%d.root"%(folder_path,source_name,epoch,isON,group,xoff_idx,yoff_idx))
                InfoTree = InputFile.Get("InfoTree")
                InfoTree.SetBranchAddress('effective_area',ROOT.AddressOf(effective_area))
                InfoTree.GetEntry(0)
                data_expo = InfoTree.exposure_hours
                if xoff_idx==0 and yoff_idx==0:
                    total_data_expo += data_expo
                HistName = "Hist_Data_AreaTime_Skymap"
                FillSkyMapHistogram(InputFile.Get(HistName),hist_areatime_skymap)
                HistName = "Hist_Data_Elev_Skymap"
                FillSkyMapHistogram(InputFile.Get(HistName),hist_elev_skymap)
                HistName = "Hist_Data_Azim_Skymap"
                FillSkyMapHistogram(InputFile.Get(HistName),hist_azim_skymap)
                HistName = "Hist_Data_NSB_Skymap"
                FillSkyMapHistogram(InputFile.Get(HistName),hist_nsb_skymap)
                for ebin in range(0,len(energy_bin)-1):
                    if energy_bin_cut_low>0:
                        if effective_area[ebin] < effective_area_cut: continue
                    HistName = "Hist_OnData_SR_Skymap_Sum_ErecS%sto%s"%(int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                    FillSkyMapHistogram(InputFile.Get(HistName),hist_real_data_skymap[ebin])
                    HistName = "Hist_OnData_CR_Skymap_%s_Sum_ErecS%sto%s"%(analysis_method,int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                    FillSkyMapHistogram(InputFile.Get(HistName),hist_real_bkgd_skymap[ebin])
                    HistName = "Hist_OnData_MSCLW_Fine_Sum_ErecS%sto%s"%(int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                    hist_real_data_msclw[ebin].Add(InputFile.Get(HistName))
                    HistName = "Hist_OnBkgd_MSCLW_Fine_Sum_ErecS%sto%s"%(int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                    hist_real_bkgd_msclw[ebin].Add(InputFile.Get(HistName))
                InputFile.Close()

hist_elev_skymap.Scale(1./n_samples)
hist_azim_skymap.Scale(1./n_samples)
hist_nsb_skymap.Scale(1./n_samples)

if doImposter:
    for imposter in range(0,n_imposters):
        n_imposter_samples = 0.
        for xoff_idx in range(0,n_xoff_bins):
            for yoff_idx in range(0,n_yoff_bins):
                for epoch in list_epoch:
                    n_groups = 0
                    file_exists = True
                    while file_exists:
                        SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_Imposter%s_G%d_X%d_Y%d.root"%(folder_path,source_name,epoch,imposter+1,n_groups,xoff_idx,yoff_idx)
                        print ('Read file: %s'%(SourceFilePath))
                        if os.path.exists(SourceFilePath):
                            n_groups += 1
                            print ('file exists.')
                        else:
                            file_exists = False
                            print ('file does not exist.')
                    
                    for group in range(0,n_groups):
                        n_imposter_samples += 1.
                        InputFile = ROOT.TFile("/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_Imposter%s_G%d_X%d_Y%d.root"%(folder_path,source_name,epoch,imposter+1,group,xoff_idx,yoff_idx))
                        InfoTree = InputFile.Get("InfoTree")
                        InfoTree.SetBranchAddress('effective_area',ROOT.AddressOf(effective_area))
                        InfoTree.GetEntry(0)
                        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
                            if energy_bin_cut_low>0:
                                if effective_area[ebin] < effective_area_cut: continue
                            HistName = "Hist_OnData_SR_Skymap_Sum_ErecS%sto%s"%(int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                            FillSkyMapHistogram(InputFile.Get(HistName),hist_imposter_data_skymap[imposter][ebin])
                            HistName = "Hist_OnData_CR_Skymap_%s_Sum_ErecS%sto%s"%(analysis_method,int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                            FillSkyMapHistogram(InputFile.Get(HistName),hist_imposter_bkgd_skymap[imposter][ebin])
                        InputFile.Close()

# bias correction
if n_imposters>3 and doBiasCorrect:
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):

        mean_bias = 0.
        for imposter in range(0,n_imposters):
            truth = hist_imposter_data_skymap[imposter][ebin].Integral()
            predict = hist_imposter_bkgd_skymap[imposter][ebin].Integral()
            if predict>0.:
                bias = truth/predict
                mean_bias += bias
        mean_bias = mean_bias/float(n_imposters)
        hist_real_bkgd_skymap[ebin].Scale(mean_bias)
        for imposter in range(0,n_imposters):
            hist_imposter_bkgd_skymap[imposter][ebin].Scale(mean_bias)

        if energy_bin[ebin]>800.: continue
        for binx in range(0,hist_imposter_data_skymap[imposter][ebin].GetNbinsX()):
            for biny in range(0,hist_imposter_data_skymap[imposter][ebin].GetNbinsY()):
                mean_bias = 0.
                for imposter in range(0,n_imposters):
                    truth = hist_imposter_data_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    predict = hist_imposter_bkgd_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    bias = truth-predict
                    mean_bias += bias
                mean_bias = mean_bias/float(n_imposters)
                original = hist_real_bkgd_skymap[ebin].GetBinContent(binx+1,biny+1)
                hist_real_bkgd_skymap[ebin].SetBinContent(binx+1,biny+1,original+mean_bias)
                for imposter in range(0,n_imposters):
                    original =  hist_imposter_bkgd_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    hist_imposter_bkgd_skymap[imposter][ebin].SetBinContent(binx+1,biny+1,original+mean_bias)


for ebin in range(0,len(energy_bin)-1):
    hist_real_norm_skymap_sum.Add(hist_real_bkgd_skymap[ebin])
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_real_diff_skymap[ebin].Add(hist_real_data_skymap[ebin])
    hist_real_diff_skymap[ebin].Add(hist_real_bkgd_skymap[ebin],-1.)
    hist_real_data_skymap_sum.Add(hist_real_data_skymap[ebin])
    hist_real_bkgd_skymap_sum.Add(hist_real_bkgd_skymap[ebin])
    hist_real_diff_skymap_sum.Add(hist_real_diff_skymap[ebin])
    hist_real_data_msclw_sum.Add(hist_real_data_msclw[ebin])
    hist_real_bkgd_msclw_sum.Add(hist_real_bkgd_msclw[ebin])
    print('hist_real_data_msclw[ebin].Integral() = %s'%(hist_real_data_msclw[ebin].Integral()))
    print('hist_real_bkgd_msclw[ebin].Integral() = %s'%(hist_real_bkgd_msclw[ebin].Integral()))
    hist_real_diff_msclw_sum.Add(hist_real_data_msclw[ebin])
    hist_real_diff_msclw_sum.Add(hist_real_bkgd_msclw[ebin],-1.)
    if ebin<energy_bin_break:
        hist_real_diff_skymap_le.Add(hist_real_diff_skymap[ebin])
        hist_real_data_skymap_le.Add(hist_real_data_skymap[ebin])
        hist_real_bkgd_skymap_le.Add(hist_real_bkgd_skymap[ebin])
    else:
        hist_real_diff_skymap_he.Add(hist_real_diff_skymap[ebin])
        hist_real_data_skymap_he.Add(hist_real_data_skymap[ebin])
        hist_real_bkgd_skymap_he.Add(hist_real_bkgd_skymap[ebin])
for imposter in range(0,n_imposters):
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_imposter_diff_skymap[imposter][ebin].Add(hist_imposter_data_skymap[imposter][ebin])
        hist_imposter_diff_skymap[imposter][ebin].Add(hist_imposter_bkgd_skymap[imposter][ebin],-1.)
        hist_imposter_data_skymap_sum[imposter].Add(hist_imposter_data_skymap[imposter][ebin])
        hist_imposter_bkgd_skymap_sum[imposter].Add(hist_imposter_bkgd_skymap[imposter][ebin])
        hist_imposter_diff_skymap_sum[imposter].Add(hist_imposter_diff_skymap[imposter][ebin])

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_real_diff_skymap_smooth = CommonPlotFunctions.Smooth2DMap(hist_real_diff_skymap[ebin],smooth_size_spectroscopy,False)
    hist_real_diff_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_smooth)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_diff_skymap_reflect,None,[hist_real_diff_skymap_reflect],fig,'RA','Dec','Excess count','SkymapExcess_E%s_%s'%(ebin,plot_tag))

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_real_data_skymap_smooth = CommonPlotFunctions.Smooth2DMap(hist_real_data_skymap[ebin],smooth_size_spectroscopy,False)
    hist_real_bkgd_skymap_smooth = CommonPlotFunctions.Smooth2DMap(hist_real_bkgd_skymap[ebin],smooth_size_spectroscopy,False)
    significance_skymap = CommonPlotFunctions.GetSignificanceMap(hist_real_data_skymap_smooth,hist_real_bkgd_skymap_smooth)
    hist_real_significance_skymap_reflect = CommonPlotFunctions.reflectXaxis(significance_skymap)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_significance_skymap_reflect,None,[],fig,'RA','Dec','Significance','SkymapSignificance_E%s_%s'%(ebin,plot_tag))

hist_real_data_skymap_le_smooth = CommonPlotFunctions.Smooth2DMap(hist_real_data_skymap_le,smooth_size_spectroscopy,False)
hist_real_bkgd_skymap_le_smooth = CommonPlotFunctions.Smooth2DMap(hist_real_bkgd_skymap_le,smooth_size_spectroscopy,False)
significance_skymap = CommonPlotFunctions.GetSignificanceMap(hist_real_data_skymap_le_smooth,hist_real_bkgd_skymap_le_smooth)
hist_real_significance_skymap_le.Add(significance_skymap)

hist_real_data_skymap_he_smooth = CommonPlotFunctions.Smooth2DMap(hist_real_data_skymap_he,smooth_size_spectroscopy,False)
hist_real_bkgd_skymap_he_smooth = CommonPlotFunctions.Smooth2DMap(hist_real_bkgd_skymap_he,smooth_size_spectroscopy,False)
significance_skymap = CommonPlotFunctions.GetSignificanceMap(hist_real_data_skymap_he_smooth,hist_real_bkgd_skymap_he_smooth)
hist_real_significance_skymap_he.Add(significance_skymap)

hist_real_data_skymap_sum_smooth = CommonPlotFunctions.Smooth2DMap(hist_real_data_skymap_sum,smooth_size_spectroscopy,False)
hist_real_bkgd_skymap_sum_smooth = CommonPlotFunctions.Smooth2DMap(hist_real_bkgd_skymap_sum,smooth_size_spectroscopy,False)
significance_skymap = CommonPlotFunctions.GetSignificanceMap(hist_real_data_skymap_sum_smooth,hist_real_bkgd_skymap_sum_smooth)
hist_real_significance_skymap_sum.Add(significance_skymap)
for imposter in range(0,n_imposters):
    hist_imposter_data_skymap_sum_smooth = CommonPlotFunctions.Smooth2DMap(hist_imposter_data_skymap_sum[imposter],smooth_size_spectroscopy,False)
    hist_imposter_bkgd_skymap_sum_smooth = CommonPlotFunctions.Smooth2DMap(hist_imposter_bkgd_skymap_sum[imposter],smooth_size_spectroscopy,False)
    significance_skymap = CommonPlotFunctions.GetSignificanceMap(hist_imposter_data_skymap_sum_smooth,hist_imposter_bkgd_skymap_sum_smooth)
    hist_imposter_significance_skymap_sum[imposter].Add(significance_skymap)

hist_real_data_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_data_skymap_sum)
hist_real_bkgd_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_bkgd_skymap_sum)
CommonPlotFunctions.BackgroundSubtractMap(fig,hist_real_data_skymap_sum_reflect,hist_real_bkgd_skymap_sum_reflect,'RA','Dec','Count','SkymapBkgSubtraction_%s'%(plot_tag))

hist_real_significance_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_significance_skymap_le)
CommonPlotFunctions.MatplotlibMap2D(hist_real_significance_skymap_reflect,None,[],fig,'RA','Dec','Significance','SkymapSignificance_LE_%s'%(plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r)
hist_real_significance_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_significance_skymap_he)
CommonPlotFunctions.MatplotlibMap2D(hist_real_significance_skymap_reflect,None,[],fig,'RA','Dec','Significance','SkymapSignificance_HE_%s'%(plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r)
hist_real_significance_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_significance_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_real_significance_skymap_reflect,None,[],fig,'RA','Dec','Significance','SkymapSignificance_Sum_%s'%(plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r)
for imposter in range(0,n_imposters):
    hist_imposter_significance_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_imposter_significance_skymap_sum[imposter])
    CommonPlotFunctions.MatplotlibMap2D(hist_imposter_significance_skymap_reflect,None,[],fig,'RA','Dec','Significance','SkymapSignificanceImposter%s_Sum_%s'%(imposter,plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r)

hist_elev_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_elev_skymap)
CommonPlotFunctions.MatplotlibMap2D(hist_elev_skymap_reflect,hist_elev_skymap_reflect,[],fig,'RA','Dec','Elevation [deg]','SkymapElev')
hist_azim_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_azim_skymap)
CommonPlotFunctions.MatplotlibMap2D(hist_azim_skymap_reflect,hist_azim_skymap_reflect,[],fig,'RA','Dec','Azimuth [deg]','SkymapAzim')
hist_nsb_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_nsb_skymap)
CommonPlotFunctions.MatplotlibMap2D(hist_nsb_skymap_reflect,hist_nsb_skymap_reflect,[],fig,'RA','Dec','NSB','SkymapNSB')

MakeFluxMap(hist_real_flux_skymap, hist_real_data_skymap, hist_real_bkgd_skymap, hist_real_norm_skymap_sum, hist_elev_skymap)
for imp in range(0,n_imposters):
    MakeFluxMap(hist_imposter_flux_skymap[imp], hist_imposter_data_skymap[imp], hist_imposter_bkgd_skymap[imp], hist_real_norm_skymap_sum, hist_elev_skymap)

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_real_flux_skymap_sum.Add(hist_real_flux_skymap[ebin])
    if ebin<energy_bin_break:
        hist_real_flux_skymap_le.Add(hist_real_flux_skymap[ebin])
    else:
        hist_real_flux_skymap_he.Add(hist_real_flux_skymap[ebin])
    for imposter in range(0,n_imposters):
        hist_imposter_flux_skymap_sum[imposter].Add(hist_imposter_flux_skymap[imposter][ebin])
        if ebin<energy_bin_break:
            hist_imposter_flux_skymap_le[imposter].Add(hist_imposter_flux_skymap[imposter][ebin])
        else:
            hist_imposter_flux_skymap_he[imposter].Add(hist_imposter_flux_skymap[imposter][ebin])

for lbin in range(0,hist_real_diff_msclw_sum.GetNbinsX()):
    for wbin in range(0,hist_real_diff_msclw_sum.GetNbinsY()):
        bkg_error = pow(abs(hist_real_bkgd_msclw_sum.GetBinContent(lbin+1,wbin+1)),0.5)
        significance = hist_real_diff_msclw_sum.GetBinContent(lbin+1,wbin+1)/max(1.,bkg_error)
        hist_real_sign_msclw_sum.SetBinContent(lbin+1,wbin+1,significance)
CommonPlotFunctions.MatplotlibHist2D(hist_real_sign_msclw_sum,fig,'scaled length','scaled width','Significance','MatrixSign_Sum_%s'%(plot_tag),zmax=5,zmin=-5)
CommonPlotFunctions.MatplotlibHist2D(hist_real_diff_msclw_sum,fig,'scaled length','scaled width','Residual','MatrixDiff_Sum_%s'%(plot_tag))

MakeExtensionProfile(region_x[0],region_y[0],region_r[0],do_fit,region_name,hist_real_flux_skymap_le,hist_imposter_flux_skymap_le,'le')
MakeExtensionProfile(region_x[0],region_y[0],region_r[0],do_fit,region_name,hist_real_flux_skymap_he,hist_imposter_flux_skymap_he,'he')
MakeExtensionProfile(region_x[0],region_y[0],region_r[0],do_fit,region_name,hist_real_flux_skymap_sum,hist_imposter_flux_skymap_sum,'sum')
MakeSpectrum(region_x,region_y,region_r,region_name,excl_region_x,excl_region_y,excl_region_r)

hist_real_flux_skymap_sum = CommonPlotFunctions.Smooth2DMap(hist_real_flux_skymap_sum,smooth_size_spectroscopy,False)
hist_real_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_reflect,None,[],fig,'RA','Dec','$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFlux_Sum_%s'%(plot_tag),psf=0.08)
SaveFITS(hist_real_flux_skymap_sum)

hist_real_diff_skymap_sum = CommonPlotFunctions.Smooth2DMap(hist_real_diff_skymap_sum,smooth_size_spectroscopy,False)
hist_real_diff_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_real_diff_skymap_reflect,None,[],fig,'RA','Dec','Excess count','SkymapExcess_Sum_%s'%(plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r,psf=0.08)
hist_real_diff_skymap_le = CommonPlotFunctions.Smooth2DMap(hist_real_diff_skymap_le,smooth_size_spectroscopy,False)
hist_real_diff_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_le)
CommonPlotFunctions.MatplotlibMap2D(hist_real_diff_skymap_reflect,None,[],fig,'RA','Dec','Excess count','SkymapExcess_LE_%s'%(plot_tag),psf=0.08)
hist_real_diff_skymap_he = CommonPlotFunctions.Smooth2DMap(hist_real_diff_skymap_he,smooth_size_spectroscopy,False)
hist_real_diff_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_he)
CommonPlotFunctions.MatplotlibMap2D(hist_real_diff_skymap_reflect,None,[],fig,'RA','Dec','Excess count','SkymapExcess_HE_%s'%(plot_tag),psf=0.08)
SaveFITS(hist_real_diff_skymap_le)
SaveFITS(hist_real_diff_skymap_he)

fig.clf()
fig.set_figheight(figsize_x)
fig.set_figwidth(figsize_x)
axbig = fig.add_subplot()
cycol = cycle('rgbcmy')
z_max_range = 6.
z_min_range = -3.
zscore_normal = GetMapNormalDistribution(hist_real_bkgd_skymap_sum)
axbig.hist(zscore_normal, bins=40, range=[z_min_range, z_max_range], density=True, color='gray', alpha=0.5, label='Normal')
for imposter in range(0,n_imposters):
    next_color = next(cycol)
    zscores = GetMapChi2Distribution(hist_imposter_significance_skymap_sum[imposter],hist_imposter_bkgd_skymap_sum[imposter])
    axbig.hist(zscores, bins=40, range=[z_min_range, z_max_range], density=True, histtype='step', color=next_color, alpha=1.0, label='Mimic #%s'%(imposter+1))
zscores = GetMapChi2Distribution(hist_real_significance_skymap_sum,hist_real_bkgd_skymap_sum)
axbig.hist(zscores, bins=40, range=[z_min_range, z_max_range], density=True, histtype='step', color='black', alpha=1.0, label='ON data')
axbig.set_xlabel('Significances')
axbig.set_ylabel('Entries')
axbig.set_yscale('log')
axbig.legend(loc='best')
fig.savefig("output_plots/SignificanceDistribution_%s.png"%(plot_tag),bbox_inches='tight')
axbig.remove()

if 'PSR_J1907_p0602' in source_name:

    Hist_Fermi = ROOT.TH2D("Hist_Fermi","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    MWL_map_file = '/home/rshang/MatrixDecompositionMethod/MWL_maps/rg_removed_residualmap_30gev-2tev_pointsource_powerlaw_2.00_residmap.fits'
    Hist_Fermi = CommonPlotFunctions.GetFITSMap(MWL_map_file, Hist_Fermi, True)
    Hist_Fermi = CommonPlotFunctions.Smooth2DMap(Hist_Fermi,0.1,False)
    #Hist_Fermi = CommonPlotFunctions.ConvertTSmapToZscore(Hist_Fermi)
    Hist_Fermi_reflect = CommonPlotFunctions.reflectXaxis(Hist_Fermi)
    CommonPlotFunctions.MatplotlibMap2D(Hist_Fermi_reflect,None,[],fig,'RA','Dec','Excess count','SkymapFermi_30GeV_%s'%(plot_tag),psf=0.1)

    hist_real_diff_skymap_le_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_le)
    hist_real_diff_skymap_he_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_he)

    hist_real_flux_skymap_sum = CommonPlotFunctions.Smooth2DMap(hist_real_flux_skymap_sum,smooth_size_spectroscopy,False)
    hist_real_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_sum)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_reflect,None,[hist_real_diff_skymap_he_reflect,hist_real_diff_skymap_le_reflect,Hist_Fermi_reflect],fig,'RA','Dec','$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFlux_Sum_%s'%(plot_tag),colormap='gray',psf=0.08)


    Hist_mc_intensity = ROOT.TH2D("Hist_mc_intensity","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    Hist_mc_column = ROOT.TH2D("Hist_mc_column","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    pc_to_cm = 3.086e+18
    CO_intensity_to_H_column_density = 2.*1e20
    # Dame, T. M.; Hartmann, Dap; Thaddeus, P., 2011, "Replication data for: First Quadrant, main survey (DHT08)", https://doi.org/10.7910/DVN/1PG9NV, Harvard Dataverse, V3
    # https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PG9NV
    FITS_correction = 1000.# the source FITS file has a mistake in velocity km/s -> m/s
    MWL_map_file = '/home/rshang/MatrixDecompositionMethod/MWL_maps/DHT08_Quad1_interp.fits' 
    CommonPlotFunctions.GetSlicedDataCubeMap(MWL_map_file, Hist_mc_intensity, 10., 40.)
    Hist_mc_column.Reset()
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,[hist_real_diff_skymap_he_reflect,hist_real_diff_skymap_le_reflect,Hist_Fermi_reflect],fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapRadioCOMap_p10p40_%s'%(plot_tag),roi_x=[287.1,287.1],roi_y=[6.5,6.5],roi_r=[0.4,0.8],colormap='gray')
    MWL_map_file = '/home/rshang/MatrixDecompositionMethod/MWL_maps/DHT08_Quad1_interp.fits' 
    CommonPlotFunctions.GetSlicedDataCubeMap(MWL_map_file, Hist_mc_intensity, 40., 70.)
    Hist_mc_column.Reset()
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,[hist_real_diff_skymap_he_reflect,hist_real_diff_skymap_le_reflect,Hist_Fermi_reflect],fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapRadioCOMap_p40p70_%s'%(plot_tag),roi_x=[287.1,287.1],roi_y=[6.5,6.5],roi_r=[0.4,0.8],colormap='gray')

    MWL_map_file = '/home/rshang/MatrixDecompositionMethod/MWL_maps/DHT08_Quad1_interp.fits' 
    vel_axis_inner, column_density_axis_inner = CommonPlotFunctions.GetVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.0, 0.4)
    vel_axis_outer, column_density_axis_outer = CommonPlotFunctions.GetVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.4, 0.8)
    column_density_axis_inner = CO_intensity_to_H_column_density*np.array(column_density_axis_inner)
    column_density_axis_outer = CO_intensity_to_H_column_density*np.array(column_density_axis_outer)
    column_density_axis_diff = column_density_axis_outer - column_density_axis_inner

    max_idx = np.argmax(column_density_axis_diff)
    print ('CO velocity of highest emission = %0.1f km/s'%(vel_axis_inner[max_idx]))
    min_idx = np.argmin(column_density_axis_diff)
    print ('CO velocity of lowest emission = %0.1f km/s'%(vel_axis_inner[min_idx]))

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.plot(vel_axis_inner, column_density_axis_inner)
    axbig.plot(vel_axis_inner, column_density_axis_outer)
    axbig.set_xlabel('$V_{LSR}$ [km/s]')
    axbig.set_ylabel('column density per channel [$1/cm^{2}/(km/s)$]')
    fig.savefig("output_plots/VelocitySpectrumCO_Cavity.png",bbox_inches='tight')
    axbig.remove()
    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.plot(vel_axis_inner, column_density_axis_diff)
    axbig.set_xlabel('$V_{LSR}$ [km/s]')
    axbig.set_ylabel('column density per channel [$1/cm^{2}/(km/s)$]')
    fig.savefig("output_plots/VelocitySpectrumCO_Diff.png",bbox_inches='tight')
    axbig.remove()

    MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/GALFA_HI_RA+DEC_284.00+02.35_N.fits' 
    CommonPlotFunctions.GetSlicedGalfaHIDataCubeMap(MWL_map_file, Hist_mc_intensity, 10.*1e3, 40.*1e3, True)
    MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/GALFA_HI_RA+DEC_284.00+10.35_N.fits' 
    CommonPlotFunctions.GetSlicedGalfaHIDataCubeMap(MWL_map_file, Hist_mc_intensity, 10.*1e3, 40.*1e3, False)
    MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/GALFA_HI_RA+DEC_292.00+02.35_N.fits' 
    CommonPlotFunctions.GetSlicedGalfaHIDataCubeMap(MWL_map_file, Hist_mc_intensity, 10.*1e3, 40.*1e3, False)
    MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/GALFA_HI_RA+DEC_292.00+10.35_N.fits' 
    CommonPlotFunctions.GetSlicedGalfaHIDataCubeMap(MWL_map_file, Hist_mc_intensity, 10.*1e3, 40.*1e3, False)
    Hist_mc_column.Reset()
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,[hist_real_diff_skymap_he_reflect,hist_real_diff_skymap_le_reflect,Hist_Fermi_reflect],fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapRadioHIMap_p10p40_%s'%(plot_tag),roi_x=[287.1,287.1],roi_y=[6.5,6.5],roi_r=[0.4,0.8],colormap='gray')

    #MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/GALFA_HI_RA+DEC_284.00+02.35_N.fits' 
    #vel_axis_inner, column_density_axis_inner = CommonPlotFunctions.GetGalfaHIVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.0, 0.4)
    #vel_axis_outer, column_density_axis_outer = CommonPlotFunctions.GetGalfaHIVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.4, 0.8)
    #column_density_axis_inner = CO_intensity_to_H_column_density*np.array(column_density_axis_inner)
    #column_density_axis_outer = CO_intensity_to_H_column_density*np.array(column_density_axis_outer)
    #MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/GALFA_HI_RA+DEC_284.00+10.35_N.fits' 
    #vel_axis_inner, column_density_axis_inner = CommonPlotFunctions.GetGalfaHIVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.0, 0.4)
    #vel_axis_outer, column_density_axis_outer = CommonPlotFunctions.GetGalfaHIVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.4, 0.8)
    #column_density_axis_inner += CO_intensity_to_H_column_density*np.array(column_density_axis_inner)
    #column_density_axis_outer += CO_intensity_to_H_column_density*np.array(column_density_axis_outer)
    #MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/GALFA_HI_RA+DEC_292.00+02.35_N.fits' 
    #vel_axis_inner, column_density_axis_inner = CommonPlotFunctions.GetGalfaHIVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.0, 0.4)
    #vel_axis_outer, column_density_axis_outer = CommonPlotFunctions.GetGalfaHIVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.4, 0.8)
    #column_density_axis_inner += CO_intensity_to_H_column_density*np.array(column_density_axis_inner)
    #column_density_axis_outer += CO_intensity_to_H_column_density*np.array(column_density_axis_outer)
    #MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/GALFA_HI_RA+DEC_292.00+10.35_N.fits' 
    #vel_axis_inner, column_density_axis_inner = CommonPlotFunctions.GetGalfaHIVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.0, 0.4)
    #vel_axis_outer, column_density_axis_outer = CommonPlotFunctions.GetGalfaHIVelocitySpectrum(MWL_map_file, 40.7, -0.8, 0.4, 0.8)
    #column_density_axis_inner += CO_intensity_to_H_column_density*np.array(column_density_axis_inner)
    #column_density_axis_outer += CO_intensity_to_H_column_density*np.array(column_density_axis_outer)

    #column_density_axis_diff = column_density_axis_outer - column_density_axis_inner

    #max_idx = np.argmax(column_density_axis_diff)
    #print ('HI velocity of highest emission = %0.1f km/s'%(vel_axis_inner[max_idx]))
    #min_idx = np.argmin(column_density_axis_diff)
    #print ('HI velocity of lowest emission = %0.1f km/s'%(vel_axis_inner[min_idx]))

    #fig.clf()
    #fig.set_figheight(figsize_y)
    #fig.set_figwidth(figsize_x)
    #axbig = fig.add_subplot()
    #axbig.plot(vel_axis_inner, column_density_axis_inner)
    #axbig.plot(vel_axis_inner, column_density_axis_outer)
    #axbig.set_xlabel('$V_{LSR}$ [km/s]')
    #axbig.set_ylabel('column density per channel [$1/cm^{2}/(km/s)$]')
    #fig.savefig("output_plots/VelocitySpectrumHI_Cavity.png",bbox_inches='tight')
    #axbig.remove()
    #fig.clf()
    #fig.set_figheight(figsize_y)
    #fig.set_figwidth(figsize_x)
    #axbig = fig.add_subplot()
    #axbig.plot(vel_axis_inner, column_density_axis_diff)
    #axbig.set_xlabel('$V_{LSR}$ [km/s]')
    #axbig.set_ylabel('column density per channel [$1/cm^{2}/(km/s)$]')
    #fig.savefig("output_plots/VelocitySpectrumHI_Diff.png",bbox_inches='tight')
    #axbig.remove()


    #Hist_Hawc = ROOT.TH2D("Hist_Hawc","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    #Hist_Hawc.Rebin2D(3,3)
    #hawc_map_list = []
    #hawc_map_list += ['cd'] # 1-3.16 TeV
    #hawc_map_list += ['ef'] # 3.16-10 TeV
    #hawc_map_list += ['gh'] # 10-31.6 TeV
    #hawc_map_list += ['ij'] # 31.6-100 TeV
    #hawc_map_list += ['kl'] # 100-316 TeV
    #hawc_map_list += ['gl'] # 10-316 TeV
    #for hfile in range(0,len(hawc_map_list)):
    #    #MWL_map_file = '/home/rshang/MatrixDecompositionMethod/MWL_maps/%s-gaussGDE.fits'%(hawc_map_list[hfile])
    #    MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/%s-pointlike.fits.gz'%(hawc_map_list[hfile])
    #    Hist_Hawc = CommonPlotFunctions.GetHealpixMap(MWL_map_file, Hist_Hawc, True)
    #    Hist_Hawc_reflect = CommonPlotFunctions.reflectXaxis(Hist_Hawc)
    #    CommonPlotFunctions.MatplotlibMap2D(Hist_Hawc_reflect,None,[],fig,'RA','Dec','Significance','SkymapHAWC_%s_%s'%(hawc_map_list[hfile],plot_tag),colormap='magma',psf=0.13)

    #Hist_Tobias = ROOT.TH2D("Hist_Tobias","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    #MWL_map_file = '/home/rshang/MatrixDecompositionMethod/MWL_maps/TobiasNewMap.fits'
    #Hist_Tobias = CommonPlotFunctions.GetFITSMap(MWL_map_file, Hist_Tobias, True)
    #Hist_Tobias_reflect = CommonPlotFunctions.reflectXaxis(Hist_Tobias)
    #CommonPlotFunctions.MatplotlibMap2D(Hist_Tobias_reflect,None,[],fig,'RA','Dec','Significance','SkymapTobias_%s'%(plot_tag))

    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('Fit 2d Gaussian (LE)')
    fit_2d_model(hist_real_data_skymap_le, hist_real_bkgd_skymap_le, 286.98, 6.04)
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('Fit 2d Gaussian (HE)')
    fit_2d_model(hist_real_data_skymap_he, hist_real_bkgd_skymap_he, 286.98, 6.04)
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('Fit 2d Gaussian (sum)')
    fit_2d_model(hist_real_data_skymap_sum, hist_real_bkgd_skymap_sum, 286.98, 6.04)
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

elif 'PSR_J2032_p4127' in source_name:

    hist_real_diff_skymap_le_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_le)
    hist_real_diff_skymap_he_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_he)
    Hist_mc_intensity = ROOT.TH2D("Hist_mc_intensity","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    Hist_mc_column = ROOT.TH2D("Hist_mc_column","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    pc_to_cm = 3.086e+18
    CO_intensity_to_H_column_density = 2.*1e20
    # Dame, T. M.; Hartmann, Dap; Thaddeus, P., 2011, "Replication data for: First Quadrant, main survey (DHT08)", https://doi.org/10.7910/DVN/1PG9NV, Harvard Dataverse, V3
    # https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PG9NV
    FITS_correction = 1000.# the source FITS file has a mistake in velocity km/s -> m/s
    MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/DHT10_Cygnus_interp.fits' 
    CommonPlotFunctions.GetSlicedDataCubeMap(MWL_map_file, Hist_mc_intensity, -20., 20.)
    Hist_mc_column.Reset()
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,[],fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapRadioCOMap_%s'%(plot_tag),colormap='gray')

    MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/DHT10_Cygnus_interp.fits' 
    vel_axis_inner, column_density_axis_inner = CommonPlotFunctions.GetVelocitySpectrum(MWL_map_file, 80.22, 1.04, 0.0, 0.4)
    vel_axis_outer, column_density_axis_outer = CommonPlotFunctions.GetVelocitySpectrum(MWL_map_file, 80.22, 1.04, 0.4, 0.8)
    column_density_axis_inner = CO_intensity_to_H_column_density*np.array(column_density_axis_inner)
    column_density_axis_outer = CO_intensity_to_H_column_density*np.array(column_density_axis_outer)
    column_density_axis_diff = column_density_axis_outer - column_density_axis_inner

    max_idx = np.argmax(column_density_axis_diff)
    print ('CO velocity of highest emission = %0.1f km/s'%(vel_axis_inner[max_idx]))
    min_idx = np.argmin(column_density_axis_diff)
    print ('CO velocity of lowest emission = %0.1f km/s'%(vel_axis_inner[min_idx]))

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.plot(vel_axis_inner, column_density_axis_inner)
    axbig.plot(vel_axis_inner, column_density_axis_outer)
    axbig.set_xlabel('$V_{LSR}$ [km/s]')
    axbig.set_ylabel('column density per channel [$1/cm^{2}/(km/s)$]')
    fig.savefig("output_plots/VelocitySpectrumCO_Cavity.png",bbox_inches='tight')
    axbig.remove()
    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.plot(vel_axis_inner, column_density_axis_diff)
    axbig.set_xlabel('$V_{LSR}$ [km/s]')
    axbig.set_ylabel('column density per channel [$1/cm^{2}/(km/s)$]')
    fig.savefig("output_plots/VelocitySpectrumCO_Diff.png",bbox_inches='tight')
    axbig.remove()


elif 'PSR_J1856_p0245' in source_name:

    hist_real_diff_skymap_le_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_le)
    hist_real_diff_skymap_he_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_he)
    Hist_mc_intensity = ROOT.TH2D("Hist_mc_intensity","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    Hist_mc_column = ROOT.TH2D("Hist_mc_column","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    pc_to_cm = 3.086e+18
    CO_intensity_to_H_column_density = 2.*1e20
    # Dame, T. M.; Hartmann, Dap; Thaddeus, P., 2011, "Replication data for: First Quadrant, main survey (DHT08)", https://doi.org/10.7910/DVN/1PG9NV, Harvard Dataverse, V3
    # https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PG9NV
    FITS_correction = 1000.# the source FITS file has a mistake in velocity km/s -> m/s
    MWL_map_file = '/home/rshang/MatrixDecompositionMethod/MWL_maps/DHT08_Quad1_interp.fits' 
    CommonPlotFunctions.GetSlicedDataCubeMap(MWL_map_file, Hist_mc_intensity, 81., 102.)
    Hist_mc_column.Reset()
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,[hist_real_diff_skymap_he_reflect,hist_real_diff_skymap_le_reflect],fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapRadioCOMap_p81p102_%s'%(plot_tag),colormap='gray')

    MWL_map_file = '/gamma_raid/userspace/rshang/MW_FITS/GALFA_HI_RA+DEC_284.00+02.35_N.fits' 
    CommonPlotFunctions.GetSlicedGalfaHIDataCubeMap(MWL_map_file, Hist_mc_intensity, 81.*1e3, 102.*1e3, True)
    Hist_mc_column.Reset()
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,[hist_real_diff_skymap_he_reflect,hist_real_diff_skymap_le_reflect],fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapRadioHIMap_p10p40_%s'%(plot_tag),colormap='gray')

elif 'SS433' in source_name:

    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('Fit 2d Gaussian')
    fit_2d_model(hist_real_data_skymap_sum, hist_real_bkgd_skymap_sum, 288.0833333, 4.9166667)
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')


print ('total_data_expo = %0.1f hrs'%(total_data_expo))
