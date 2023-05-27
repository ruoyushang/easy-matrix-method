
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

energy_bin = CommonPlotFunctions.energy_bin
energy_bin_cut_low = int(sys.argv[4])
energy_bin_cut_up = int(sys.argv[5])
doImposter = int(sys.argv[6])
source_name = sys.argv[1]
input_epoch = sys.argv[2] # 'V5' or 'V6' or 'V5V6'
isON = sys.argv[3]  # 'ON' or 'OFF'

analysis_method = CommonPlotFunctions.analysis_method
n_xoff_bins = CommonPlotFunctions.n_xoff_bins
smooth_size_spectroscopy = CommonPlotFunctions.smooth_size_spectroscopy

n_imposters = 5
if not doImposter:
    n_imposters = 0

plot_tag = source_name
plot_tag += '_'+analysis_method
plot_tag += '_E'+sys.argv[4]+'_'+sys.argv[5]

list_epoch = []
if 'V5' in input_epoch:
    list_epoch += ['V5']
if 'V6' in input_epoch:
    list_epoch += ['V6']

total_data_expo = 0.

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

    str_flux_calibration_el80 = ['6.63e-11', '1.35e-11', '5.50e-12', '2.59e-12', '1.30e-12', '5.71e-13', '1.86e-13', '5.01e-14', '1.39e-14', '7.30e-15', '3.37e-15']
    str_flux_calibration_el70 = ['6.71e-11', '1.44e-11', '6.42e-12', '2.68e-12', '1.24e-12', '5.17e-13', '1.55e-13', '4.46e-14', '1.47e-14', '5.40e-15', '1.27e-15']
    str_flux_calibration_el60 = ['6.73e-11', '1.81e-11', '8.77e-12', '3.81e-12', '1.47e-12', '6.10e-13', '1.58e-13', '3.99e-14', '1.43e-14', '4.17e-15', '1.36e-15']
    str_flux_calibration_el50 = ['1.28e-10', '1.81e-11', '1.05e-11', '4.90e-12', '2.46e-12', '7.91e-13', '1.45e-13', '3.92e-14', '8.53e-15', '3.36e-15', '1.30e-15']

    flux_calibration_el80 = []
    flux_calibration_el70 = []
    flux_calibration_el60 = []
    flux_calibration_el50 = []
    for string in str_flux_calibration_el80:
        flux_calibration_el80.append(float(string))
    for string in str_flux_calibration_el70:
        flux_calibration_el70.append(float(string))
    for string in str_flux_calibration_el60:
        flux_calibration_el60.append(float(string))
    for string in str_flux_calibration_el50:
        flux_calibration_el50.append(float(string))

    xp = [55.,65.,75.,85.]
    fp = [flux_calibration_el50[energy],flux_calibration_el60[energy],flux_calibration_el70[energy],flux_calibration_el80[energy]]

    return np.interp(elev, xp, fp)

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
        energy_edge_lo += [pow(10,energy_mean_log[eb]-energy_log_delta)]
        energy_edge_hi += [pow(10,energy_mean_log[eb]+energy_log_delta)]
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

def GetHawcSaraFluxJ1908():

    energies = [1.38,2.58,4.58,6.89,10.88,18.23,35.56,61.32,107.69,186.29]
    fluxes = [2.7211e-11,2.7457e-11,2.2287e-11,1.7933e-11,1.3624e-11,9.4346e-12,5.9615e-12,2.5984e-12,1.4372e-12,5.7010e-13]
    flux_errs = [2.7211e-11,2.7457e-11,2.2287e-11,1.7933e-11,1.3624e-11,9.4346e-12,5.9615e-12,2.5984e-12,1.4372e-12,5.7010e-13]
    flux_errs_up = [+1.3091e-12,+1.1941e-12,+1.0624e-12,+7.6946e-13,+5.3155e-13,+5.0018e-13,+3.1416e-13,+2.5557e-13,+2.3385e-13,+1.7388e-13]
    flux_errs_low = [-1.2902e-12,-1.2274e-12,-1.0946e-12,-7.6676e-13,-5.2608e-13,-5.1950e-13,-3.2465e-13,-2.4536e-13,-2.3011e-13,-1.7487e-13]

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

def GetVeritasTobiasFluxJ1908():

    energy_edges = [794,1580,3160,6310,12600]
    energies = []
    for edge in range(0,len(energy_edges)-1):
        energies += [0.5*(energy_edges[edge]+energy_edges[edge+1])]
    fluxes = [8.96e-12, 8.99e-12, 7.49e-12, 4.41e-12]
    flux_errs = []
    fluxes_imp = []
    fluxes_imp += [[4.79e-13, 1.23e-12, 2.00e-13, -9.25e-13]]
    fluxes_imp += [[-5.14e-13, 8.54e-14, 4.01e-13, -3.21e-13]]
    fluxes_imp += [[5.16e-13, 7.70e-13, 1.16e-12, -5.53e-13]]
    fluxes_imp += [[-2.54e-13, 2.10e-12, 6.83e-13, -4.54e-13]]
    fluxes_imp += [[-5.66e-13, -4.65e-13, -5.27e-13, -3.10e-13]]

    for entry in range(0,len(energies)):
        syst_err = 0.
        for imp in range(0,len(fluxes_imp)):
            syst_err += pow(fluxes_imp[imp][entry],2)
        syst_err = pow(syst_err/float(len(fluxes_imp)-1),0.5)
        flux_errs += [syst_err]

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
            if doImposter==1 and imposter_bkgd_list[imposter][ebin]>0. and imposter_data_avg[ebin]>0.:
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
    axbig = fig.add_subplot()
    if 'Crab' in source_name:
        axbig.plot(xdata, ydata_crab,'r-',label='1508.06442', zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_stat_err, bottom=real_flux-real_flux_stat_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',zorder=3)
    elif source_name=='PSR_J1907_p0602':
        log_energy = np.linspace(log10(1e2),log10(1e5),50)
        xdata_ref = pow(10.,log_energy)
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHawcDiffusionFluxJ1908()
        Sara_energies, Sara_fluxes, Sara_flux_errs = GetHawcSaraFluxJ1908()
        HESS_energies, HESS_fluxes, HESS_flux_errs = GetHessFluxJ1908()
        OldV_energies, OldV_fluxes, OldV_flux_errs = GetVeritasFluxJ1908()
        Jordan_energies, Jordan_fluxes, Jordan_flux_errs = GetFermiJordanFluxJ1908()
        Tobias_energies, Tobias_fluxes, Tobias_flux_errs = GetVeritasTobiasFluxJ1908()
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ1908()

        axbig.errorbar(Jordan_energies,Jordan_fluxes,Jordan_flux_errs,color='b',marker='s',ls='none',label='Fermi-LAT',zorder=4)
        axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='r',marker='s',ls='none',label='HAWC (2021 paper)',zorder=3)
        #axbig.errorbar(Sara_energies,Sara_fluxes,Sara_flux_errs,color='c',marker='s',ls='none',label='HAWC (Sara)',zorder=3)
        axbig.errorbar(HESS_energies,HESS_fluxes,HESS_flux_errs,color='g',marker='s',ls='none',label='HESS',zorder=2)
        #axbig.errorbar(OldV_energies,OldV_fluxes,OldV_flux_errs,color='orange',marker='s',ls='none',label='VERITAS (2014)',zorder=1)
        axbig.errorbar(Tobias_energies,Tobias_fluxes,Tobias_flux_errs,color='orange',marker='s',ls='none',label='VERITAS (Tobias)',zorder=1)
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='m',marker='s',ls='none',label='LHAASO',zorder=7)

        axbig.bar(energy_axis, 2.*real_flux_stat_err, bottom=real_flux-real_flux_stat_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS')

        PrintSpectralDataForNaima(Tobias_energies,Tobias_fluxes,Tobias_flux_errs,'Tobias')
        PrintSpectralDataForNaima(HESS_energies,HESS_fluxes,HESS_flux_errs,'HESS')
        PrintSpectralDataForNaima(Jordan_energies,Jordan_fluxes,Jordan_flux_errs,'Fermi')
        PrintSpectralDataForNaima(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,'LHAASO')
        PrintSpectralDataForNaima(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,'HAWC')

    else:
        axbig.bar(energy_axis, 2.*real_flux_stat_err, bottom=real_flux-real_flux_stat_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_total_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS')
        PrintSpectralDataForNaima(energy_axis,real_flux,real_flux_total_err,'VERITAS')

    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    axbig.set_xscale('log')
    axbig.set_yscale('log')
    axbig.legend(loc='best')
    plotname = 'RealSpectrum_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    zscore = []
    for eb in range(0,len(energy_axis)):
        zscore += [real_flux[eb]/pow(pow(real_flux_stat_err[eb],2)+pow(real_flux_syst_err[eb],2),0.5)]
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

                #correction = GetFluxCalibration(ebin,elev_content)

                norm_ratio = norm_content/norm_content_max
                norm_weight = 1./(1.+np.exp(-(norm_ratio-0.3)/0.1))

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
                stat_data_err = pow(max(data_content,0.),0.5)
                flux_stat_err = max(stat_data_err,1.)/local_areatime*pow(energy_bin[ebin]/1e3,2)/(100.*100.*3600.)/delta_E
                flux_map[ebin].SetBinContent(binx+1,biny+1,flux_content*norm_weight)
                flux_map[ebin].SetBinError(binx+1,biny+1,flux_stat_err*norm_weight)


def MakeExtensionProfile(roi_x,roi_y,roi_r,fit_profile,roi_name,real_map,imposter_maps,erange_tag):

    if real_map.Integral()==0.:
        return

    real_profile, real_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension(real_map,roi_x,roi_y,3.0*roi_r)
    imposter_profile_list = []
    imposter_profile_err_list = []
    for imposter in range(0,n_imposters):
        imposter_profile, imposter_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension(imposter_maps[imposter],roi_x,roi_y,3.0*roi_r)
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
    #axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    #axbig.errorbar(theta2,real_profile,real_profile_stat_err,color='k',marker='s',ls='none',label='ON data')
    axbig.errorbar(theta2,real_profile,real_profile_total_err,color='k',marker='s',ls='none',label='ON data')
    axbig.bar(theta2, 2.*real_profile_stat_err, bottom=-real_profile_stat_err+real_profile, width=1.*theta2_err, color='b', align='center', alpha=0.2)
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

InputFile.Close()

hist_areatime_skymap = ROOT.TH2D("hist_areatime_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_elev_skymap = ROOT.TH2D("hist_elev_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_azim_skymap = ROOT.TH2D("hist_azim_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_nsb_skymap = ROOT.TH2D("hist_nsb_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)

hist_real_flux_skymap_sum = ROOT.TH2D("hist_real_flux_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_data_skymap_sum = ROOT.TH2D("hist_real_data_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_bkgd_skymap_sum = ROOT.TH2D("hist_real_bkgd_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_norm_skymap_sum = ROOT.TH2D("hist_real_norm_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_diff_skymap_sum = ROOT.TH2D("hist_real_diff_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_significance_skymap_sum = ROOT.TH2D("hist_real_significance_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_flux_skymap = []
hist_real_data_skymap = []
hist_real_bkgd_skymap = []
hist_real_diff_skymap = []
for ebin in range(0,len(energy_bin)-1):
    hist_real_flux_skymap += [ROOT.TH2D("hist_real_flux_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_data_skymap += [ROOT.TH2D("hist_real_data_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_bkgd_skymap += [ROOT.TH2D("hist_real_bkgd_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_diff_skymap += [ROOT.TH2D("hist_real_diff_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]

hist_imposter_data_skymap_sum = []
hist_imposter_bkgd_skymap_sum = []
hist_imposter_diff_skymap_sum = []
hist_imposter_flux_skymap_sum = []
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
    for yoff_idx in range(0,n_xoff_bins):
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
                        if effective_area[ebin] < 30000.: continue
                    HistName = "Hist_OnData_SR_Skymap_Sum_ErecS%sto%s"%(int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                    FillSkyMapHistogram(InputFile.Get(HistName),hist_real_data_skymap[ebin])
                    HistName = "Hist_OnData_CR_Skymap_%s_Sum_ErecS%sto%s"%(analysis_method,int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                    FillSkyMapHistogram(InputFile.Get(HistName),hist_real_bkgd_skymap[ebin])
                    hist_real_diff_skymap[ebin].Add(hist_real_data_skymap[ebin])
                    hist_real_diff_skymap[ebin].Add(hist_real_bkgd_skymap[ebin],-1.)
                InputFile.Close()

hist_elev_skymap.Scale(1./n_samples)
hist_azim_skymap.Scale(1./n_samples)
hist_nsb_skymap.Scale(1./n_samples)

if doImposter:
    for imposter in range(0,n_imposters):
        n_imposter_samples = 0.
        for xoff_idx in range(0,n_xoff_bins):
            for yoff_idx in range(0,n_xoff_bins):
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
                                if effective_area[ebin] < 30000.: continue
                            HistName = "Hist_OnData_SR_Skymap_Sum_ErecS%sto%s"%(int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                            FillSkyMapHistogram(InputFile.Get(HistName),hist_imposter_data_skymap[imposter][ebin])
                            HistName = "Hist_OnData_CR_Skymap_%s_Sum_ErecS%sto%s"%(analysis_method,int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                            FillSkyMapHistogram(InputFile.Get(HistName),hist_imposter_bkgd_skymap[imposter][ebin])
                            hist_imposter_diff_skymap[imposter][ebin].Add(hist_imposter_data_skymap[imposter][ebin])
                            hist_imposter_diff_skymap[imposter][ebin].Add(hist_imposter_bkgd_skymap[imposter][ebin],-1.)
                        InputFile.Close()

for ebin in range(0,len(energy_bin)-1):
    hist_real_norm_skymap_sum.Add(hist_real_bkgd_skymap[ebin])
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_real_data_skymap_sum.Add(hist_real_data_skymap[ebin])
    hist_real_bkgd_skymap_sum.Add(hist_real_bkgd_skymap[ebin])
    hist_real_diff_skymap_sum.Add(hist_real_diff_skymap[ebin])
for imposter in range(0,n_imposters):
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
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

hist_real_significance_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_significance_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_real_significance_skymap_reflect,None,[],fig,'RA','Dec','Significance','SkymapSignificance_Sum_%s'%(plot_tag))

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
    for imposter in range(0,n_imposters):
        hist_imposter_flux_skymap_sum[imposter].Add(hist_imposter_flux_skymap[imposter][ebin])

excl_region_x = MapCenter_x
excl_region_y = MapCenter_y
excl_region_r = 0.0
region_x = MapCenter_x
region_y = MapCenter_y
region_r = 2.0
region_name = 'Center'
do_fit = 0
if 'Crab' in source_name:
    region_x = MapCenter_x
    region_y = MapCenter_y
    region_r = CommonPlotFunctions.calibration_radius
    region_name = 'Center'
elif 'PSR_J1907_p0602' in source_name:
    #3HWC J1908+063, 287.05, 6.39 
    region_x = 287.05
    region_y = 6.39
    region_r = 1.2
    region_name = '3HWC'
elif 'PSR_J2021_p4026' in source_name:
    region_x = 305.0200000
    region_y = 40.7572222
    region_r = 1.0
    region_name = 'Center'
MakeExtensionProfile(region_x,region_y,region_r,do_fit,region_name,hist_real_flux_skymap_sum,hist_imposter_flux_skymap_sum,'sum')
MakeSpectrum(region_x,region_y,region_r,region_name,excl_region_x,excl_region_y,excl_region_r)

hist_real_flux_skymap_sum = CommonPlotFunctions.Smooth2DMap(hist_real_flux_skymap_sum,smooth_size_spectroscopy,False)
hist_real_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_reflect,None,[],fig,'RA','Dec','$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFlux_Sum_%s'%(plot_tag))

hist_real_diff_skymap_sum = CommonPlotFunctions.Smooth2DMap(hist_real_diff_skymap_sum,smooth_size_spectroscopy,False)
hist_real_diff_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_real_diff_skymap_reflect,None,[],fig,'RA','Dec','Excess count','SkymapExcess_Sum_%s'%(plot_tag))

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

print ('total_data_expo = %0.1f hrs'%(total_data_expo))
