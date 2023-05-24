
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

n_imposters = 5
if not doImposter:
    n_imposters = 0

plot_tag = source_name
plot_tag += '_'+analysis_method

list_epoch = []
if 'V5' in input_epoch:
    list_epoch += ['V5']
if 'V6' in input_epoch:
    list_epoch += ['V6']

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

def GetFluxCalibration(energy,elev):

    #return 1.

    # The energy threshold needs to be as low as 100 GeV for this method to work.

    str_flux_calibration_el70 = ['2.04e-10', '3.78e-11', '2.41e-11', '1.78e-11', '1.28e-11', '5.13e-12', '3.06e-12', '1.93e-12', '1.32e-12', '8.38e-13', '5.83e-13']
    str_flux_calibration_el60 = ['2.51e-10', '4.15e-11', '2.38e-11', '1.67e-11', '1.19e-11', '4.86e-12', '2.88e-12', '1.76e-12', '1.24e-12', '7.66e-13', '5.32e-13']
    str_flux_calibration_el50 = ['1.22e-09', '8.26e-11', '2.73e-11', '1.42e-11', '8.48e-12', '3.48e-12', '2.04e-12', '1.19e-12', '8.39e-13', '4.99e-13', '3.50e-13']

    flux_calibration_el70 = []
    flux_calibration_el60 = []
    flux_calibration_el50 = []
    for string in str_flux_calibration_el70:
        flux_calibration_el70.append(float(string))
    for string in str_flux_calibration_el60:
        flux_calibration_el60.append(float(string))
    for string in str_flux_calibration_el50:
        flux_calibration_el50.append(float(string))

    xp = [62.,74.,77.]
    fp = [flux_calibration_el50[energy],flux_calibration_el60[energy],flux_calibration_el70[energy]]

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
    print ('data_name = %s'%(data_name))
    for eb in range(0,len(energy_axis)):
        print ('%.2f %.2f %.2f %.2e %.2e %s'%(energy_mean[eb],energy_edge_lo[eb],energy_edge_hi[eb],flux_mean[eb],flux_error[eb],0))

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
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',zorder=3)
    else:
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS')
        real_flux_total_err = pow(real_flux_syst_err*real_flux_syst_err+real_flux_stat_err*real_flux_stat_err,0.5)
        PrintSpectralDataForNaima(energy_axis,real_flux,real_flux_total_err,'VERITAS')

    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    axbig.set_xscale('log')
    axbig.set_yscale('log')
    axbig.legend(loc='best')
    plotname = 'RealSpectrum_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()



def MakeFluxMap(flux_map, data_map, bkgd_map, norm_map, elev_map):

    skymap_bin_size_x = data_map[0].GetXaxis().GetBinCenter(2)-data_map[0].GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = data_map[0].GetYaxis().GetBinCenter(2)-data_map[0].GetYaxis().GetBinCenter(1)
    calibration_radius = CommonPlotFunctions.calibration_radius
    for ebin in range(0,len(energy_bin)-1):
        flux_map[ebin].Reset()
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

                norm_ratio = norm_content/norm_content_max
                norm_weight = 1./(1.+np.exp(-(norm_ratio-0.3)/0.1))
                correction = correction*norm_weight

                stat_data_err = pow(max(data_content,0.),0.5)
                flux_stat_err = max(stat_data_err,1.)/norm_content*correction*pow(energy_bin[ebin]/1e3,2)
                flux_content = (data_content-bkgd_content)/norm_content*correction*pow(energy_bin[ebin]/1e3,2)
                flux_map[ebin].SetBinContent(binx+1,biny+1,flux_content)
                flux_map[ebin].SetBinError(binx+1,biny+1,flux_stat_err)


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
    real_profile_total_err = []
    for ubin in range(0,len(theta2)):
        stat_err = real_profile_stat_err[ubin]
        syst_err = real_profile_syst_err[ubin]
        real_profile_total_err += [pow(stat_err*stat_err+syst_err*syst_err,0.5)]
    real_profile_total_err = np.array(real_profile_total_err)

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

    fig.clf()
    fig.set_figheight(figsize_y*0.8)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    axbig.errorbar(theta2,real_profile,real_profile_stat_err,color='k',marker='s',ls='none',label='ON data')
    #if fit_profile!=0:
    #    axbig.plot(theta2,diffusion_func(theta2,*popt),color='r')
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
    axbig.set_xlabel('angular distance from source [degree]')
    axbig.legend(loc='best')
    #axbig.set_ylim([-0.2e-11, 1.2e-11])
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
InputFile.Close()

hist_elev_skymap = ROOT.TH2D("hist_elev_skymap","",20,MapEdge_left,MapEdge_right,20,MapEdge_lower,MapEdge_upper)
hist_azim_skymap = ROOT.TH2D("hist_azim_skymap","",20,MapEdge_left,MapEdge_right,20,MapEdge_lower,MapEdge_upper)
hist_nsb_skymap = ROOT.TH2D("hist_nsb_skymap","",20,MapEdge_left,MapEdge_right,20,MapEdge_lower,MapEdge_upper)

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
hist_real_significance_skymap = []
for ebin in range(0,len(energy_bin)-1):
    hist_real_flux_skymap += [ROOT.TH2D("hist_real_flux_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_data_skymap += [ROOT.TH2D("hist_real_data_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_bkgd_skymap += [ROOT.TH2D("hist_real_bkgd_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_diff_skymap += [ROOT.TH2D("hist_real_diff_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_significance_skymap += [ROOT.TH2D("hist_real_significance_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]

hist_imposter_data_skymap_sum = []
hist_imposter_bkgd_skymap_sum = []
hist_imposter_diff_skymap_sum = []
hist_imposter_flux_skymap_sum = []
hist_imposter_data_skymap = []
hist_imposter_bkgd_skymap = []
hist_imposter_diff_skymap = []
for imposter in range(0,n_imposters):
    hist_imposter_data_skymap_sum += [ROOT.TH2D("hist_imposter_data_skymap_sum_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_bkgd_skymap_sum += [ROOT.TH2D("hist_imposter_bkgd_skymap_sum_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_diff_skymap_sum += [ROOT.TH2D("hist_imposter_diff_skymap_sum_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_flux_skymap_sum += [ROOT.TH2D("hist_imposter_diff_skymap_sum_I%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap_sublist = []
    hist_imposter_bkgd_skymap_sublist = []
    hist_imposter_diff_skymap_sublist = []
    for ebin in range(0,len(energy_bin)-1):
        hist_imposter_data_skymap_sublist += [ROOT.TH2D("hist_imposter_data_skymap_I%s_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_bkgd_skymap_sublist += [ROOT.TH2D("hist_imposter_bkgd_skymap_I%s_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_diff_skymap_sublist += [ROOT.TH2D("hist_imposter_diff_skymap_I%s_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap += [hist_imposter_data_skymap_sublist]
    hist_imposter_bkgd_skymap += [hist_imposter_bkgd_skymap_sublist]
    hist_imposter_diff_skymap += [hist_imposter_diff_skymap_sublist]


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
                HistName = "Hist_Data_Elev_Skymap"
                FillSkyMapHistogram(InputFile.Get(HistName),hist_elev_skymap)
                HistName = "Hist_Data_Azim_Skymap"
                FillSkyMapHistogram(InputFile.Get(HistName),hist_azim_skymap)
                HistName = "Hist_Data_NSB_Skymap"
                FillSkyMapHistogram(InputFile.Get(HistName),hist_nsb_skymap)
                for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
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
                        HistName = "Hist_Data_Elev_Skymap"
                        FillSkyMapHistogram(InputFile.Get(HistName),hist_elev_skymap)
                        HistName = "Hist_Data_Azim_Skymap"
                        FillSkyMapHistogram(InputFile.Get(HistName),hist_azim_skymap)
                        HistName = "Hist_Data_NSB_Skymap"
                        FillSkyMapHistogram(InputFile.Get(HistName),hist_nsb_skymap)
                        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
                            if effective_area[ebin] < 30000.: continue
                            HistName = "Hist_OnData_SR_Skymap_Sum_ErecS%sto%s"%(int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                            FillSkyMapHistogram(InputFile.Get(HistName),hist_imposter_data_skymap[imposter][ebin])
                            HistName = "Hist_OnData_CR_Skymap_%s_Sum_ErecS%sto%s"%(analysis_method,int(energy_bin[ebin]),int(energy_bin[ebin+1]))
                            FillSkyMapHistogram(InputFile.Get(HistName),hist_imposter_bkgd_skymap[imposter][ebin])
                            hist_imposter_diff_skymap[imposter][ebin].Add(hist_imposter_data_skymap[imposter][ebin])
                            hist_imposter_diff_skymap[imposter][ebin].Add(hist_imposter_bkgd_skymap[imposter][ebin],-1.)
                        InputFile.Close()

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    significance_skymap = CommonPlotFunctions.GetSignificanceMap(hist_real_data_skymap[ebin],hist_real_bkgd_skymap[ebin])
    hist_real_significance_skymap[ebin].Add(significance_skymap)

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_real_diff_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap[ebin])
    CommonPlotFunctions.MatplotlibMap2D(hist_real_diff_skymap_reflect,None,[hist_real_diff_skymap_reflect],fig,'RA','Dec','Excess count','SkymapExcess_E%s_%s'%(ebin,plot_tag))
    hist_real_significance_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_significance_skymap[ebin])
    CommonPlotFunctions.MatplotlibMap2D(hist_real_significance_skymap_reflect,None,[],fig,'RA','Dec','Significance','SkymapSignificance_E%s_%s'%(ebin,plot_tag))

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

hist_real_data_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_data_skymap_sum)
hist_real_bkgd_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_bkgd_skymap_sum)
CommonPlotFunctions.BackgroundSubtractMap(fig,hist_real_data_skymap_sum_reflect,hist_real_bkgd_skymap_sum_reflect,'RA','Dec','Count','SkymapBkgSubtraction_%s'%(plot_tag))

significance_skymap = CommonPlotFunctions.GetSignificanceMap(hist_real_data_skymap_sum,hist_real_bkgd_skymap_sum)
hist_real_significance_skymap_sum.Add(significance_skymap)
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

excl_region_x = MapCenter_x
excl_region_y = MapCenter_y
excl_region_r = 0.0
region_x = MapCenter_x
region_y = MapCenter_y
region_r = 2.0
region_name = 'Center'
do_fit = 0
MakeExtensionProfile(region_x,region_y,region_r,do_fit,region_name,hist_real_flux_skymap_sum,hist_imposter_flux_skymap_sum,'sum')
MakeSpectrum(region_x,region_y,region_r,region_name,excl_region_x,excl_region_y,excl_region_r)

smooth_size_spectroscopy = 0.07
hist_real_flux_skymap_sum = CommonPlotFunctions.Smooth2DMap(hist_real_flux_skymap_sum,smooth_size_spectroscopy,False)
hist_real_diff_skymap_sum = CommonPlotFunctions.Smooth2DMap(hist_real_diff_skymap_sum,smooth_size_spectroscopy,False)

hist_real_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_reflect,None,[],fig,'RA','Dec','$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFlux_Sum_%s'%(plot_tag))
hist_real_diff_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_diff_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_real_diff_skymap_reflect,None,[],fig,'RA','Dec','Excess count','SkymapExcess_Sum_%s'%(plot_tag))
