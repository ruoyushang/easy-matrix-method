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
from scipy import interpolate
import scipy.stats as st
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import cycle
from matplotlib import ticker

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

import CommonPlotFunctions

folder_path = CommonPlotFunctions.folder_path
folder_tag = CommonPlotFunctions.folder_tag
energy_bin = CommonPlotFunctions.energy_bin
n_xoff_bins = CommonPlotFunctions.n_xoff_bins
n_yoff_bins = CommonPlotFunctions.n_yoff_bins

fig, ax = plt.subplots()
figsize_x = 8
figsize_y = 6
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")
np.set_printoptions(precision=4)

measurement_rebin = 1

elev_tag = 'incl'
elev_range = [40.,90.]
#elev_tag = 'sza'
#elev_range = [60.,90.]
#elev_tag = 'lza'
#elev_range = [40.,60.]

total_data_expo = 0.
expo_sum_all_energies = 0.
total_n_measurements = 0.

def GetMeanRMSProfile(array_x, array_y, start_x, end_x, binsize_x):

    nbins_x = int((end_x-start_x)/binsize_x)

    profile_x = []
    profile_mean = []
    profile_rms = []
    for binx in range(0,nbins_x):
        x_low = start_x+binx*binsize_x
        x_up = start_x+(binx+1)*binsize_x
        new_array_y = []
        for entry in range(0,len(array_x)):
            if array_x[entry]<x_low: continue
            if array_x[entry]>x_up: continue
            new_array_y += [array_y[entry]]
        if len(new_array_y)>2:
            new_array_y_mean = np.mean(new_array_y)
            new_array_y_rms =  np.sqrt(np.mean(np.square(new_array_y)))
        else:
            new_array_y_mean = 0.
            new_array_y_rms = 0.
        profile_x += [0.5*(x_low+x_up)]
        profile_mean += [new_array_y_mean]
        profile_rms += [new_array_y_rms]

    return np.array(profile_x), np.array(profile_mean), np.array(profile_rms)


def GetRunInfo(file_path):

    InputFile = ROOT.TFile(file_path)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    total_cr_count = InfoTree.total_cr_count
    data_expo = InfoTree.exposure_hours
    elev_mean = InfoTree.Elev_mean
    azim_mean = InfoTree.Azim_mean
    nsb_mean = InfoTree.NSB_mean
    avg_diff_nsb = InfoTree.avg_diff_nsb
    avg_diff_el = InfoTree.avg_diff_el
    avg_diff_az = InfoTree.avg_diff_az
    InputFile.Close()

    return data_expo, total_cr_count, elev_mean, azim_mean, nsb_mean, avg_diff_el, avg_diff_az, avg_diff_nsb

def Find_s1_on_diagonal_correlation(s0_input,s1_input,e_min,e_max):

    n_eqn = 0
    s0_array = []
    s1_array = []
    for e_idx in range(e_min,e_max):
        n_eqn += len(s1_input[e_idx])
        for entry in range(0,len(s1_input[e_idx])):
            s0_array += [s0_input[e_idx][entry]]
            s1_array += [s1_input[e_idx][entry]]

    n_var = 1
    if n_eqn<n_var:
        print ('Not enough equations.')
        return

    b = np.zeros(n_eqn)
    A = np.zeros((n_eqn,n_var))

    for eqn in range(0,n_eqn):
            A[eqn,0] = s0_array[eqn]
            b[eqn] = s1_array[eqn]

    # Now compute the SVD of  A
    U, sigma, VT = la.svd(A)
    # Make a matrix Sigma of the correct size
    Sigma = np.zeros(A.shape)
    Sigma[:n_var,:n_var] = np.diag(sigma)

    # Now define Sigma_pinv as the "pseudo-"inverse of Sigma, where "pseudo" means "don't divide by zero
    Sigma_pinv = np.zeros(A.shape).T
    Sigma_pinv[:n_var,:n_var] = np.diag(1/sigma[:n_var])

    # Now compute the SVD-based solution for the least-squares problem
    x_svd = VT.T.dot(Sigma_pinv).dot(U.T).dot(b)
    chi2 = la.norm(A.dot(x_svd)-b, 2)/float(n_eqn)

    txt_info = 'chi2 = %0.3e \t,'%(chi2)
    txt_info += 'x0 = %0.3f \t,'%(x_svd[0])
    print (txt_info)

    if np.isnan(x_svd[0]):
        return 0
    else:
        return x_svd[0]

def Find_t11_off_diagonal_correlation(t00_input,t01_input,t10_input,t11_input,e_min,e_max):

    n_eqn = 0
    t00_array = []
    t01_array = []
    t10_array = []
    t11_array = []
    for e_idx in range(e_min,e_max):
        n_eqn += len(t11_input[e_idx])
        for entry in range(0,len(t11_input[e_idx])):
            t00_array += [t00_input[e_idx][entry]]
            t01_array += [t01_input[e_idx][entry]]
            t10_array += [t10_input[e_idx][entry]]
            t11_array += [t11_input[e_idx][entry]]

    n_var = 3
    if n_eqn<n_var:
        print ('Not enough equations.')
        return

    b = np.zeros(n_eqn)
    A = np.zeros((n_eqn,n_var))

    for eqn in range(0,n_eqn):
            A[eqn,0] = t00_array[eqn]
            A[eqn,1] = t01_array[eqn]
            A[eqn,2] = t10_array[eqn]
            b[eqn] = t11_array[eqn]

    # Now compute the SVD of  A
    U, sigma, VT = la.svd(A)
    # Make a matrix Sigma of the correct size
    Sigma = np.zeros(A.shape)
    Sigma[:n_var,:n_var] = np.diag(sigma)

    # Now define Sigma_pinv as the "pseudo-"inverse of Sigma, where "pseudo" means "don't divide by zero
    Sigma_pinv = np.zeros(A.shape).T
    Sigma_pinv[:n_var,:n_var] = np.diag(1/sigma[:n_var])

    # Now compute the SVD-based solution for the least-squares problem
    x_svd = VT.T.dot(Sigma_pinv).dot(U.T).dot(b)
    chi2 = la.norm(A.dot(x_svd)-b, 2)/float(n_eqn)

    txt_info = 'chi2 = %0.3e \t,'%(chi2)
    txt_info += 'x00 = %0.3f \t,'%(x_svd[0])
    txt_info += 'x01 = %0.3f \t,'%(x_svd[1])
    txt_info += 'x10 = %0.3f \t,'%(x_svd[2])
    print (txt_info)

    if np.isnan(x_svd[0]):
        return 0, 0, 0
    else:
        return x_svd[0], x_svd[1], x_svd[2]

def GetGammaCounts(file_path,ebin):

    effective_area = ROOT.std.vector("double")(20)
    s0_truth = ROOT.std.vector("double")(20)
    s1_truth = ROOT.std.vector("double")(20)
    s2_truth = ROOT.std.vector("double")(20)
    s0_recon = ROOT.std.vector("double")(20)
    s1_recon = ROOT.std.vector("double")(20)
    s2_recon = ROOT.std.vector("double")(20)
    t00_truth = ROOT.std.vector("double")(20)
    t01_truth = ROOT.std.vector("double")(20)
    t10_truth = ROOT.std.vector("double")(20)
    t11_truth = ROOT.std.vector("double")(20)
    t02_truth = ROOT.std.vector("double")(20)
    t20_truth = ROOT.std.vector("double")(20)
    t00_recon = ROOT.std.vector("double")(20)
    t01_recon = ROOT.std.vector("double")(20)
    t10_recon = ROOT.std.vector("double")(20)
    t11_recon = ROOT.std.vector("double")(20)
    t02_recon = ROOT.std.vector("double")(20)
    t20_recon = ROOT.std.vector("double")(20)
    data_count = ROOT.std.vector("double")(20)
    ratio_bkgd_count = ROOT.std.vector("double")(20)
    regression_bkgd_count = ROOT.std.vector("double")(20)
    init_perturbation_bkgd_count = ROOT.std.vector("double")(20)
    perturbation_bkgd_count = ROOT.std.vector("double")(20)
    combined_bkgd_count = ROOT.std.vector("double")(20)

    InputFile = ROOT.TFile(file_path)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.SetBranchAddress('effective_area',ROOT.AddressOf(effective_area))
    InfoTree.SetBranchAddress('s0_truth',ROOT.AddressOf(s0_truth))
    InfoTree.SetBranchAddress('s1_truth',ROOT.AddressOf(s1_truth))
    InfoTree.SetBranchAddress('s2_truth',ROOT.AddressOf(s2_truth))
    InfoTree.SetBranchAddress('s0_recon',ROOT.AddressOf(s0_recon))
    InfoTree.SetBranchAddress('s1_recon',ROOT.AddressOf(s1_recon))
    InfoTree.SetBranchAddress('s2_recon',ROOT.AddressOf(s2_recon))
    InfoTree.SetBranchAddress('t00_truth',ROOT.AddressOf(t00_truth))
    InfoTree.SetBranchAddress('t01_truth',ROOT.AddressOf(t01_truth))
    InfoTree.SetBranchAddress('t10_truth',ROOT.AddressOf(t10_truth))
    InfoTree.SetBranchAddress('t11_truth',ROOT.AddressOf(t11_truth))
    InfoTree.SetBranchAddress('t02_truth',ROOT.AddressOf(t02_truth))
    InfoTree.SetBranchAddress('t20_truth',ROOT.AddressOf(t20_truth))
    InfoTree.SetBranchAddress('t00_recon',ROOT.AddressOf(t00_recon))
    InfoTree.SetBranchAddress('t01_recon',ROOT.AddressOf(t01_recon))
    InfoTree.SetBranchAddress('t10_recon',ROOT.AddressOf(t10_recon))
    InfoTree.SetBranchAddress('t11_recon',ROOT.AddressOf(t11_recon))
    InfoTree.SetBranchAddress('t02_recon',ROOT.AddressOf(t02_recon))
    InfoTree.SetBranchAddress('t20_recon',ROOT.AddressOf(t20_recon))
    InfoTree.SetBranchAddress('data_count',ROOT.AddressOf(data_count))
    InfoTree.SetBranchAddress('ratio_bkgd_count',ROOT.AddressOf(ratio_bkgd_count))
    InfoTree.SetBranchAddress('regression_bkgd_count',ROOT.AddressOf(regression_bkgd_count))
    InfoTree.SetBranchAddress('init_perturbation_bkgd_count',ROOT.AddressOf(init_perturbation_bkgd_count))
    InfoTree.SetBranchAddress('perturbation_bkgd_count',ROOT.AddressOf(perturbation_bkgd_count))
    InfoTree.SetBranchAddress('combined_bkgd_count',ROOT.AddressOf(combined_bkgd_count))
    InfoTree.GetEntry(0)
    InputFile.Close()

    return effective_area[ebin], data_count[ebin], ratio_bkgd_count[ebin], regression_bkgd_count[ebin], init_perturbation_bkgd_count[ebin], perturbation_bkgd_count[ebin], combined_bkgd_count[ebin], t00_truth[ebin], t01_truth[ebin], t10_truth[ebin], t11_truth[ebin], t02_truth[ebin], t20_truth[ebin], t00_recon[ebin], t01_recon[ebin], t10_recon[ebin], t11_recon[ebin], t02_recon[ebin], t20_recon[ebin], s0_truth[ebin], s1_truth[ebin], s2_truth[ebin], s0_recon[ebin], s1_recon[ebin], s2_recon[ebin]

def MakeMultipleFitPlot(ax_input,Hists,legends,title_x,title_y):

    hist_xdata = []
    hist_ydata = []
    hist_error = []
    for entry in range(0,len(Hists)):
        xdata = []
        ydata = []
        error = []
        for binx in range(0,Hists[entry].GetNbinsX()):
            xdata += [Hists[entry].GetBinCenter(binx+1)]
            ydata += [Hists[entry].GetBinContent(binx+1)]
            error += [max(1.,Hists[entry].GetBinError(binx+1))]
        hist_xdata += [xdata]
        hist_ydata += [ydata]
        hist_error += [error]

    cycol = cycle('brgcmk')
    for entry in range(0,len(Hists)):
        next_color = next(cycol)
        bin_width = hist_xdata[entry][1]-hist_xdata[entry][0]
        ax_input.bar(hist_xdata[entry], hist_ydata[entry], color=next_color, ls='solid', label='%s'%(legends[entry]),width=bin_width,alpha=0.5)

    ax_input.legend(loc='best')
    ax_input.set_xlabel(title_x)
    ax_input.set_ylabel(title_y)
    return(ax_input)


ONOFF_tag_sample = 'OFF'
source_of_interest = ''
#source_of_interest = '1ES0229'
#source_of_interest = 'H1426'
#source_of_interest = 'UrsaMinor'

sample_list = []

sample_list += ['OJ287_V6_OFF']
sample_list += ['Segue1_V6_OFF']
sample_list += ['3C264_V6_OFF']
sample_list += ['3C273_V6_OFF']
sample_list += ['PG1553_V6_OFF']
sample_list += ['1ES0229_V6_OFF']
sample_list += ['PKS1424_V6_OFF']
sample_list += ['RGB_J0710_p591_V6_OFF'] # north 
sample_list += ['UrsaMajorII_V6_OFF'] # north
sample_list += ['UrsaMinor_V6_OFF'] # north
sample_list += ['H1426_V6_OFF'] # north
sample_list += ['NGC1275_V6_OFF'] # north
sample_list += ['Draco_V6_OFF'] # north
sample_list += ['BLLac_V6_OFF'] # north
sample_list += ['1ES0502_V6_OFF'] # north
sample_list += ['M82_V6_OFF'] # north
sample_list += ['1ES0414_V6_OFF'] # north
sample_list += ['1ES1011_V6_OFF'] # north
sample_list += ['1ES0647_V6_OFF']
sample_list += ['OJ287_V5_OFF']
sample_list += ['Segue1_V5_OFF']
sample_list += ['3C264_V5_OFF']
sample_list += ['3C273_V5_OFF']
sample_list += ['PG1553_V5_OFF']
sample_list += ['1ES0229_V5_OFF']
sample_list += ['PKS1424_V5_OFF']
sample_list += ['RGB_J0710_p591_V5_OFF'] # north 
sample_list += ['UrsaMajorII_V5_OFF'] # north
sample_list += ['UrsaMinor_V5_OFF'] # north
sample_list += ['H1426_V5_OFF'] # north
sample_list += ['NGC1275_V5_OFF'] # north
sample_list += ['Draco_V5_OFF'] # north
sample_list += ['BLLac_V5_OFF'] # north
sample_list += ['1ES0502_V5_OFF'] # north
sample_list += ['M82_V5_OFF'] # north
sample_list += ['1ES0414_V5_OFF'] # north
sample_list += ['1ES1011_V5_OFF'] # north
sample_list += ['1ES0647_V5_OFF']

sample_list += ['CrabNebula_elev_80_90_V6_OFF']
sample_list += ['CrabNebula_elev_70_80_V6_OFF']
sample_list += ['CrabNebula_elev_60_70_V6_OFF']
sample_list += ['CrabNebula_elev_50_60_V6_OFF']
sample_list += ['CrabNebula_elev_40_50_V6_OFF']
sample_list += ['CrabNebula_elev_80_90_V5_OFF']
sample_list += ['CrabNebula_elev_70_80_V5_OFF']
sample_list += ['CrabNebula_elev_60_70_V5_OFF']
sample_list += ['CrabNebula_elev_50_60_V5_OFF']
sample_list += ['CrabNebula_elev_40_50_V5_OFF']

nbins = 21
hist_limit = 30.0
Hist_SystErrDist_Ratio = []
Hist_SystErrDist_Regression = []
Hist_SystErrDist_Init_Perturbation = []
Hist_SystErrDist_Perturbation = []
Hist_SystErrDist_Combined = []
array_cr_count = []
array_elev_mean = []
array_azim_mean = []
array_nsb_mean = []
array_elev_diff = []
array_azim_diff = []
array_nsb_diff = []
array_s0_truth = []
array_s1_truth = []
array_s2_truth = []
array_s0_recon = []
array_s1_recon = []
array_s2_recon = []
array_t00_truth = []
array_t01_truth = []
array_t10_truth = []
array_t11_truth = []
array_t02_truth = []
array_t20_truth = []
array_t00_recon = []
array_t01_recon = []
array_t10_recon = []
array_t11_recon = []
array_t02_recon = []
array_t20_recon = []
array_syst_err_ratio = []
array_syst_err_regression = []
array_syst_err_init_perturbation = []
array_syst_err_perturbation = []
array_syst_err_combined = []
array_rebin_stat_err = []
array_rebin_syst_err_ratio = []
array_rebin_syst_err_regression = []
array_rebin_syst_err_init_perturbation = []
array_rebin_syst_err_perturbation = []
array_rebin_syst_err_combined = []
for energy_idx in range(0,len(energy_bin)-1):
    Hist_SystErrDist_Ratio += [ROOT.TH1D("Hist_SystErrDist_Ratio_E%s"%(energy_idx),"",nbins,-hist_limit,hist_limit)]
    Hist_SystErrDist_Regression += [ROOT.TH1D("Hist_SystErrDist_Regression_E%s"%(energy_idx),"",nbins,-hist_limit,hist_limit)]
    Hist_SystErrDist_Init_Perturbation += [ROOT.TH1D("Hist_SystErrDist_Init_Perturbation_E%s"%(energy_idx),"",nbins,-hist_limit,hist_limit)]
    Hist_SystErrDist_Perturbation += [ROOT.TH1D("Hist_SystErrDist_Perturbation_E%s"%(energy_idx),"",nbins,-hist_limit,hist_limit)]
    Hist_SystErrDist_Combined += [ROOT.TH1D("Hist_SystErrDist_Combined_E%s"%(energy_idx),"",nbins,-hist_limit,hist_limit)]

txt_warning = ''
for energy_idx in range(0,len(energy_bin)-1):
    array_per_energy_cr_count = []
    array_per_energy_elev_mean = []
    array_per_energy_azim_mean = []
    array_per_energy_nsb_mean = []
    array_per_energy_elev_diff = []
    array_per_energy_azim_diff = []
    array_per_energy_s0_truth = []
    array_per_energy_s1_truth = []
    array_per_energy_s2_truth = []
    array_per_energy_s0_recon = []
    array_per_energy_s1_recon = []
    array_per_energy_s2_recon = []
    array_per_energy_t00_truth = []
    array_per_energy_t01_truth = []
    array_per_energy_t10_truth = []
    array_per_energy_t11_truth = []
    array_per_energy_t02_truth = []
    array_per_energy_t20_truth = []
    array_per_energy_t00_recon = []
    array_per_energy_t01_recon = []
    array_per_energy_t10_recon = []
    array_per_energy_t11_recon = []
    array_per_energy_t02_recon = []
    array_per_energy_t20_recon = []
    array_per_energy_nsb_diff = []
    array_syst_err_per_energy_ratio = []
    array_syst_err_per_energy_regression = []
    array_syst_err_per_energy_init_perturbation = []
    array_syst_err_per_energy_perturbation = []
    array_syst_err_per_energy_combined = []
    array_rebin_stat_err_per_energy = []
    array_rebin_syst_err_per_energy_ratio = []
    array_rebin_syst_err_per_energy_regression = []
    array_rebin_syst_err_per_energy_init_perturbation = []
    array_rebin_syst_err_per_energy_perturbation = []
    array_rebin_syst_err_per_energy_combined = []
    for src in range(0,len(sample_list)):
        if not source_of_interest=='':
            if not source_of_interest in sample_list[src]: continue
        xoff_idx = 0
        yoff_idx = 0
        for xoff_idx in range(0,n_xoff_bins):
            for yoff_idx in range(0,n_yoff_bins):
                n_groups = 0
                file_exists = True
                while file_exists:
                    SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_G%d_X%d_Y%d.root"%(folder_path,sample_list[src],n_groups,xoff_idx,yoff_idx)
                    if os.path.exists(SourceFilePath):
                        n_groups += 1
                        #print ('Read file: %s'%(SourceFilePath))
                        #print ('file exists.')
                    else:
                        file_exists = False
                        print ('Read file: %s'%(SourceFilePath))
                        print ('file does not exist.')
                total_data_truth = 0.
                total_ratio_bkgd = 0.
                total_regression_bkgd = 0.
                total_init_perturbation_bkgd = 0.
                total_perturbation_bkgd = 0.
                total_combined_bkgd = 0.
                n_rebin = 0
                for group in range(0,n_groups):
                    SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_G%d_X%d_Y%d.root"%(folder_path,sample_list[src],group,xoff_idx,yoff_idx)
                    if not os.path.exists(SourceFilePath): continue
                    print ('Read file: %s'%(SourceFilePath))
                    eff_area, data_truth, ratio_bkgd, regression_bkgd, init_perturbation_bkgd, perturbation_bkgd, combined_bkgd, t00_truth, t01_truth, t10_truth, t11_truth, t02_truth, t20_truth, t00_recon, t01_recon, t10_recon, t11_recon, t02_recon, t20_recon, s0_truth, s1_truth, s2_truth, s0_recon, s1_recon, s2_recon = GetGammaCounts(SourceFilePath,energy_idx)
                    data_expo, total_cr_count, elev_mean, azim_mean, nsb_mean, avg_diff_el, avg_diff_az, avg_diff_nsb = GetRunInfo(SourceFilePath)
                    if elev_mean<elev_range[0] or elev_mean>elev_range[1]: continue
                    if xoff_idx==0 and yoff_idx==0:
                        expo_sum_all_energies += data_expo
                        if energy_idx==1:
                            total_data_expo += data_expo
                    if data_truth==0.: continue
                    total_data_truth += data_truth
                    total_ratio_bkgd += ratio_bkgd
                    total_regression_bkgd += regression_bkgd
                    total_init_perturbation_bkgd += init_perturbation_bkgd
                    total_perturbation_bkgd += perturbation_bkgd
                    total_combined_bkgd += combined_bkgd
                    array_per_energy_cr_count += [total_cr_count]
                    array_per_energy_elev_mean += [elev_mean]
                    array_per_energy_azim_mean += [azim_mean]
                    array_per_energy_nsb_mean += [nsb_mean]
                    array_per_energy_elev_diff += [avg_diff_el]
                    array_per_energy_azim_diff += [avg_diff_az]
                    array_per_energy_nsb_diff += [avg_diff_nsb]
                    array_syst_err_per_energy_ratio += [-(ratio_bkgd-data_truth)/pow(data_truth,0.5)]
                    array_syst_err_per_energy_regression += [-(regression_bkgd-data_truth)/pow(data_truth,0.5)]
                    array_syst_err_per_energy_init_perturbation += [-(init_perturbation_bkgd-data_truth)/pow(data_truth,0.5)]
                    array_syst_err_per_energy_perturbation += [-(perturbation_bkgd-data_truth)/pow(data_truth,0.5)]
                    array_syst_err_per_energy_combined += [-(combined_bkgd-data_truth)/pow(data_truth,0.5)]
                    perturbation_error = -(perturbation_bkgd-data_truth)/pow(data_truth,0.5)
                    if abs(perturbation_error)>15.:
                        txt_warning += '++++++++++++++++++++++++++++++++++++ \n'
                        txt_warning += '!!! Error > 15 sigma !!!! \n'
                        txt_warning += 'sample: %s \n '%(sample_list[src])
                        txt_warning += 'elev_mean = %s \n'%(elev_mean)
                        txt_warning += 'energy index = %s \n'%(energy_idx)
                        txt_warning += 'error = %0.1f \n '%(perturbation_error)
                        txt_warning += 'truth = %0.1f \n '%(data_truth)
                    if s0_truth==0.: continue
                    #if avg_diff_el<0.: continue
                    #if avg_diff_el>0.: continue
                    #if avg_diff_az<0.: continue
                    #if avg_diff_az>0.: continue
                    #if avg_diff_nsb<0.: continue
                    #if avg_diff_nsb>0.: continue
                    #if azim_mean>90.: continue
                    array_per_energy_s0_truth += [s0_truth]
                    array_per_energy_s1_truth += [s1_truth]
                    array_per_energy_s2_truth += [s2_truth]
                    array_per_energy_s0_recon += [s0_recon]
                    array_per_energy_s1_recon += [s1_recon]
                    array_per_energy_s2_recon += [s2_recon]
                    array_per_energy_t00_truth += [t00_truth]
                    array_per_energy_t01_truth += [t01_truth]
                    array_per_energy_t10_truth += [t10_truth]
                    array_per_energy_t11_truth += [t11_truth]
                    array_per_energy_t02_truth += [t02_truth]
                    array_per_energy_t20_truth += [t20_truth]
                    array_per_energy_t00_recon += [t00_recon]
                    array_per_energy_t01_recon += [t01_recon]
                    array_per_energy_t10_recon += [t10_recon]
                    array_per_energy_t11_recon += [t11_recon]
                    array_per_energy_t02_recon += [t02_recon]
                    array_per_energy_t20_recon += [t20_recon]
            n_rebin += 1
            if n_rebin == measurement_rebin:
                if total_data_truth>0.:
                    Hist_SystErrDist_Ratio[energy_idx].Fill((total_ratio_bkgd-total_data_truth)/pow(total_data_truth,0.5))
                    Hist_SystErrDist_Regression[energy_idx].Fill((total_regression_bkgd-total_data_truth)/pow(total_data_truth,0.5))
                    Hist_SystErrDist_Init_Perturbation[energy_idx].Fill((total_init_perturbation_bkgd-total_data_truth)/pow(total_data_truth,0.5))
                    Hist_SystErrDist_Perturbation[energy_idx].Fill((total_perturbation_bkgd-total_data_truth)/pow(total_data_truth,0.5))
                    Hist_SystErrDist_Combined[energy_idx].Fill((total_combined_bkgd-total_data_truth)/pow(total_data_truth,0.5))
                    array_rebin_stat_err_per_energy += [1.]
                    array_rebin_syst_err_per_energy_ratio += [(total_ratio_bkgd-total_data_truth)/pow(total_data_truth,0.5)]
                    array_rebin_syst_err_per_energy_regression += [(total_regression_bkgd-total_data_truth)/pow(total_data_truth,0.5)]
                    array_rebin_syst_err_per_energy_init_perturbation += [(total_init_perturbation_bkgd-total_data_truth)/pow(total_data_truth,0.5)]
                    array_rebin_syst_err_per_energy_perturbation += [(total_perturbation_bkgd-total_data_truth)/pow(total_data_truth,0.5)]
                    array_rebin_syst_err_per_energy_combined += [(total_combined_bkgd-total_data_truth)/pow(total_data_truth,0.5)]
                total_data_truth = 0.
                total_ratio_bkgd = 0.
                total_regression_bkgd = 0.
                total_init_perturbation_bkgd = 0.
                total_perturbation_bkgd = 0.
                total_combined_bkgd = 0.
                n_rebin = 0
                total_n_measurements += 1.
    array_cr_count += [array_per_energy_cr_count]
    array_elev_mean += [array_per_energy_elev_mean]
    array_azim_mean += [array_per_energy_azim_mean]
    array_nsb_mean += [array_per_energy_nsb_mean]
    array_elev_diff += [array_per_energy_elev_diff]
    array_azim_diff += [array_per_energy_azim_diff]
    array_nsb_diff += [array_per_energy_nsb_diff]
    array_s0_truth += [array_per_energy_s0_truth]
    array_s1_truth += [array_per_energy_s1_truth]
    array_s2_truth += [array_per_energy_s2_truth]
    array_s0_recon += [array_per_energy_s0_recon]
    array_s1_recon += [array_per_energy_s1_recon]
    array_s2_recon += [array_per_energy_s2_recon]
    array_t00_truth += [array_per_energy_t00_truth]
    array_t01_truth += [array_per_energy_t01_truth]
    array_t10_truth += [array_per_energy_t10_truth]
    array_t11_truth += [array_per_energy_t11_truth]
    array_t02_truth += [array_per_energy_t02_truth]
    array_t20_truth += [array_per_energy_t20_truth]
    array_t00_recon += [array_per_energy_t00_recon]
    array_t01_recon += [array_per_energy_t01_recon]
    array_t10_recon += [array_per_energy_t10_recon]
    array_t11_recon += [array_per_energy_t11_recon]
    array_t02_recon += [array_per_energy_t02_recon]
    array_t20_recon += [array_per_energy_t20_recon]
    array_syst_err_ratio += [array_syst_err_per_energy_ratio]
    array_syst_err_regression += [array_syst_err_per_energy_regression]
    array_syst_err_init_perturbation += [array_syst_err_per_energy_init_perturbation]
    array_syst_err_perturbation += [array_syst_err_per_energy_perturbation]
    array_syst_err_combined += [array_syst_err_per_energy_combined]
    array_rebin_stat_err += [array_rebin_stat_err_per_energy]
    array_rebin_syst_err_ratio += [array_rebin_syst_err_per_energy_ratio]
    array_rebin_syst_err_regression += [array_rebin_syst_err_per_energy_regression]
    array_rebin_syst_err_init_perturbation += [array_rebin_syst_err_per_energy_init_perturbation]
    array_rebin_syst_err_perturbation += [array_rebin_syst_err_per_energy_perturbation]
    array_rebin_syst_err_combined += [array_rebin_syst_err_per_energy_combined]



for energy_idx in range(0,len(energy_bin)-1):

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_s0_truth[energy_idx],array_s1_truth[energy_idx],color='k',alpha=0.5)
    axbig.scatter(array_s0_recon[energy_idx],array_s1_recon[energy_idx],color='r',alpha=0.5)
    axbig.set_xlabel('$s_{0}$')
    axbig.set_ylabel('$s_{1}$')
    fig.savefig("output_plots/s0_vs_s1_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_s0_truth[energy_idx],array_t00_truth[energy_idx],color='k',alpha=0.5)
    axbig.set_xlabel('$s_{0}$')
    axbig.set_ylabel('$t_{00}$')
    fig.savefig("output_plots/s0_vs_t00_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_t00_truth[energy_idx],array_t11_truth[energy_idx],color='k',alpha=0.5)
    axbig.scatter(array_t00_recon[energy_idx],array_t11_recon[energy_idx],color='r',alpha=0.5)
    axbig.set_xlabel('$t_{00}$')
    axbig.set_ylabel('$t_{11}$')
    fig.savefig("output_plots/t00_vs_t11_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_t01_truth[energy_idx],array_t11_truth[energy_idx],color='k',alpha=0.5)
    axbig.scatter(array_t01_recon[energy_idx],array_t11_recon[energy_idx],color='r',alpha=0.5)
    axbig.set_xlabel('$t_{01}$')
    axbig.set_ylabel('$t_{11}$')
    fig.savefig("output_plots/t01_vs_t11_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_t10_truth[energy_idx],array_t11_truth[energy_idx],color='k',alpha=0.5)
    axbig.scatter(array_t10_recon[energy_idx],array_t11_recon[energy_idx],color='r',alpha=0.5)
    axbig.set_xlabel('$t_{10}$')
    axbig.set_ylabel('$t_{11}$')
    fig.savefig("output_plots/t10_vs_t11_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_cr_count[energy_idx],array_syst_err_perturbation[energy_idx],color='k',alpha=0.5)
    axbig.set_xlabel('Total CR count')
    axbig.set_ylabel('perturbation method $\epsilon$')
    fig.savefig("output_plots/CR_count_vs_Error_Perturbation_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_elev_diff[energy_idx],array_syst_err_perturbation[energy_idx],color='k',alpha=0.5)
    axbig.set_xlabel('Elevation difference [deg]')
    axbig.set_ylabel('perturbation method $\epsilon$')
    fig.savefig("output_plots/ElevDiff_vs_Error_Perturbation_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_azim_diff[energy_idx],array_syst_err_perturbation[energy_idx],color='k',alpha=0.5)
    axbig.set_xlabel('Azimuth difference [deg]')
    axbig.set_ylabel('perturbation method $\epsilon$')
    fig.savefig("output_plots/AzimDiff_vs_Error_Perturbation_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_nsb_diff[energy_idx],array_syst_err_perturbation[energy_idx],color='k',alpha=0.5)
    axbig.set_xlabel('NSB difference')
    axbig.set_ylabel('perturbation method $\epsilon$')
    fig.savefig("output_plots/NSBDiff_vs_Error_Perturbation_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_cr_count[energy_idx],array_syst_err_regression[energy_idx],color='k',alpha=0.5)
    axbig.set_xlabel('Total CR count')
    axbig.set_ylabel('regression method $\epsilon$')
    fig.savefig("output_plots/CR_count_vs_Error_Regression_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_elev_diff[energy_idx],array_syst_err_regression[energy_idx],color='k',alpha=0.5)
    axbig.set_xlabel('Elevation difference [deg]')
    axbig.set_ylabel('regression method $\epsilon$')
    fig.savefig("output_plots/ElevDiff_vs_Error_Regression_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_azim_diff[energy_idx],array_syst_err_regression[energy_idx],color='k',alpha=0.5)
    axbig.set_xlabel('Azimuth difference [deg]')
    axbig.set_ylabel('regression method $\epsilon$')
    fig.savefig("output_plots/AzimDiff_vs_Error_Regression_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_nsb_diff[energy_idx],array_syst_err_regression[energy_idx],color='k',alpha=0.5)
    axbig.set_xlabel('NSB difference')
    axbig.set_ylabel('regression method $\epsilon$')
    fig.savefig("output_plots/NSBDiff_vs_Error_Regression_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()


    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    x_start = 55.
    x_end = 90.
    x_delta = 5.
    baseline_xaxis = np.linspace(x_start,x_end,100)
    baseline_yaxis = [0. for i in range(0,len(baseline_xaxis))]
    axbig.scatter(array_elev_mean[energy_idx],array_syst_err_ratio[energy_idx],color='k',alpha=0.1)
    axbig.scatter(array_elev_mean[energy_idx],array_syst_err_perturbation[energy_idx],color='r',alpha=0.1)
    axbig.plot(baseline_xaxis, baseline_yaxis, color='gray', ls='dashed')
    xaxis_elev, yaxis_ratio_mean, yaxis_ratio_rms = GetMeanRMSProfile(array_elev_mean[energy_idx], array_syst_err_ratio[energy_idx], x_start, x_end, x_delta)
    axbig.errorbar(xaxis_elev-0.2*x_delta,yaxis_ratio_mean,yaxis_ratio_rms,color='k',marker='_',ls='none',linewidth=2)
    xaxis_elev, yaxis_perturbation_mean, yaxis_perturbation_rms = GetMeanRMSProfile(array_elev_mean[energy_idx], array_syst_err_perturbation[energy_idx], x_start, x_end, x_delta)
    axbig.errorbar(xaxis_elev+0.2*x_delta,yaxis_perturbation_mean,yaxis_perturbation_rms,color='r',marker='_',ls='none',linewidth=2)
    axbig.set_xlabel('Elevation [deg]')
    axbig.set_ylabel('Error $\epsilon$')
    fig.savefig("output_plots/Elev_vs_Error_Perturbation_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    x_start = 0.
    x_end = 180.
    x_delta = 20.
    baseline_xaxis = np.linspace(x_start,x_end,100)
    baseline_yaxis = [0. for i in range(0,len(baseline_xaxis))]
    axbig.scatter(array_azim_mean[energy_idx],array_syst_err_ratio[energy_idx],color='k',alpha=0.1)
    axbig.scatter(array_azim_mean[energy_idx],array_syst_err_perturbation[energy_idx],color='r',alpha=0.1)
    axbig.plot(baseline_xaxis, baseline_yaxis, color='gray', ls='dashed')
    xaxis_azim, yaxis_ratio_mean, yaxis_ratio_rms = GetMeanRMSProfile(array_azim_mean[energy_idx], array_syst_err_ratio[energy_idx], x_start, x_end, x_delta)
    axbig.errorbar(xaxis_azim-0.2*x_delta,yaxis_ratio_mean,yaxis_ratio_rms,color='k',marker='_',ls='none',linewidth=2)
    xaxis_azim, yaxis_perturbation_mean, yaxis_perturbation_rms = GetMeanRMSProfile(array_azim_mean[energy_idx], array_syst_err_perturbation[energy_idx], x_start, x_end, x_delta)
    axbig.errorbar(xaxis_azim+0.2*x_delta,yaxis_perturbation_mean,yaxis_perturbation_rms,color='r',marker='_',ls='none',linewidth=2)
    axbig.set_xlabel('Azimuth [deg]')
    axbig.set_ylabel('Error $\epsilon$')
    fig.savefig("output_plots/Azim_vs_Error_Perturbation_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    x_start = 2.5
    x_end = 9.5
    x_delta = 1.
    baseline_xaxis = np.linspace(x_start,x_end,100)
    baseline_yaxis = [0. for i in range(0,len(baseline_xaxis))]
    axbig.scatter(array_nsb_mean[energy_idx],array_syst_err_ratio[energy_idx],color='k',alpha=0.1)
    axbig.scatter(array_nsb_mean[energy_idx],array_syst_err_perturbation[energy_idx],color='r',alpha=0.1)
    axbig.plot(baseline_xaxis, baseline_yaxis, color='gray', ls='dashed')
    xaxis_nsb, yaxis_ratio_mean, yaxis_ratio_rms = GetMeanRMSProfile(array_nsb_mean[energy_idx], array_syst_err_ratio[energy_idx], x_start, x_end, x_delta)
    axbig.errorbar(xaxis_nsb-0.2*x_delta,yaxis_ratio_mean,yaxis_ratio_rms,color='k',marker='_',ls='none',linewidth=2)
    xaxis_nsb, yaxis_perturbation_mean, yaxis_perturbation_rms = GetMeanRMSProfile(array_nsb_mean[energy_idx], array_syst_err_perturbation[energy_idx], x_start, x_end, x_delta)
    axbig.errorbar(xaxis_nsb+0.2*x_delta,yaxis_perturbation_mean,yaxis_perturbation_rms,color='r',marker='_',ls='none',linewidth=2)
    axbig.set_xlabel('NSB')
    axbig.set_ylabel('Error $\epsilon$')
    fig.savefig("output_plots/NSB_vs_Error_Perturbation_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    x_start = 2.5
    x_end = 9.5
    x_delta = 1.
    baseline_xaxis = np.linspace(x_start,x_end,100)
    baseline_yaxis = [0. for i in range(0,len(baseline_xaxis))]
    axbig.scatter(array_nsb_mean[energy_idx],array_syst_err_ratio[energy_idx],color='k',alpha=0.1)
    axbig.scatter(array_nsb_mean[energy_idx],array_syst_err_regression[energy_idx],color='r',alpha=0.1)
    axbig.plot(baseline_xaxis, baseline_yaxis, color='gray', ls='dashed')
    xaxis_nsb, yaxis_ratio_mean, yaxis_ratio_rms = GetMeanRMSProfile(array_nsb_mean[energy_idx], array_syst_err_ratio[energy_idx], x_start, x_end, x_delta)
    axbig.errorbar(xaxis_nsb-0.2*x_delta,yaxis_ratio_mean,yaxis_ratio_rms,color='k',marker='_',ls='none',linewidth=2)
    xaxis_nsb, yaxis_regression_mean, yaxis_regression_rms = GetMeanRMSProfile(array_nsb_mean[energy_idx], array_syst_err_regression[energy_idx], x_start, x_end, x_delta)
    axbig.errorbar(xaxis_nsb+0.2*x_delta,yaxis_regression_mean,yaxis_regression_rms,color='r',marker='_',ls='none',linewidth=2)
    axbig.set_xlabel('NSB')
    axbig.set_ylabel('Error $\epsilon$')
    fig.savefig("output_plots/NSB_vs_Error_Regression_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()



for energy_idx in range(0,len(energy_bin)-1):
    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_syst_err_regression[energy_idx],array_syst_err_perturbation[energy_idx],color='b',alpha=0.5)
    axbig.set_xlabel('Regression method $\epsilon$')
    axbig.set_ylabel('Perturbation method $\epsilon$')
    fig.savefig("output_plots/Regression_vs_Perturbation_Correlation_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

for energy_idx in range(0,len(energy_bin)-1):
    fig.clf()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    axbig = fig.add_subplot()
    axbig.scatter(array_syst_err_init_perturbation[energy_idx],array_syst_err_perturbation[energy_idx],color='b',alpha=0.5)
    axbig.set_xlabel('Off-diagonal method $\epsilon$')
    axbig.set_ylabel('Perturbation method $\epsilon$')
    fig.savefig("output_plots/Init_vs_Perturbation_Correlation_E%s_%s.png"%(energy_idx,folder_tag))
    axbig.remove()

for energy_idx in range(0,len(energy_bin)-1):
    Hists = []
    legends = []
    Hists += [Hist_SystErrDist_Ratio[energy_idx]]
    legends += ['Simple Scaling']
    Hists += [Hist_SystErrDist_Perturbation[energy_idx]]
    legends += ['Low-rank Perturbation']

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    MakeMultipleFitPlot(axbig,Hists,legends,'relative error $\epsilon$','number of entries')
    fig.savefig("output_plots/SystErrDist_E%s_Perturbation_%s.png"%(energy_idx,folder_tag))
    axbig.remove()



energy_dep_stat_err = []
energy_dep_syst_err_ratio_mean = []
energy_dep_syst_err_regression_mean = []
energy_dep_syst_err_init_perturbation_mean = []
energy_dep_syst_err_perturbation_mean = []
energy_dep_syst_err_combined_mean = []
energy_dep_syst_err_ratio_rms = []
energy_dep_syst_err_regression_rms = []
energy_dep_syst_err_init_perturbation_rms = []
energy_dep_syst_err_perturbation_rms = []
energy_dep_syst_err_combined_rms = []
for energy_idx in range(0,len(energy_bin)-1):
    array_stat_err = np.mean(array_rebin_stat_err[energy_idx])
    array_syst_err_ratio_mean = np.mean(array_rebin_syst_err_ratio[energy_idx])
    array_syst_err_regression_mean = np.mean(array_rebin_syst_err_regression[energy_idx])
    array_syst_err_init_perturbation_mean = np.mean(array_rebin_syst_err_init_perturbation[energy_idx])
    array_syst_err_perturbation_mean = np.mean(array_rebin_syst_err_perturbation[energy_idx])
    array_syst_err_combined_mean = np.mean(array_rebin_syst_err_combined[energy_idx])
    array_syst_err_ratio_rms = np.sqrt(np.mean(np.square(array_rebin_syst_err_ratio[energy_idx])))
    array_syst_err_regression_rms = np.sqrt(np.mean(np.square(array_rebin_syst_err_regression[energy_idx])))
    array_syst_err_init_perturbation_rms = np.sqrt(np.mean(np.square(array_rebin_syst_err_init_perturbation[energy_idx])))
    array_syst_err_perturbation_rms = np.sqrt(np.mean(np.square(array_rebin_syst_err_perturbation[energy_idx])))
    array_syst_err_combined_rms = np.sqrt(np.mean(np.square(array_rebin_syst_err_combined[energy_idx])))

    if math.isnan(array_stat_err):
        array_stat_err = 1.
    if math.isnan(array_syst_err_ratio_mean):
        array_syst_err_ratio_mean = 1.
    if math.isnan(array_syst_err_regression_mean):
        array_syst_err_regression_mean = 1.
    if math.isnan(array_syst_err_init_perturbation_mean):
        array_syst_err_init_perturbation_mean = 1.
    if math.isnan(array_syst_err_perturbation_mean):
        array_syst_err_perturbation_mean = 1.
    if math.isnan(array_syst_err_combined_mean):
        array_syst_err_combined_mean = 1.
    if math.isnan(array_syst_err_ratio_rms):
        array_syst_err_ratio_rms = 1.
    if math.isnan(array_syst_err_regression_rms):
        array_syst_err_regression_rms = 1.
    if math.isnan(array_syst_err_init_perturbation_rms):
        array_syst_err_init_perturbation_rms = 1.
    if math.isnan(array_syst_err_perturbation_rms):
        array_syst_err_perturbation_rms = 1.
    if math.isnan(array_syst_err_combined_rms):
        array_syst_err_combined_rms = 1.

    energy_dep_stat_err += [array_stat_err]
    energy_dep_syst_err_ratio_mean += [array_syst_err_ratio_mean]
    energy_dep_syst_err_regression_mean += [array_syst_err_regression_mean]
    energy_dep_syst_err_init_perturbation_mean += [array_syst_err_init_perturbation_mean]
    energy_dep_syst_err_perturbation_mean += [array_syst_err_perturbation_mean]
    energy_dep_syst_err_combined_mean += [array_syst_err_combined_mean]
    energy_dep_syst_err_ratio_rms += [array_syst_err_ratio_rms]
    energy_dep_syst_err_regression_rms += [array_syst_err_regression_rms]
    energy_dep_syst_err_init_perturbation_rms += [array_syst_err_init_perturbation_rms]
    energy_dep_syst_err_perturbation_rms += [array_syst_err_perturbation_rms]
    energy_dep_syst_err_combined_rms += [array_syst_err_combined_rms]

    print ('================================================================================================')
    print ('Energy = %s'%(energy_bin[energy_idx]))
    print ('stat. error of data = %0.3f'%(array_stat_err))
    print ('mean of syst. error of simple scaling method = %0.3f'%(array_syst_err_ratio_mean))
    print ('rms of syst. error of simple scaling method = %0.3f'%(array_syst_err_ratio_rms))
    print ('mean of syst. error of regression method = %0.3f'%(array_syst_err_regression_mean))
    print ('rms of syst. error of regression method = %0.3f'%(array_syst_err_regression_rms))
    print ('mean of syst. error of init_perturbation method = %0.3f'%(array_syst_err_init_perturbation_mean))
    print ('rms of syst. error of init_perturbation method = %0.3f'%(array_syst_err_init_perturbation_rms))
    print ('mean of syst. error of perturbation method = %0.3f'%(array_syst_err_perturbation_mean))
    print ('rms of syst. error of perturbation method = %0.3f'%(array_syst_err_perturbation_rms))
    print ('mean of syst. error of combined method = %0.3f'%(array_syst_err_combined_mean))
    print ('rms of syst. error of combined method = %0.3f'%(array_syst_err_combined_rms))

print ('================================================================================================')
txt_string = 'double method_ratio_rms[N_energy_bins] =              {'
for energy_idx in range(0,len(energy_bin)-1):
    txt_string += '%0.3f,'%(energy_dep_syst_err_ratio_rms[energy_idx])
txt_string += '};'
print (txt_string)
txt_string = 'double method_init_perturbation_rms[N_energy_bins] =  {'
for energy_idx in range(0,len(energy_bin)-1):
    txt_string += '%0.3f,'%(energy_dep_syst_err_init_perturbation_rms[energy_idx])
txt_string += '};'
print (txt_string)
txt_string = 'double method_perturbation_rms[N_energy_bins] =       {'
for energy_idx in range(0,len(energy_bin)-1):
    txt_string += '%0.3f,'%(energy_dep_syst_err_perturbation_rms[energy_idx])
txt_string += '};'
print (txt_string)
txt_string = 'double method_regression_rms[N_energy_bins] =         {'
for energy_idx in range(0,len(energy_bin)-1):
    txt_string += '%0.3f,'%(energy_dep_syst_err_regression_rms[energy_idx])
txt_string += '};'
print (txt_string)
txt_string = 'double method_combined_rms[N_energy_bins] =           {'
for energy_idx in range(0,len(energy_bin)-1):
    txt_string += '%0.3f,'%(energy_dep_syst_err_combined_rms[energy_idx])
txt_string += '};'
print (txt_string)
print ('================================================================================================')
for energy_idx in range(0,len(energy_bin)-1):
    txt_string = 'Energy %0.1f GeV\t, '%(energy_bin[energy_idx])
    txt_string += '%0.3f \t, '%(energy_dep_syst_err_ratio_rms[energy_idx])
    txt_string += '%0.3f \t, '%(energy_dep_syst_err_init_perturbation_rms[energy_idx])
    txt_string += '%0.3f \t, '%(energy_dep_syst_err_perturbation_rms[energy_idx])
    txt_string += '%0.3f \t, '%(energy_dep_syst_err_regression_rms[energy_idx])
    txt_string += '%0.3f \t, '%(energy_dep_syst_err_combined_rms[energy_idx])
    print (txt_string)

print ('================================================================================================')
print (txt_warning)


txt_info_x00 = 'double coefficient_t11xt00_%s[N_energy_bins] = {'%(elev_tag)
txt_info_x01 = 'double coefficient_t11xt01_%s[N_energy_bins] = {'%(elev_tag)
txt_info_x10 = 'double coefficient_t11xt10_%s[N_energy_bins] = {'%(elev_tag)
print ('================================================================================================')
for energy_idx in range(0,len(energy_bin)-1):
    x00, x01, x10 = Find_t11_off_diagonal_correlation(array_t00_truth,array_t01_truth,array_t10_truth,array_t11_truth,energy_idx,energy_idx+1)
    txt_info_x00 += '%0.3f'%(x00)
    txt_info_x01 += '%0.3f'%(x01)
    txt_info_x10 += '%0.3f'%(x10)
    if energy_idx<len(energy_bin)-2:
        txt_info_x00 += ','
        txt_info_x01 += ','
        txt_info_x10 += ','
txt_info_x00 += '};'
txt_info_x01 += '};'
txt_info_x10 += '};'
print (txt_info_x00)
print (txt_info_x01)
print (txt_info_x10)
print ('================================================================================================')

txt_info_x0 = 'double coefficient_s1xs0_%s[N_energy_bins] = {'%(elev_tag)
print ('================================================================================================')
for energy_idx in range(0,len(energy_bin)-1):
    x0 = Find_s1_on_diagonal_correlation(array_s0_truth,array_s1_truth,energy_idx,energy_idx+1)
    txt_info_x0 += '%0.3f'%(x0)
    if energy_idx<len(energy_bin)-2:
        txt_info_x0 += ','
txt_info_x0 += '};'
print (txt_info_x0)
print ('================================================================================================')

print ('total_data_expo = %0.1f hrs'%(total_data_expo))
print ('avg expo per measurement = %0.1f'%(expo_sum_all_energies/total_n_measurements))
