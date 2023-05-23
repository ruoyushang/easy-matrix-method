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
energy_bin = CommonPlotFunctions.energy_bin

fig, ax = plt.subplots()
figsize_x = 6.5
figsize_y = 4.5
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")
np.set_printoptions(precision=4)

def GetGammaCounts(file_path,ebin):

    effective_area = ROOT.std.vector("double")(20)
    data_count = ROOT.std.vector("double")(20)
    ratio_bkgd_count = ROOT.std.vector("double")(20)
    regression_bkgd_count = ROOT.std.vector("double")(20)
    perturbation_bkgd_count = ROOT.std.vector("double")(20)
    combined_bkgd_count = ROOT.std.vector("double")(20)

    InputFile = ROOT.TFile(file_path)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.SetBranchAddress('effective_area',ROOT.AddressOf(effective_area))
    InfoTree.SetBranchAddress('data_count',ROOT.AddressOf(data_count))
    InfoTree.SetBranchAddress('ratio_bkgd_count',ROOT.AddressOf(ratio_bkgd_count))
    InfoTree.SetBranchAddress('regression_bkgd_count',ROOT.AddressOf(regression_bkgd_count))
    InfoTree.SetBranchAddress('perturbation_bkgd_count',ROOT.AddressOf(perturbation_bkgd_count))
    InfoTree.SetBranchAddress('combined_bkgd_count',ROOT.AddressOf(combined_bkgd_count))
    InfoTree.GetEntry(0)
    data_expo = InfoTree.exposure_hours
    InputFile.Close()

    return effective_area[ebin], data_count[ebin], ratio_bkgd_count[ebin], regression_bkgd_count[ebin], perturbation_bkgd_count[ebin], combined_bkgd_count[ebin]

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

nbins = 21
hist_limit = 0.2
Hist_SystErrDist_Ratio = []
Hist_SystErrDist_Regression = []
Hist_SystErrDist_Perturbation = []
Hist_SystErrDist_Combined = []
array_syst_err_ratio = []
array_syst_err_regression = []
array_syst_err_perturbation = []
array_syst_err_combined = []
for energy_idx in range(0,len(energy_bin)-1):
    Hist_SystErrDist_Ratio += [ROOT.TH1D("Hist_SystErrDist_Ratio_E%s"%(energy_idx),"",nbins,-hist_limit,hist_limit)]
    Hist_SystErrDist_Regression += [ROOT.TH1D("Hist_SystErrDist_Regression_E%s"%(energy_idx),"",nbins,-hist_limit,hist_limit)]
    Hist_SystErrDist_Perturbation += [ROOT.TH1D("Hist_SystErrDist_Perturbation_E%s"%(energy_idx),"",nbins,-hist_limit,hist_limit)]
    Hist_SystErrDist_Combined += [ROOT.TH1D("Hist_SystErrDist_Combined_E%s"%(energy_idx),"",nbins,-hist_limit,hist_limit)]

for energy_idx in range(0,len(energy_bin)-1):
    array_syst_err_per_energy_ratio = []
    array_syst_err_per_energy_regression = []
    array_syst_err_per_energy_perturbation = []
    array_syst_err_per_energy_combined = []
    for src in range(0,len(sample_list)):
        n_groups = 1
        file_exists = True
        while file_exists:
            SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_G%d.root"%(folder_path,sample_list[src],n_groups)
            print ('Read file: %s'%(SourceFilePath))
            if os.path.exists(SourceFilePath):
                n_groups += 1
                print ('file exists.')
            else:
                file_exists = False
                print ('file does not exist.')
        for group in range(1,n_groups):
            SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_G%d.root"%(folder_path,sample_list[src],group)
            eff_area, data_truth, ratio_bkgd, regression_bkgd, perturbation_bkgd, combined_bkgd = GetGammaCounts(SourceFilePath,energy_idx)
            if eff_area<10000.: continue
            if data_truth<10.: continue
            print ('data_truth = %s'%(data_truth))
            Hist_SystErrDist_Ratio[energy_idx].Fill((ratio_bkgd-data_truth)/data_truth)
            Hist_SystErrDist_Regression[energy_idx].Fill((regression_bkgd-data_truth)/data_truth)
            Hist_SystErrDist_Perturbation[energy_idx].Fill((perturbation_bkgd-data_truth)/data_truth)
            Hist_SystErrDist_Combined[energy_idx].Fill((combined_bkgd-data_truth)/data_truth)
            array_syst_err_per_energy_ratio += [(ratio_bkgd-data_truth)/data_truth]
            array_syst_err_per_energy_regression += [(regression_bkgd-data_truth)/data_truth]
            array_syst_err_per_energy_perturbation += [(perturbation_bkgd-data_truth)/data_truth]
            array_syst_err_per_energy_combined += [(combined_bkgd-data_truth)/data_truth]
    array_syst_err_ratio += [array_syst_err_per_energy_ratio]
    array_syst_err_regression += [array_syst_err_per_energy_regression]
    array_syst_err_perturbation += [array_syst_err_per_energy_perturbation]
    array_syst_err_combined += [array_syst_err_per_energy_combined]

for energy_idx in range(0,len(energy_bin)-1):
    array_syst_err_ratio_mean = np.mean(array_syst_err_ratio[energy_idx])
    array_syst_err_regression_mean = np.mean(array_syst_err_regression[energy_idx])
    array_syst_err_perturbation_mean = np.mean(array_syst_err_perturbation[energy_idx])
    array_syst_err_combined_mean = np.mean(array_syst_err_combined[energy_idx])
    array_syst_err_ratio_rms = np.sqrt(np.mean(np.square(array_syst_err_ratio[energy_idx])))
    array_syst_err_regression_rms = np.sqrt(np.mean(np.square(array_syst_err_regression[energy_idx])))
    array_syst_err_perturbation_rms = np.sqrt(np.mean(np.square(array_syst_err_perturbation[energy_idx])))
    array_syst_err_combined_rms = np.sqrt(np.mean(np.square(array_syst_err_combined[energy_idx])))
    print ('================================================================================================')
    print ('Energy = %s'%(energy_bin[energy_idx]))
    print ('rms of syst. error of simple scaling method = %0.3f'%(array_syst_err_ratio_rms))
    print ('rms of syst. error of regression method = %0.3f'%(array_syst_err_regression_rms))
    print ('rms of syst. error of perturbation method = %0.3f'%(array_syst_err_perturbation_rms))
    print ('rms of syst. error of combined method = %0.3f'%(array_syst_err_combined_rms))

#for energy_idx in range(0,len(energy_bin)-1):
#    Hists = []
#    legends = []
#    Hists += [Hist_SystErrDist_Ratio[energy_idx]]
#    legends += ['Simple Scaling']
#    Hists += [Hist_SystErrDist_Combined[energy_idx]]
#    legends += ['Combined']
#
#    fig.clf()
#    figsize_x = 6.5
#    figsize_y = 4.5
#    fig.set_figheight(figsize_y)
#    fig.set_figwidth(figsize_x)
#    axbig = fig.add_subplot()
#    MakeMultipleFitPlot(axbig,Hists,legends,'relative error $\epsilon$','number of entries')
#    fig.savefig("output_plots/SystErrDist_E%s_Combined.png"%(energy_idx))
#    axbig.remove()

#for energy_idx in range(0,len(energy_bin)-1):
#    Hists = []
#    legends = []
#    Hists += [Hist_SystErrDist_Ratio[energy_idx]]
#    legends += ['Simple Scaling']
#    Hists += [Hist_SystErrDist_Perturbation[energy_idx]]
#    legends += ['Low-rank Perturbation']
#
#    fig.clf()
#    figsize_x = 6.5
#    figsize_y = 4.5
#    fig.set_figheight(figsize_y)
#    fig.set_figwidth(figsize_x)
#    axbig = fig.add_subplot()
#    MakeMultipleFitPlot(axbig,Hists,legends,'relative error $\epsilon$','number of entries')
#    fig.savefig("output_plots/SystErrDist_E%s_Perturbation.png"%(energy_idx))
#    axbig.remove()
#
for energy_idx in range(0,len(energy_bin)-1):
    Hists = []
    legends = []
    Hists += [Hist_SystErrDist_Ratio[energy_idx]]
    legends += ['Simple Scaling']
    Hists += [Hist_SystErrDist_Regression[energy_idx]]
    legends += ['Linear Regression']

    fig.clf()
    figsize_x = 6.5
    figsize_y = 4.5
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    MakeMultipleFitPlot(axbig,Hists,legends,'relative error $\epsilon$','number of entries')
    fig.savefig("output_plots/SystErrDist_E%s_Regression.png"%(energy_idx))
    axbig.remove()

