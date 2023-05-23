
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

folder_path = CommonPlotFunctions.folder_path

epoch_idx = 0
SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_%s_G1.root"%(folder_path,source_name,list_epoch[0],isON)
if os.path.exists(SourceFilePath):
    epoch_idx = 0
else:
    epoch_idx = 1
InputFile = ROOT.TFile("/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_%s_G1.root"%(folder_path,source_name,list_epoch[epoch_idx],isON))
HistName = "Hist_OnData_SR_Skymap_Sum_ErecS100to200"
nbins_x = InputFile.Get(HistName).GetNbinsX()
nbins_y = InputFile.Get(HistName).GetNbinsY()
MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
InputFile.Close()

hist_real_data_skymap_sum = ROOT.TH2D("hist_real_data_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_bkgd_skymap_sum = ROOT.TH2D("hist_real_bkgd_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_significance_skymap_sum = ROOT.TH2D("hist_real_significance_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_data_skymap = []
hist_real_bkgd_skymap = []
hist_real_excess_skymap = []
hist_real_significance_skymap = []
for ebin in range(0,len(energy_bin)-1):
    hist_real_data_skymap += [ROOT.TH2D("hist_real_data_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_bkgd_skymap += [ROOT.TH2D("hist_real_bkgd_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_excess_skymap += [ROOT.TH2D("hist_real_excess_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_significance_skymap += [ROOT.TH2D("hist_real_significance_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]

hist_elev_skymap = ROOT.TH2D("hist_elev_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_azim_skymap = ROOT.TH2D("hist_azim_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_nsb_skymap = ROOT.TH2D("hist_nsb_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)

effective_area = ROOT.std.vector("double")(20)

for epoch in list_epoch:
    n_groups = 1
    file_exists = True
    while file_exists:
        SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_%s_G%d.root"%(folder_path,source_name,epoch,isON,n_groups)
        print ('Read file: %s'%(SourceFilePath))
        if os.path.exists(SourceFilePath):
            n_groups += 1
            print ('file exists.')
        else:
            file_exists = False
            print ('file does not exist.')
    
    for group in range(1,n_groups):
        InputFile = ROOT.TFile("/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_%s_%s_%s_G%d.root"%(folder_path,source_name,epoch,isON,group))
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
            hist_real_excess_skymap[ebin].Add(hist_real_data_skymap[ebin])
            hist_real_excess_skymap[ebin].Add(hist_real_bkgd_skymap[ebin],-1.)
        InputFile.Close()

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    smooth_size_spectroscopy = 0.07
    hist_real_excess_skymap[ebin] = CommonPlotFunctions.Smooth2DMap(hist_real_excess_skymap[ebin],smooth_size_spectroscopy,False)
    hist_real_data_skymap[ebin] = CommonPlotFunctions.Smooth2DMap(hist_real_data_skymap[ebin],smooth_size_spectroscopy,False)
    hist_real_bkgd_skymap[ebin] = CommonPlotFunctions.Smooth2DMap(hist_real_bkgd_skymap[ebin],smooth_size_spectroscopy,False)
    significance_skymap = CommonPlotFunctions.GetSignificanceMap(hist_real_data_skymap[ebin],hist_real_bkgd_skymap[ebin])
    hist_real_significance_skymap[ebin].Add(significance_skymap)

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_real_excess_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_excess_skymap[ebin])
    CommonPlotFunctions.MatplotlibMap2D(hist_real_excess_skymap_reflect,None,[hist_real_excess_skymap_reflect],fig,'RA','Dec','Excess count','SkymapExcess_E%s_%s'%(ebin,plot_tag))
    hist_real_significance_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_significance_skymap[ebin])
    CommonPlotFunctions.MatplotlibMap2D(hist_real_significance_skymap_reflect,None,[],fig,'RA','Dec','Significance','SkymapSignificance_E%s_%s'%(ebin,plot_tag))

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_real_data_skymap_sum.Add(hist_real_data_skymap[ebin])
    hist_real_bkgd_skymap_sum.Add(hist_real_bkgd_skymap[ebin])

hist_real_data_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_data_skymap_sum)
hist_real_bkgd_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_bkgd_skymap_sum)
CommonPlotFunctions.BackgroundSubtractMap(fig,hist_real_data_skymap_sum_reflect,hist_real_bkgd_skymap_sum_reflect,'RA','Dec','Count','SkymapBkgSubtraction_%s'%(plot_tag))

significance_skymap = CommonPlotFunctions.GetSignificanceMap(hist_real_data_skymap_sum,hist_real_bkgd_skymap_sum)
hist_real_significance_skymap_sum.Add(significance_skymap)
hist_real_significance_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_significance_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_real_significance_skymap_reflect,None,[],fig,'RA','Dec','Significance','SkymapSignificance_Sum_%s'%(plot_tag))

hist_elev_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_elev_skymap)
CommonPlotFunctions.MatplotlibMap2D(hist_elev_skymap_reflect,hist_real_significance_skymap_reflect,[],fig,'RA','Dec','Elevation [deg]','SkymapElev')
hist_azim_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_azim_skymap)
CommonPlotFunctions.MatplotlibMap2D(hist_azim_skymap_reflect,hist_real_significance_skymap_reflect,[],fig,'RA','Dec','Azimuth [deg]','SkymapAzim')
hist_nsb_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_nsb_skymap)
CommonPlotFunctions.MatplotlibMap2D(hist_nsb_skymap_reflect,hist_real_significance_skymap_reflect,[],fig,'RA','Dec','NSB','SkymapNSB')

