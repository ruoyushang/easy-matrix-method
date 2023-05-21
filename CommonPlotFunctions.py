
import os
import sys,ROOT
import array
import math
import csv
from array import *
from ROOT import *
from astropy import units as my_unit
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
import astropy.utils as utils
from astropy.nddata import Cutout2D
from scipy import special
from scipy import interpolate
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import cycle
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from matplotlib import ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import NullFormatter
from operator import itemgetter, attrgetter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits import mplot3d
import healpy as hp
from spectral_cube import SpectralCube

# Great examples of matplotlib plots: https://atmamani.github.io/cheatsheets/matplotlib/matplotlib_2/

folder_path = 'output_nuclear'

#analysis_method = 'FoV'
#analysis_method = 'Ratio'
analysis_method = 'Regression'
#analysis_method = 'Perturbation'

energy_bin = [100.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.]

def reflectXaxis(hist):

    # taken from VPlotAnasumHistograms.cpp
	
    # temporary histogram
    new_hist_name = "%s_REFLECTED"%(hist.GetName())
    if ROOT.gDirectory.FindObject(new_hist_name):
        hT = ROOT.gDirectory.Get(new_hist_name)
    else:
        hT = ROOT.TH2D( new_hist_name, "", hist.GetNbinsX(), -1.*hist.GetXaxis().GetXmax(), -1.*hist.GetXaxis().GetXmin(), hist.GetNbinsY(), hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax() )
    hT.SetStats( 0 )
    hT.SetXTitle( hist.GetXaxis().GetTitle() )
    hT.SetYTitle( hist.GetYaxis().GetTitle() )
	
    for binx in range(1,hist.GetNbinsX()+1):
        for biny in range(1,hist.GetNbinsX()+1):
            hT.SetBinContent( hist.GetNbinsX() + 1 - binx, biny, hist.GetBinContent( binx, biny ) )
            hT.SetBinError( hist.GetNbinsX() + 1 - binx, biny, hist.GetBinError( binx, biny ) )
    return hT

def GetGammaSourceInfo():

    other_stars = []
    other_stars_type = []
    other_star_coord = []

    return other_stars, other_stars_type, other_star_coord

    near_source_cut = 0.1

    drawBrightStar = False
    drawPulsar = True
    drawSNR = True
    drawFermi = True
    drawHAWC = False
    drawTeV = False

    if drawBrightStar:
        star_name, star_ra, star_dec = ReadBrightStarListFromFile()
        for src in range(0,len(star_name)):
            src_ra = star_ra[src]
            src_dec = star_dec[src]
            if doGalacticCoord:
                src_ra, src_dec = ConvertRaDecToGalactic(src_ra,src_dec)
            other_stars += [star_name[src]]
            other_stars_type += ['Star']
            other_star_coord += [[src_ra,src_dec,0.]]

    if drawSNR:
        target_snr_name, target_snr_ra, target_snr_dec, target_snr_size = ReadSNRTargetListFromCSVFile()
        for src in range(0,len(target_snr_name)):
            gamma_source_name = target_snr_name[src]
            gamma_source_ra = target_snr_ra[src]
            gamma_source_dec = target_snr_dec[src]
            gamma_source_size = target_snr_size[src]
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['SNR']
                other_star_coord += [[gamma_source_ra,gamma_source_dec,gamma_source_size]]

    if drawPulsar:
        target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age = ReadATNFTargetListFromFile('ATNF_pulsar_full_list.txt')
        for src in range(0,len(target_psr_name)):
            gamma_source_name = target_psr_name[src]
            gamma_source_ra = target_psr_ra[src]
            gamma_source_dec = target_psr_dec[src]
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                if target_psr_age[src]/1000.>1e6:
                        other_stars_type += ['MSP']
                else:
                        other_stars_type += ['PSR']
                other_star_coord += [[gamma_source_ra,gamma_source_dec,0.]]

    if drawFermi:
        fermi_name, fermi_ra, fermi_dec = ReadFermiCatelog()
        for src in range(0,len(fermi_name)):
            gamma_source_name = fermi_name[src]
            gamma_source_ra = fermi_ra[src]
            gamma_source_dec = fermi_dec[src]
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['Fermi']
                other_star_coord += [[gamma_source_ra,gamma_source_dec,0.]]

    if drawHAWC:
        target_hwc_name, target_hwc_ra, target_hwc_dec = ReadHAWCTargetListFromFile('Cat_3HWC.txt')
        for src in range(0,len(target_hwc_name)):
            gamma_source_name = target_hwc_name[src]
            gamma_source_ra = target_hwc_ra[src]
            gamma_source_dec = target_hwc_dec[src]
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['HAWC']
                other_star_coord += [[gamma_source_ra,gamma_source_dec,0.]]

    if drawTeV:
        inputFile = open('TeVCat_RaDec_w_Names.txt')
        for line in inputFile:
            gamma_source_name = line.split(',')[0]
            gamma_source_ra = float(line.split(',')[1])
            gamma_source_dec = float(line.split(',')[2])
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source and not '%' in gamma_source_name:
                other_stars += [gamma_source_name]
                other_stars_type += ['TeV']
                other_star_coord += [[gamma_source_ra,gamma_source_dec,0.]]

    return other_stars, other_stars_type, other_star_coord


def MatplotlibMap2D(hist_map,hist_tone,hist_contour,fig,label_x,label_y,label_z,plotname,roi_x=0.,roi_y=0.,roi_r=0.,rotation_angle=0.):

    print ('Making plot %s...'%(plotname))

    colormap = 'coolwarm'

    map_nbins_x = hist_map.GetNbinsX()
    map_nbins_y = hist_map.GetNbinsY()
    MapEdge_left = hist_map.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_map.GetXaxis().GetBinLowEdge(hist_map.GetNbinsX()+1)
    MapEdge_lower = hist_map.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_map.GetYaxis().GetBinLowEdge(hist_map.GetNbinsY()+1)
    MapCenter_x = (MapEdge_left+MapEdge_right)/2.
    MapCenter_y = (MapEdge_lower+MapEdge_upper)/2.

    deg_per_bin = (MapEdge_right-MapEdge_left)/map_nbins_x
    nbins_per_deg = map_nbins_x/(MapEdge_right-MapEdge_left)
    x_axis = np.linspace(MapEdge_left,MapEdge_right,map_nbins_x)
    prelim_x_axis_sparse = np.linspace(MapEdge_left,MapEdge_right,6)
    x_axis_sparse = []
    for i in prelim_x_axis_sparse:
        if int(i)>MapEdge_left and int(i)<MapEdge_right:
            x_axis_sparse += [int(i)]
    x_axis_reflect = ["{:6.1f}".format(-1.*i) for i in x_axis_sparse]
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,map_nbins_y)

    grid_z = np.zeros((map_nbins_y, map_nbins_x))
    max_z = 0.
    min_z = 0.
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            hist_bin_x = xbin+1
            hist_bin_y = ybin+1
            if hist_bin_x<1: continue
            if hist_bin_y<1: continue
            if hist_bin_x>hist_map.GetNbinsX(): continue
            if hist_bin_y>hist_map.GetNbinsY(): continue
            grid_z[ybin,xbin] = hist_map.GetBinContent(hist_bin_x,hist_bin_y)
            if max_z<hist_map.GetBinContent(hist_bin_x,hist_bin_y):
                max_z = hist_map.GetBinContent(hist_bin_x,hist_bin_y)
            if min_z>hist_map.GetBinContent(hist_bin_x,hist_bin_y):
                min_z = hist_map.GetBinContent(hist_bin_x,hist_bin_y)

    if not hist_tone==None:
        max_zscore = hist_tone.GetMaximum()
        mid_class_z = 0.
        mid_class_z_bins = 0.
        low_class_z = 0.
        low_class_z_bins = 0.
        high_class_z = 0.
        high_class_z_bins = 0.
        for ybin in range(0,len(y_axis)):
            for xbin in range(0,len(x_axis)):
                #hist_bin_x = hist_map.GetXaxis().FindBin(x_axis[xbin])
                #hist_bin_y = hist_map.GetYaxis().FindBin(y_axis[ybin])
                hist_bin_x = xbin+1
                hist_bin_y = ybin+1
                if hist_bin_x<1: continue
                if hist_bin_y<1: continue
                if hist_bin_x>hist_map.GetNbinsX(): continue
                if hist_bin_y>hist_map.GetNbinsY(): continue
                zscore = hist_tone.GetBinContent(hist_bin_x,hist_bin_y)
                if zscore>1.-0.5 and zscore<1.+0.5:
                    mid_class_z += hist_map.GetBinContent(hist_bin_x,hist_bin_y)
                    mid_class_z_bins += 1.
                if zscore>0.-0.5 and zscore<0.+0.5:
                    low_class_z += hist_map.GetBinContent(hist_bin_x,hist_bin_y)
                    low_class_z_bins += 1.
                if zscore>max_zscore-0.5:
                    high_class_z += hist_map.GetBinContent(hist_bin_x,hist_bin_y)
                    high_class_z_bins += 1.
        if mid_class_z_bins>0.:
            mid_class_z = mid_class_z/mid_class_z_bins
        if low_class_z_bins>0.:
            low_class_z = low_class_z/low_class_z_bins
        if high_class_z_bins>0.:
            high_class_z = high_class_z/high_class_z_bins
        min_z = low_class_z-2.*(mid_class_z-low_class_z)
        max_z = max(high_class_z*1.1,mid_class_z+2.*(mid_class_z-low_class_z))

    list_levels = []
    list_levels += [np.arange(0.5, 1.0, 0.3)]
    list_levels += [np.arange(0.5, 1.0, 0.3)]
    list_levels += [np.arange(0.5, 1.0, 0.3)]
    list_grid_contour = []
    if not hist_contour==None:
        for ctr in range(0,len(hist_contour)):
            if max_z==0.: continue
            if hist_contour[ctr].GetMaximum()==0.: continue
            grid_contour = np.zeros((hist_contour[ctr].GetNbinsY(),hist_contour[ctr].GetNbinsX()))
            #max_z_contour = list_levels[ctr][1]
            #min_z_contour = list_levels[ctr][0]
            #delta_z = list_levels[ctr][2]
            max_z_contour = 1.0*hist_contour[ctr].GetMaximum()
            min_z_contour = 0.7*hist_contour[ctr].GetMaximum()
            delta_z = 0.3*hist_contour[ctr].GetMaximum()
            list_levels[ctr] = np.arange(min_z_contour*max_z/max_z_contour, max_z_contour*max_z/max_z_contour, delta_z*max_z/max_z_contour)
            contour_x_axis = np.linspace(MapEdge_left,MapEdge_right,hist_contour[ctr].GetNbinsX())
            contour_y_axis = np.linspace(MapEdge_lower,MapEdge_upper,hist_contour[ctr].GetNbinsY())
            for ybin in range(0,len(contour_y_axis)):
                for xbin in range(0,len(contour_x_axis)):
                    hist_bin_x = xbin+1
                    hist_bin_y = ybin+1
                    if hist_bin_x<1: continue
                    if hist_bin_y<1: continue
                    if hist_bin_x>hist_contour[ctr].GetNbinsX(): continue
                    if hist_bin_y>hist_contour[ctr].GetNbinsY(): continue
                    grid_contour[ybin,xbin] = hist_contour[ctr].GetBinContent(hist_bin_x,hist_bin_y)
                    if label_z!='Significance' and max_z_contour>0.:
                        grid_contour[ybin,xbin] = hist_contour[ctr].GetBinContent(hist_bin_x,hist_bin_y)*max_z/max_z_contour
            list_grid_contour += [grid_contour]

    other_stars, other_star_type, other_star_coord = GetGammaSourceInfo() 

    other_star_labels = []
    other_star_types = []
    other_star_markers = []
    star_range = 0.7*(MapEdge_right-MapEdge_left)/2.
    source_ra = (MapEdge_left+MapEdge_right)/2.
    source_dec = (MapEdge_lower+MapEdge_upper)/2.
    n_stars = 0
    for star in range(0,len(other_stars)):
        if doGalacticCoord:
            if abs(source_ra+other_star_coord[star][0])>star_range: continue
            if abs(source_dec-other_star_coord[star][1])>0.5*star_range: continue
        else:
            #if pow(source_ra+other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
            if abs(source_ra+other_star_coord[star][0])>star_range: continue
            if abs(source_dec-other_star_coord[star][1])>star_range: continue
        #if '#' in other_stars[star]: continue
        other_star_markers += [[-other_star_coord[star][0],other_star_coord[star][1],other_star_coord[star][2]]]
        other_star_labels += ['(%s) %s'%(n_stars,other_stars[star])]
        other_star_types += [other_star_type[star]]
        n_stars += 1

    fig.clf()
    figsize_x = 14
    figsize_y = 7
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    label_x = 'RA'
    label_y = 'Dec'
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    if label_z=='Significance':
        max_z = 5.
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=-max_z,vmax=max_z,zorder=0)
    else:
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=min_z,vmax=max_z,zorder=0)

    list_colors = ['orange','lime','deepskyblue']
    list_styles = ['solid','solid','solid']
    for ctr in range(0,len(list_grid_contour)):
        axbig.contour(list_grid_contour[len(list_grid_contour)-1-ctr], list_levels[len(list_grid_contour)-1-ctr], linestyles=list_styles[len(list_grid_contour)-1-ctr], colors=list_colors[len(list_grid_contour)-1-ctr], extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=1)
    favorite_color = 'red'
    for star in range(0,len(other_star_markers)):
        marker_size = 60
        if other_star_types[star]=='PSR':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c=favorite_color, marker='^', label=other_star_labels[star])
        if other_star_types[star]=='SNR':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c=favorite_color, marker='+', label=other_star_labels[star])
            mycircle = plt.Circle( (other_star_markers[star][0], other_star_markers[star][1]), other_star_markers[star][2], fill = False, color=favorite_color)
            axbig.add_patch(mycircle)
        if other_star_types[star]=='HAWC':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='violet', marker='+', label=other_star_labels[star])
        if other_star_types[star]=='Fermi':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='cyan', marker='+', label=other_star_labels[star])
        if other_star_types[star]=='MSP':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='tomato', marker='+', label=other_star_labels[star])
        if other_star_types[star]=='TeV':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='lime', marker='+', label=other_star_labels[star])
        if other_star_types[star]=='Star':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='k', marker='o', label=other_star_labels[star])
        text_offset_x = 0.07
        text_offset_y = 0.07

    divider = make_axes_locatable(axbig)
    cax = divider.append_axes("bottom", size="5%", pad=0.7)
    cbar = fig.colorbar(im,orientation="horizontal",cax=cax)
    cbar.set_label(label_z)
    axbig.set_xticks(x_axis_sparse)
    axbig.set_xticklabels(x_axis_reflect)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')

    if 'SkymapZscore_Sum' in plotname:
        if not doGalacticCoord:
            if roi_r>0.:
                mycircle = plt.Circle((-roi_x, roi_y), roi_r, color='b', fill=False)
                axbig.add_patch(mycircle)
        axbig.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0, fontsize=7)
        fig.savefig("output_plots/%s_legend.png"%(plotname),bbox_inches='tight')
    axbig.remove()

    if 'SkymapZscore_Sum' in plotname:
        for star in range(0,len(other_star_markers)):
            print ('%s, (%0.2f,%0.2f)'%(other_star_labels[star],-other_star_markers[star][0],other_star_markers[star][1]))

    figsize_x = 6.4
    figsize_y = 4.8
    fig.set_figheight(figsize_y)

def Smooth2DMap(Hist_Old,smooth_size,normalized):

    nbins = Hist_Old.GetNbinsX()
    MapEdge_left = Hist_Old.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Old.GetXaxis().GetBinLowEdge(Hist_Old.GetNbinsX()+1)
    map_size = (MapEdge_right-MapEdge_left)/2.
    Hist_Kernel = ROOT.TH2D("Hist_Kernel","",nbins,-map_size,map_size,nbins,-map_size,map_size)
    Hist_Kernel.Reset()
    for bx1 in range(1,Hist_Old.GetNbinsX()+1):
        for by1 in range(1,Hist_Old.GetNbinsY()+1):
            cell_x = Hist_Kernel.GetXaxis().GetBinCenter(bx1)
            cell_y = Hist_Kernel.GetYaxis().GetBinCenter(by1)
            distance = pow(cell_x*cell_x+cell_y*cell_y,0.5)
            bin_content = ROOT.TMath.Gaus(distance,0,smooth_size)
            Hist_Kernel.SetBinContent(bx1,by1,bin_content)
    #print ('Hist_Kernel.Integral() = %s'%(Hist_Kernel.Integral()))

    Hist_Smooth = Hist_Old.Clone()

    bin_size = Hist_Old.GetXaxis().GetBinCenter(2)-Hist_Old.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(2*smooth_size/bin_size) + 1
    central_bin = int(nbins/2) + 1
    for bx1 in range(1,Hist_Old.GetNbinsX()+1):
        for by1 in range(1,Hist_Old.GetNbinsY()+1):
            old_content = Hist_Old.GetBinContent(bx1,by1)
            old_error = Hist_Old.GetBinError(bx1,by1)
            bin_content = 0
            bin_error = 0
            for bx2 in range(bx1-nbin_smooth,bx1+nbin_smooth):
                for by2 in range(by1-nbin_smooth,by1+nbin_smooth):
                    if bx2<1: 
                        continue
                    if by2<1: 
                        continue
                    if bx2>Hist_Old.GetNbinsX(): 
                        continue
                    if by2>Hist_Old.GetNbinsY(): 
                        continue
                    bin_content += Hist_Kernel.GetBinContent(bx2-bx1+central_bin,by2-by1+central_bin)*Hist_Old.GetBinContent(bx2,by2)
                    bin_error += Hist_Kernel.GetBinContent(bx2-bx1+central_bin,by2-by1+central_bin)*pow(Hist_Old.GetBinError(bx2,by2),2)
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
            Hist_Smooth.SetBinError(bx1,by1,pow(bin_error,0.5))

    if normalized:
        Hist_Smooth.Scale(1./Hist_Kernel.Integral())

    return Hist_Smooth

def GetSignificanceMap(Hist_SR,Hist_Bkg):

    Hist_Skymap = Hist_SR.Clone()
    MapEdge_left = Hist_Skymap.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Skymap.GetXaxis().GetBinLowEdge(Hist_Skymap.GetNbinsX()+1)
    MapEdge_lower = Hist_Skymap.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Skymap.GetYaxis().GetBinLowEdge(Hist_Skymap.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapSize_x = (MapEdge_right-MapEdge_left)/2.
    MapSize_y = (MapEdge_upper-MapEdge_lower)/2.
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_SR.GetBinContent(bx+1,by+1)==0: continue
            if Hist_Bkg.GetBinContent(bx+1,by+1)<=0: continue
            #NSR = Hist_SR.GetBinContent(bx+1,by+1)
            #NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            #NBkgRaw = Hist_Raw_Bkg.GetBinContent(bx+1,by+1)
            #Syst_Err = Hist_Syst.GetBinContent(bx+1,by+1)
            #Syst_Err = 0.
            #alpha = NBkg/NBkgRaw
            #delta_alpha = Syst_Err/NBkg
            #TS = pow((NSR-NBkg)/pow(NSR+pow(Syst_Err*Syst_Err,2),0.5),2)
            #if NSR<15.:
            #    TS = 2.*(NSR*math.log((alpha*(delta_alpha+1)+1)/(alpha*(delta_alpha+1))*(NSR/(NSR+NBkgRaw)))+NBkgRaw*math.log((alpha*(delta_alpha+1)+1)*(NBkgRaw/(NSR+NBkgRaw))))
            #if (NSR-NBkg)>0.:
            #    Hist_Skymap.SetBinContent(bx+1,by+1,pow(TS,0.5))
            #else:
            #    Hist_Skymap.SetBinContent(bx+1,by+1,-1.*pow(TS,0.5))
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            zscore = (NSR-NBkg)/pow(NSR_err*NSR_err,0.5)
            Hist_Skymap.SetBinContent(bx+1,by+1,zscore)
    return Hist_Skymap

def FindCountProjection(Hist_Data_input,proj_type="Y"):

    #n_bins_y = Hist_Data_input.GetNbinsY()
    #n_bins_x = Hist_Data_input.GetNbinsX()
    n_bins_y = min(50,Hist_Data_input.GetNbinsY())
    n_bins_x = min(50,Hist_Data_input.GetNbinsX())
    MapEdge_left = Hist_Data_input.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Data_input.GetXaxis().GetBinLowEdge(Hist_Data_input.GetNbinsX()+1)
    MapEdge_lower = Hist_Data_input.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Data_input.GetYaxis().GetBinLowEdge(Hist_Data_input.GetNbinsY()+1)

    MapCenter_x = (MapEdge_left+MapEdge_right)/2.
    if proj_type=="X":
        MapCenter_x = (MapEdge_upper+MapEdge_lower)/2.

    hist_nbins = n_bins_y
    hist_edge_lower = MapEdge_lower
    hist_edge_upper = MapEdge_upper
    if proj_type=="X":
        hist_nbins = n_bins_x
        hist_edge_lower = MapEdge_left
        hist_edge_upper = MapEdge_right
    Hist_Profile_Y = ROOT.TH1D("Hist_Profile_Y","",hist_nbins,hist_edge_lower,hist_edge_upper)
    for br in range(0,Hist_Profile_Y.GetNbinsX()):
        range_limit = Hist_Profile_Y.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_Y.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        total_error_weight = 0.
        total_cell_size = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                if proj_type=="X":
                    cell_y = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                    cell_x = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                if abs(MapCenter_x-cell_x)>2.0: continue
                if cell_y>=range_limit_previous and cell_y<range_limit:
                    if not data_error==0.:
                        total_cell_size += 1.
                        slice_data += data_content
                        slice_data_err += data_error*data_error
        if total_cell_size==0.: 
            slice_data = 0.
            slice_data_err = 0.
        else:
            slice_data = slice_data
            slice_data_err = pow(slice_data_err,0.5)
        Hist_Profile_Y.SetBinContent(br+1,slice_data)
        Hist_Profile_Y.SetBinError(br+1,slice_data_err)

    profile = []
    profile_err = []
    theta2 = []
    theta2_err = []
    for binx in range(0,Hist_Profile_Y.GetNbinsX()):
        center = Hist_Profile_Y.GetBinCenter(binx+1)
        range_limit = Hist_Profile_Y.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Profile_Y.GetBinLowEdge(binx+1)
        profile_content = Hist_Profile_Y.GetBinContent(binx+1)
        profile_error = Hist_Profile_Y.GetBinError(binx+1)
        theta2 += [center]
        theta2_err += [range_limit-range_limit_previous]
        profile += [profile_content]
        profile_err += [profile_error]

    return profile, profile_err, theta2, theta2_err

def BackgroundSubtractMap(fig,hist_data,hist_bkgd,label_x,label_y,label_z,plotname):

    colormap = 'coolwarm'
    #colormap = 'viridis'

    Old_map_nbins_x = hist_data.GetNbinsX()
    Old_map_nbins_y = hist_data.GetNbinsY()
    Old_MapEdge_left = hist_data.GetXaxis().GetBinLowEdge(1)
    Old_MapEdge_right = hist_data.GetXaxis().GetBinLowEdge(hist_data.GetNbinsX()+1)
    Old_MapEdge_lower = hist_data.GetYaxis().GetBinLowEdge(1)
    Old_MapEdge_upper = hist_data.GetYaxis().GetBinLowEdge(hist_data.GetNbinsY()+1)
    Old_MapEdge_center_x = 0.5*(Old_MapEdge_left+Old_MapEdge_right)
    Old_MapEdge_center_y = 0.5*(Old_MapEdge_lower+Old_MapEdge_upper)
    Old_MapEdge_size_x = Old_MapEdge_right-Old_MapEdge_center_x
    Old_MapEdge_size_y = Old_MapEdge_upper-Old_MapEdge_center_y

    map_bin_size = 2.*Old_MapEdge_size_x/float(Old_map_nbins_x)
    map_nbins_x = Old_map_nbins_x
    map_nbins_y = Old_map_nbins_y
    MapEdge_left = Old_MapEdge_center_x-int(map_nbins_x/2)*map_bin_size
    MapEdge_right = Old_MapEdge_center_x+int(map_nbins_x/2)*map_bin_size
    MapEdge_lower = Old_MapEdge_center_y-int(map_nbins_y/2)*map_bin_size
    MapEdge_upper = Old_MapEdge_center_y+int(map_nbins_y/2)*map_bin_size

    deg_per_bin = (MapEdge_right-MapEdge_left)/map_nbins_x
    nbins_per_deg = map_nbins_x/(MapEdge_right-MapEdge_left)
    old_x_axis = np.linspace(Old_MapEdge_left,Old_MapEdge_right,Old_map_nbins_x)
    old_y_axis = np.linspace(Old_MapEdge_lower,Old_MapEdge_upper,Old_map_nbins_y)
    x_axis = np.linspace(MapEdge_left,MapEdge_right,map_nbins_x)
    prelim_x_axis_sparse = np.linspace(MapEdge_left,MapEdge_right,6)
    x_axis_sparse = []
    for i in prelim_x_axis_sparse:
        if int(i)>MapEdge_left and int(i)<MapEdge_right:
            x_axis_sparse += [int(i)]
    x_axis_reflect = ["{:6.1f}".format(-1.*i) for i in x_axis_sparse]
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,map_nbins_y)

    grid_data = np.zeros((map_nbins_y, map_nbins_x))
    grid_bkgd = np.zeros((map_nbins_y, map_nbins_x))
    grid_excs = np.zeros((map_nbins_y, map_nbins_x))
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            #hist_bin_x = hist_data.GetXaxis().FindBin(x_axis[xbin])
            #hist_bin_y = hist_data.GetYaxis().FindBin(y_axis[ybin])
            hist_bin_x = xbin+1
            hist_bin_y = ybin+1
            if hist_bin_x<1: continue
            if hist_bin_y<1: continue
            if hist_bin_x>hist_data.GetNbinsX(): continue
            if hist_bin_y>hist_data.GetNbinsY(): continue
            grid_data[ybin,xbin] = hist_data.GetBinContent(hist_bin_x,hist_bin_y)
            grid_bkgd[ybin,xbin] = hist_bkgd.GetBinContent(hist_bin_x,hist_bin_y)
            grid_excs[ybin,xbin] = hist_data.GetBinContent(hist_bin_x,hist_bin_y)-hist_bkgd.GetBinContent(hist_bin_x,hist_bin_y)

    # Define the locations for the axes
    left, width = 0.12, 0.6
    bottom, height = 0.12, 0.6
    bottom_h = bottom+height+0.03
    left_h = left+width+0.03
     
    # Set up the geometry of the three plots
    rect_temperature = [left, bottom, width, height] # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.20] # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.20, height] # dimensions of y-histogram

    # Set up the size of the figure
    #fig = plt.figure(1, figsize=(9.5,9))
    fig = plt.figure(1)
    fig.set_figheight(8)
    fig.set_figwidth(8)

    fig.clf()
    # Make the three plots
    axTemperature = plt.axes(rect_temperature) # temperature plot
    axHistx = plt.axes(rect_histx) # x histogram
    axHisty = plt.axes(rect_histy) # y histogram

    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Plot the temperature data
    cax = (axTemperature.imshow(grid_data, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()), origin='lower', cmap=colormap))

    #Plot the axes labels
    axTemperature.set_xlabel(label_x)
    axTemperature.set_ylabel(label_y)
    axTemperature.set_xticks(x_axis_sparse)
    axTemperature.set_xticklabels(x_axis_reflect)

    x_profile, x_profile_stat_err, x_ax, x_ax_err = FindCountProjection(hist_data,proj_type='X')
    y_profile, y_profile_stat_err, y_ax, y_ax_err = FindCountProjection(hist_data,proj_type='Y')
    x_bkgd_profile, x_bkgd_profile_stat_err, x_ax, x_ax_err = FindCountProjection(hist_bkgd,proj_type='X')
    y_bkgd_profile, y_bkgd_profile_stat_err, y_ax, y_ax_err = FindCountProjection(hist_bkgd,proj_type='Y')
    #Plot the histograms
    axHistx.errorbar(x_ax,x_profile,yerr=x_profile_stat_err,color='k',marker='.',ls='none')
    axHisty.errorbar(y_profile,y_ax,xerr=y_profile_stat_err,color='k',marker='.',ls='none')
    axHistx.plot(x_ax,x_bkgd_profile,color='r',ls='solid')
    axHisty.plot(y_bkgd_profile,y_ax,color='r',ls='solid')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1,1))
    axHistx.yaxis.set_major_formatter(formatter)
    axHisty.xaxis.set_major_formatter(formatter)

    fig.savefig("output_plots/Before_%s.png"%(plotname),bbox_inches='tight')

    fig.clf()
    # Make the three plots
    axTemperature = plt.axes(rect_temperature) # temperature plot
    axHistx = plt.axes(rect_histx) # x histogram
    axHisty = plt.axes(rect_histy) # y histogram

    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Plot the temperature data
    cax = (axTemperature.imshow(grid_excs, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()), origin='lower', cmap=colormap))

    #Plot the axes labels
    axTemperature.set_xlabel(label_x)
    axTemperature.set_ylabel(label_y)
    axTemperature.set_xticks(x_axis_sparse)
    axTemperature.set_xticklabels(x_axis_reflect)

    nbins = hist_data.GetNbinsX()
    MapEdge_left = hist_data.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_data.GetXaxis().GetBinLowEdge(hist_data.GetNbinsX()+1)
    MapEdge_bottom = hist_data.GetYaxis().GetBinLowEdge(1)
    MapEdge_top = hist_data.GetYaxis().GetBinLowEdge(hist_data.GetNbinsY()+1)
    hist_excs = ROOT.TH2D("hist_excs","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_bottom,MapEdge_top)
    hist_excs.Reset()
    hist_excs.Add(hist_data)
    hist_excs.Add(hist_bkgd,-1.)

    x_profile, x_profile_stat_err, x_ax, x_ax_err = FindCountProjection(hist_excs,proj_type='X')
    y_profile, y_profile_stat_err, y_ax, y_ax_err = FindCountProjection(hist_excs,proj_type='Y')
    #Plot the histograms
    axHistx.errorbar(x_ax,x_profile,yerr=x_profile_stat_err,color='k',marker='.',ls='none')
    axHisty.errorbar(y_profile,y_ax,xerr=y_profile_stat_err,color='k',marker='.',ls='none')

    fig.savefig("output_plots/After_%s.png"%(plotname),bbox_inches='tight')

