
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

