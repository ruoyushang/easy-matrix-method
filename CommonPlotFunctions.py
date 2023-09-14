
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

input_path = '/gamma_raid/userspace/rshang/SMI_output'

#folder_path = 'output_j1908_paper'
#folder_path = 'output_nuclear_v487'
#folder_tag = 'paper'
#energy_bin = [100.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.]

folder_path = 'output_test_2'
folder_tag = '_test'
#folder_path = 'output_nuclear_3tel'
#folder_tag = '_3tel'

#folder_path = 'output_weight_log_m0p0'
#folder_path = 'output_weight_log_m0p5'
#folder_path = 'output_weight_log_m1p0'
#folder_path = 'output_weight_log_m1p5'
#folder_path = 'output_weight_log_m2p0'

#analysis_method = 'FoV'
#analysis_method = 'Ratio'
#analysis_method = 'Regression'
#analysis_method = 'Init_Perturbation'
analysis_method = 'Perturbation'
#analysis_method = 'Combined'

energy_bin = [100.,159.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.]
#energy_bin = [100.,167.,300.,538.,965.,1732.,3107.,5574.,10000.]

smooth_size_spectroscopy = 0.07
calibration_radius = 0.15 # need to be larger than the PSF and smaller than the integration radius
#calibration_radius = 0.2 # need to be larger than the PSF and smaller than the integration radius
#calibration_radius = 0.3 # need to be larger than the PSF and smaller than the integration radius
n_xoff_bins = 1
n_yoff_bins = 1

MSCW_lower_blind = -0.5
MSCL_lower_blind = -0.7
MSCW_upper_blind = 0.7
MSCL_upper_blind = 0.5
n_extra_lower_bins = 0
n_extra_upper_bins = 4
mtx_dim_w_fine = 4
mtx_dim_l_fine = 4

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

def ConvertRaDecToGalactic(ra, dec):
    delta = dec*ROOT.TMath.Pi()/180.
    delta_G = 27.12825*ROOT.TMath.Pi()/180.
    alpha = ra*ROOT.TMath.Pi()/180.
    alpha_G = 192.85948*ROOT.TMath.Pi()/180.
    l_NCP = 122.93192*ROOT.TMath.Pi()/180.
    sin_b = ROOT.TMath.Sin(delta)*ROOT.TMath.Sin(delta_G)+ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(delta_G)*ROOT.TMath.Cos(alpha-alpha_G)
    cos_b = ROOT.TMath.Cos(ROOT.TMath.ASin(sin_b))
    sin_l_NCP_m_l = ROOT.TMath.Cos(delta)*ROOT.TMath.Sin(alpha-alpha_G)/cos_b
    cos_l_NCP_m_l = (ROOT.TMath.Cos(delta_G)*ROOT.TMath.Sin(delta)-ROOT.TMath.Sin(delta_G)*ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(alpha-alpha_G))/cos_b
    b = (ROOT.TMath.ASin(sin_b))*180./ROOT.TMath.Pi()
    l = (l_NCP-ROOT.TMath.ATan2(sin_l_NCP_m_l,cos_l_NCP_m_l))*180./ROOT.TMath.Pi()
    return l, b

def ReadBrightStarListFromFile():
    star_name = []
    star_ra = []
    star_dec = []
    inputFile = open('Hipparcos_MAG8_1997.dat')
    for line in inputFile:
        if line[0]=="#": continue
        if line[0]=="*": continue
        if len(line.split())<5: continue
        target_ra = float(line.split()[0])
        target_dec = float(line.split()[1])
        target_brightness = float(line.split()[3])+float(line.split()[4])
        if target_brightness>7.: continue
        #if target_brightness>6.: continue
        #if target_brightness>5.: continue
        star_ra += [target_ra]
        star_dec += [target_dec]
        star_name += ['bmag = %0.1f'%(target_brightness)]
    #print ('Found %s bright stars.'%(len(star_name)))
    return star_name, star_ra, star_dec

def HMS2deg(ra='', dec=''):
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split(':')]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
    if ra:
        H, M, S = [float(i) for i in ra.split(':')]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)           
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC

def ReadSNRTargetListFromCSVFile():
    source_name = []
    source_ra = []
    source_dec = []
    source_size = []
    with open('SNRcat20221001-SNR.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            if len(row)==0: continue
            if '#' in row[0]: continue
            target_name = row[0]
            target_min_dist = row[13]
            if target_min_dist=='':
                target_min_dist = '1000'
            #if float(target_min_dist)>6.: continue
            target_size = row[15]
            if target_size=='':
                target_size = 0.
            target_ra = row[19]
            target_dec = row[20]
            source_name += [target_name]
            source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
            source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
            source_size += [0.5*float(target_size)/60.]
            #print('target_min_dist = %s'%(target_min_dist))
            #print('source_name = %s'%(source_name[len(source_name)-1]))
            #print('source_ra = %0.2f'%(source_ra[len(source_ra)-1]))
            #print('source_dec = %0.2f'%(source_dec[len(source_dec)-1]))
            #print(row)
    return source_name, source_ra, source_dec, source_size

def ReadATNFTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_age = []
    source_edot = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        target_name = line.split(',')[0].strip(" ")
        if target_name=="\n": continue
        target_ra = line.split(',')[1].strip(" ")
        target_dec = line.split(',')[2].strip(" ")
        target_dist = line.split(',')[3].strip(" ")
        target_age = line.split(',')[4].strip(" ")
        target_edot = line.split(',')[5].strip(" ")
        if target_dist=='*': continue
        if target_age=='*': continue
        if target_edot=='*': continue
        target_brightness = float(target_edot)/pow(float(target_dist),2)

        if float(target_edot)<1e35: continue
        #if float(target_dist)>6.: continue

        #ra_deg = float(HMS2deg(target_ra,target_dec)[0])
        #dec_deg = float(HMS2deg(target_ra,target_dec)[1])
        #gal_l, gal_b = ConvertRaDecToGalactic(ra_deg,dec_deg)
        #if abs(gal_b)<5.: continue

        source_name += [target_name]
        source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
        source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
        source_dist += [float(target_dist)]
        source_age += [float(target_age)]
        source_edot += [float(target_edot)]
    return source_name, source_ra, source_dec, source_dist, source_age

def ReadLhaasoListFromFile():
    source_name = []
    source_ra = []
    source_dec = []
    file_path = 'LHAASO_1st_Catalog.txt'
    inputFile = open(file_path)
    for line in inputFile:
        source_name += ['%s %s'%(line.split()[0],line.split()[1])]
        source_ra += [float(line.split()[3])]
        source_dec += [float(line.split()[4])]
    return source_name, source_ra, source_dec

def ReadFermiCatelog():
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open('gll_psc_v26.xml')
    target_name = ''
    target_type = ''
    target_info = ''
    target_flux = ''
    target_ra = ''
    target_dec = ''
    for line in inputFile:
        if line.split(' ')[0]=='<source':
            for block in range(0,len(line.split(' '))):
                if 'Unc_' in line.split(' ')[block]: continue
                if 'name=' in line.split(' ')[block]:
                    target_name = line.split('name="')[1].split('"')[0]
                if 'type=' in line.split(' ')[block]:
                    target_type = line.split('type="')[1].split('"')[0]
                if 'Flux1000=' in line.split(' ')[block]:
                    target_flux = line.split('Flux1000="')[1].split('"')[0]
                if 'Energy_Flux100=' in line.split(' ')[block]:
                    target_info = line.split(' ')[block]
                    target_info = target_info.strip('>\n')
                    target_info = target_info.strip('"')
                    target_info = target_info.lstrip('Energy_Flux100="')
        if '<parameter' in line and 'name="RA"' in line:
            for block in range(0,len(line.split(' '))):
                if 'value=' in line.split(' ')[block]:
                    target_ra = line.split(' ')[block].split('"')[1]
        if '<parameter' in line and 'name="DEC"' in line:
            for block in range(0,len(line.split(' '))):
                if 'value=' in line.split(' ')[block]:
                    target_dec = line.split(' ')[block].split('"')[1]
        if 'source>' in line:
            if target_ra=='': 
                target_name = ''
                target_type = ''
                target_info = ''
                target_ra = ''
                target_dec = ''
                continue
            #if target_type=='PointSource': 
            #if float(target_info)<1e-5: 
            if float(target_flux)<1e-9: 
                target_name = ''
                target_type = ''
                target_info = ''
                target_ra = ''
                target_dec = ''
                continue
            #print ('target_name = %s'%(target_name))
            #print ('target_type = %s'%(target_type))
            #print ('target_ra = %s'%(target_ra))
            #print ('target_dec = %s'%(target_dec))
            source_name += [target_name]
            #source_name += ['%0.2e'%(float(target_info))]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            target_name = ''
            target_type = ''
            target_ra = ''
            target_dec = ''
    return source_name, source_ra, source_dec

def ReadHAWCTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        if '- name:' in line:
            target_name = line.lstrip('   - name: ')
            target_name = target_name.strip('\n')
        if 'RA:' in line:
            target_ra = line.lstrip('     RA: ')
        if 'Dec:' in line:
            target_dec = line.lstrip('     Dec: ')
        if 'flux measurements:' in line:
            source_name += [target_name]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            target_name = ''
            target_ra = ''
            target_dec = ''
    return source_name, source_ra, source_dec

def GetGammaSourceInfo():

    other_stars = []
    other_stars_type = []
    other_star_coord = []

    #return other_stars, other_stars_type, other_star_coord

    near_source_cut = 0.1

    drawBrightStar = False
    drawPulsar = True
    drawSNR = False
    drawLHAASO = False
    drawFermi = False
    drawHAWC = False
    drawTeV = False

    if drawBrightStar:
        star_name, star_ra, star_dec = ReadBrightStarListFromFile()
        for src in range(0,len(star_name)):
            src_ra = star_ra[src]
            src_dec = star_dec[src]
            other_stars += [star_name[src]]
            other_stars_type += ['Star']
            other_star_coord += [[src_ra,src_dec,0.]]

    if drawPulsar:
        target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age = ReadATNFTargetListFromFile('ATNF_pulsar_full_list.txt')
        for src in range(0,len(target_psr_name)):
            gamma_source_name = target_psr_name[src]
            gamma_source_ra = target_psr_ra[src]
            gamma_source_dec = target_psr_dec[src]
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

    if drawSNR:
        target_snr_name, target_snr_ra, target_snr_dec, target_snr_size = ReadSNRTargetListFromCSVFile()
        for src in range(0,len(target_snr_name)):
            gamma_source_name = target_snr_name[src]
            gamma_source_ra = target_snr_ra[src]
            gamma_source_dec = target_snr_dec[src]
            gamma_source_size = target_snr_size[src]
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['SNR']
                other_star_coord += [[gamma_source_ra,gamma_source_dec,gamma_source_size]]

    if drawHAWC:
        target_hwc_name, target_hwc_ra, target_hwc_dec = ReadHAWCTargetListFromFile('Cat_3HWC.txt')
        for src in range(0,len(target_hwc_name)):
            gamma_source_name = target_hwc_name[src]
            gamma_source_ra = target_hwc_ra[src]
            gamma_source_dec = target_hwc_dec[src]
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['HAWC']
                other_star_coord += [[gamma_source_ra,gamma_source_dec,0.]]

    if drawLHAASO:
        target_lhs_name, target_lhs_ra, target_lhs_dec = ReadLhaasoListFromFile()
        for src in range(0,len(target_lhs_name)):
            gamma_source_name = target_lhs_name[src]
            gamma_source_ra = target_lhs_ra[src]
            gamma_source_dec = target_lhs_dec[src]
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['LHAASO']
                other_star_coord += [[gamma_source_ra,gamma_source_dec,0.]]

    if drawFermi:
        fermi_name, fermi_ra, fermi_dec = ReadFermiCatelog()
        for src in range(0,len(fermi_name)):
            gamma_source_name = fermi_name[src]
            gamma_source_ra = fermi_ra[src]
            gamma_source_dec = fermi_dec[src]
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['Fermi']
                other_star_coord += [[gamma_source_ra,gamma_source_dec,0.]]

    if drawTeV:
        inputFile = open('TeVCat_RaDec_w_Names.txt')
        for line in inputFile:
            gamma_source_name = line.split(',')[0]
            gamma_source_ra = float(line.split(',')[1])
            gamma_source_dec = float(line.split(',')[2])
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


def MatplotlibMap2D(hist_map,hist_tone,hist_contour,fig,label_x,label_y,label_z,plotname,roi_x=[],roi_y=[],roi_r=[],rotation_angle=0.,colormap='coolwarm',psf=0.):

    print ('Making plot %s...'%(plotname))

    #colormap = 'coolwarm'
    #if 'SkymapCOMap' in plotname:
    #    colormap = 'gray'

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

        #skymap_mean = hist_tone.GetMean()
        #skymap_rms = hist_tone.GetRMS()
        #print ('skymap_mean = %s'%(skymap_mean))
        #print ('skymap_rms = %s'%(skymap_rms))
        #min_z = skymap_mean-2.*skymap_rms
        #max_z = skymap_mean+2.*skymap_rms

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
            max_z_contour = 1.0*hist_contour[ctr].GetMaximum()
            min_z_contour = 0.5*hist_contour[ctr].GetMaximum()
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
                    if label_z!='significance' and max_z_contour>0.:
                        grid_contour[ybin,xbin] = hist_contour[ctr].GetBinContent(hist_bin_x,hist_bin_y)*max_z/max_z_contour
            list_grid_contour += [grid_contour]

    other_stars, other_star_type, other_star_coord = GetGammaSourceInfo() 

    other_star_labels = []
    other_star_types = []
    other_star_markers = []
    star_range = 0.8*(MapEdge_right-MapEdge_left)/2.
    source_ra = (MapEdge_left+MapEdge_right)/2.
    source_dec = (MapEdge_lower+MapEdge_upper)/2.
    n_stars = 0
    for star in range(0,len(other_stars)):
        if abs(source_ra+other_star_coord[star][0])>star_range: continue
        if abs(source_dec-other_star_coord[star][1])>star_range: continue
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
    if label_z=='significance' and not 'SkymapHAWC' in plotname:
        max_z = 5.
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=-max_z,vmax=max_z,zorder=0)
    elif 'SkymapHAWC' in plotname:
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=0,zorder=0)
    else:
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=min_z,vmax=max_z,zorder=0)

    list_colors = ['tomato','lime','deepskyblue']
    list_styles = ['solid','solid','solid']
    for ctr in range(0,len(list_grid_contour)):
        axbig.contour(list_grid_contour[len(list_grid_contour)-1-ctr], list_levels[len(list_grid_contour)-1-ctr], linestyles=list_styles[len(list_grid_contour)-1-ctr], colors=list_colors[len(list_grid_contour)-1-ctr], extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=1)
    favorite_color = 'k'
    if colormap=='gray':
        favorite_color = 'r'
    for star in range(0,len(other_star_markers)):
        marker_size = 60
        if other_star_types[star]=='PSR':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=1.5*marker_size, c=favorite_color, marker='+', label=other_star_labels[star])
        if other_star_types[star]=='SNR':
            #axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c=favorite_color, marker='^', label=other_star_labels[star])
            mycircle = plt.Circle( (other_star_markers[star][0], other_star_markers[star][1]), other_star_markers[star][2], fill = False, color=favorite_color)
            axbig.add_patch(mycircle)
        if other_star_types[star]=='LHAASO':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c=favorite_color, marker='^', label=other_star_labels[star])
        if other_star_types[star]=='HAWC':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c=favorite_color, marker='o', label=other_star_labels[star])
        if other_star_types[star]=='Fermi':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c=favorite_color, marker='s', label=other_star_labels[star])
        if other_star_types[star]=='MSP':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='tomato', marker='+', label=other_star_labels[star])
        if other_star_types[star]=='TeV':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='lime', marker='+', label=other_star_labels[star])
        if other_star_types[star]=='Star':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='k', marker='o', label=other_star_labels[star])
        text_offset_x = 0.07
        text_offset_y = 0.07

    if psf>0.:
        circle_center_x = MapCenter_x-0.8*(MapEdge_right-MapEdge_left)/2.
        circle_center_y = MapCenter_y-0.8*(MapEdge_upper-MapEdge_lower)/2.
        if psf>0.4:
            circle_center_x = MapCenter_x-0.6*(MapEdge_right-MapEdge_left)/2.
            circle_center_y = MapCenter_y-0.6*(MapEdge_upper-MapEdge_lower)/2.
        mycircle = plt.Circle((circle_center_x, circle_center_y), psf, color='white', fill=False)
        axbig.add_patch(mycircle)

    divider = make_axes_locatable(axbig)
    cax = divider.append_axes("bottom", size="5%", pad=0.7)
    cbar = fig.colorbar(im,orientation="horizontal",cax=cax)
    cbar.set_label(label_z)
    axbig.set_xticks(x_axis_sparse)
    axbig.set_xticklabels(x_axis_reflect)

    for roi in range(0,len(roi_r)):
        mycircle = plt.Circle((-roi_x[roi], roi_y[roi]), roi_r[roi], color='w', linestyle='dashed', fill=False)
        axbig.add_patch(mycircle)

    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')

    #if len(roi_r)>0:
    #for roi in range(0,len(roi_r)):
    #    mycircle = plt.Circle((-roi_x[roi], roi_y[roi]), roi_r[roi], color='w', linestyle='dashed', fill=False)
    #    axbig.add_patch(mycircle)
    #axbig.legend(bbox_to_anchor=(0.7, 1), borderaxespad=0, fontsize=7)
    #fig.savefig("output_plots/%s_legend.png"%(plotname),bbox_inches='tight')

    axbig.remove()

    if 'SkymapSignificance_Sum' in plotname:
        for star in range(0,len(other_star_markers)):
            print ('%s, (%0.2f,%0.2f)'%(other_star_labels[star],-other_star_markers[star][0],other_star_markers[star][1]))

    figsize_x = 6.4
    figsize_y = 4.8
    fig.set_figheight(figsize_y)

def Smooth2DMap(Hist_Old,smooth_size,normalized):

    #return Hist_Old

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
            NSR_err = max(1.,NSR_err)
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

def FindExtension(Hist_Data_input,roi_x,roi_y,integration_range):

    global calibration_radius

    n_bins_2d = Hist_Data_input.GetNbinsX()
    n_bins_1d = int(integration_range/0.1)

    n_bins_y = Hist_Data_input.GetNbinsY()
    n_bins_x = Hist_Data_input.GetNbinsX()
    MapEdge_left = Hist_Data_input.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Data_input.GetXaxis().GetBinLowEdge(Hist_Data_input.GetNbinsX()+1)
    MapEdge_lower = Hist_Data_input.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Data_input.GetYaxis().GetBinLowEdge(Hist_Data_input.GetNbinsY()+1)
    cell_size = (MapEdge_right-MapEdge_left)/float(n_bins_x)*(MapEdge_upper-MapEdge_lower)/float(n_bins_y)

    Hist_Profile_Theta2 = ROOT.TH1D("Hist_Profile_Theta2","",n_bins_1d,0,integration_range)
    for br in range(0,Hist_Profile_Theta2.GetNbinsX()):
        range_limit = Hist_Profile_Theta2.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_Theta2.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        total_error_weight = 0.
        total_cell_size = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                distance_sq = pow(cell_x-roi_x,2)+pow(cell_y-roi_y,2)
                #data_content = Hist_Data_input.GetBinContent(bx+1,by+1)/cell_size
                #data_error = Hist_Data_input.GetBinError(bx+1,by+1)/cell_size
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                if distance_sq>=pow(range_limit_previous,2) and distance_sq<pow(range_limit,2):
                    if not data_error==0.:
                        slice_data += data_content
                        slice_data_err += data_error*data_error
                        total_cell_size += cell_size
        if total_cell_size==0.: 
            slice_data = 0.
            slice_data_err = 0.
        else:
            slice_data = slice_data/total_cell_size
            slice_data_err = pow(slice_data_err,0.5)/total_cell_size
        Hist_Profile_Theta2.SetBinContent(br+1,slice_data)
        Hist_Profile_Theta2.SetBinError(br+1,slice_data_err)

    profile = []
    profile_err = []
    theta2 = []
    theta2_err = []
    for binx in range(0,Hist_Profile_Theta2.GetNbinsX()):
        center = Hist_Profile_Theta2.GetBinCenter(binx+1)
        range_limit = Hist_Profile_Theta2.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Profile_Theta2.GetBinLowEdge(binx+1)
        profile_content = Hist_Profile_Theta2.GetBinContent(binx+1)
        profile_error = Hist_Profile_Theta2.GetBinError(binx+1)
        if center>integration_range: continue
        theta2 += [center]
        theta2_err += [range_limit-range_limit_previous]
        profile += [profile_content]
        profile_err += [profile_error]

    return profile, profile_err, theta2, theta2_err

def GetRegionIntegral(hist_data_skymap,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r):

    flux_sum = 0.
    flux_stat_err = 0.
    for bx in range(0,hist_data_skymap.GetNbinsX()):
        for by in range(0,hist_data_skymap.GetNbinsY()):
            bin_ra = hist_data_skymap.GetXaxis().GetBinCenter(bx+1)
            bin_dec = hist_data_skymap.GetYaxis().GetBinCenter(by+1)
            keep_event = False
            for roi in range(0,len(roi_x)):
                distance = pow(pow(bin_ra-roi_x[roi],2) + pow(bin_dec-roi_y[roi],2),0.5)
                if distance<roi_r[roi]: 
                    keep_event = True
            for roi in range(0,len(excl_roi_x)):
                excl_distance = pow(pow(bin_ra-excl_roi_x[roi],2) + pow(bin_dec-excl_roi_y[roi],2),0.5)
                if excl_distance<excl_roi_r[roi]: 
                    keep_event = True
            if keep_event:
                flux_sum += hist_data_skymap.GetBinContent(bx+1,by+1)
                flux_stat_err += pow(hist_data_skymap.GetBinError(bx+1,by+1),2)
    flux_stat_err = pow(flux_stat_err,0.5)
    return flux_sum, flux_stat_err

def GetRegionSpectrum(hist_data_skymap,ebin_low,ebin_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r):

    x_axis = []
    x_error = []
    y_axis = []
    y_error = []
    for ebin in range(ebin_low,ebin_up):
        flux_sum = 0.
        flux_stat_err = 0.
        flux_syst_err = 0.
        flux_sum, flux_stat_err = GetRegionIntegral(hist_data_skymap[ebin],roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
        x_axis += [0.5*(energy_bin[ebin]+energy_bin[ebin+1])]
        x_error += [0.5*(energy_bin[ebin+1]-energy_bin[ebin])]
        y_axis += [flux_sum]
        y_error += [flux_stat_err]

    return x_axis, x_error, y_axis, y_error

def GetFITSMap(map_file, hist_map, isRaDec):

    hist_map.Reset()
    nbins_x = hist_map.GetNbinsX()
    nbins_y = hist_map.GetNbinsY()
    MapEdge_left = hist_map.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_map.GetXaxis().GetBinLowEdge(hist_map.GetNbinsX()+1)
    MapEdge_lower = hist_map.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_map.GetYaxis().GetBinLowEdge(hist_map.GetNbinsY()+1)
    MapCenter_x = (MapEdge_left+MapEdge_right)/2.
    MapCenter_y = (MapEdge_lower+MapEdge_upper)/2.

    hdu = fits.open(map_file)[0]
    wcs = WCS(hdu.header)
    image_data = hdu.data

    print ("wcs")
    print (wcs)

    image_shape = np.shape(image_data)
    for binx in range(0,nbins_x):
        for biny in range(0,nbins_y):
            ra = hist_map.GetXaxis().GetBinCenter(binx+1)
            dec = hist_map.GetYaxis().GetBinCenter(biny+1)
            pixs = wcs.all_world2pix(ra,dec,1)
            if int(pixs[0])<0: continue
            if int(pixs[1])<0: continue
            if int(pixs[0])>=image_shape[0]: continue
            if int(pixs[1])>=image_shape[1]: continue
            fits_data = image_data[int(pixs[1]),int(pixs[0])]
            hist_map.SetBinContent(binx+1,biny+1,fits_data)

    return hist_map

def GetSlicedGalfaHIDataCubeMap(map_file, hist_map, vel_low, vel_up, reset_hist):

    if reset_hist:
        hist_map.Reset()
    nbinsx = hist_map.GetNbinsX()
    nbinsy = hist_map.GetNbinsY()
    map_center_ra = hist_map.GetXaxis().GetBinCenter(int(float(nbinsx)/2.))
    map_center_dec = hist_map.GetYaxis().GetBinCenter(int(float(nbinsy)/2.))

    filename = map_file

    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)
    image_data = hdu.data
    print ("wcs")
    print (wcs)

    pixs_start = wcs.all_world2pix(map_center_ra,map_center_dec,vel_low,1)
    pixs_end = wcs.all_world2pix(map_center_ra,map_center_dec,vel_up,1)
    vel_idx_start = int(pixs_start[2])
    vel_idx_end = int(pixs_end[2])

    image_data_reduced_z = np.full((image_data[vel_idx_start, :, :].shape),0.)
    for idx in range(vel_idx_start,vel_idx_end):
        world_coord = wcs.all_pix2world(0,0,idx,1) 
        velocity = world_coord[2]
        world_coord = wcs.all_pix2world(0,0,idx+1,1) 
        velocity_next = world_coord[2]
        delta_vel = velocity_next - velocity
        image_data_reduced_z += image_data[idx, :, :]*delta_vel

    for binx in range(0,nbinsx):
        for biny in range(0,nbinsy):
            map_ra = hist_map.GetXaxis().GetBinCenter(binx+1)
            map_dec = hist_map.GetYaxis().GetBinCenter(biny+1)
            map_pixs = wcs.all_world2pix(map_ra, map_dec, vel_low, 1)
            pix_ra = int(map_pixs[0])
            pix_dec = int(map_pixs[1])
            if pix_ra<0: continue
            if pix_dec<0: continue
            if pix_ra>=image_data_reduced_z[:,:].shape[1]: continue
            if pix_dec>=image_data_reduced_z[:,:].shape[0]: continue
            hist_map.SetBinContent(binx+1,biny+1,image_data_reduced_z[pix_dec,pix_ra])

def GetSlicedDataCubeMap(map_file, hist_map, vel_low, vel_up):

    hist_map.Reset()
    nbinsx = hist_map.GetNbinsX()
    nbinsy = hist_map.GetNbinsY()
    map_center_ra = hist_map.GetXaxis().GetBinCenter(int(float(nbinsx)/2.))
    map_center_dec = hist_map.GetYaxis().GetBinCenter(int(float(nbinsy)/2.))
    map_center_lon, map_center_lat = ConvertRaDecToGalactic(map_center_ra, map_center_dec)

    filename = map_file

    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)
    image_data = hdu.data
    print ("wcs")
    print (wcs)

    pixs_start = wcs.all_world2pix(vel_low,map_center_lon,map_center_lat,1)
    pixs_end = wcs.all_world2pix(vel_up,map_center_lon,map_center_lat,1)
    vel_idx_start = int(pixs_start[0])
    vel_idx_end = int(pixs_end[0])

    image_data_reduced_z = np.full((image_data[:, :, vel_idx_start].shape),0.)
    for idx in range(vel_idx_start,vel_idx_end):
        world_coord = wcs.all_pix2world(idx,0,0,1) 
        velocity = world_coord[0]
        world_coord = wcs.all_pix2world(idx+1,0,0,1) 
        velocity_next = world_coord[0]
        delta_vel = velocity_next - velocity
        image_data_reduced_z += image_data[:, :, idx]*delta_vel

    for binx in range(0,nbinsx):
        for biny in range(0,nbinsy):
            map_ra = hist_map.GetXaxis().GetBinCenter(binx+1)
            map_dec = hist_map.GetYaxis().GetBinCenter(biny+1)
            map_lon, map_lat = ConvertRaDecToGalactic(map_ra, map_dec)
            map_pixs = wcs.all_world2pix(vel_low, map_lon, map_lat, 1)
            pix_lon = int(map_pixs[1])
            pix_lat = int(map_pixs[2])
            hist_map.SetBinContent(binx+1,biny+1,image_data_reduced_z[pix_lat,pix_lon])

def GetHealpixMap(map_file, hist_map, isRaDec):

    hist_map.Reset()
    nbins_x = hist_map.GetNbinsX()
    nbins_y = hist_map.GetNbinsY()
    MapEdge_left = hist_map.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_map.GetXaxis().GetBinLowEdge(hist_map.GetNbinsX()+1)
    MapEdge_lower = hist_map.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_map.GetYaxis().GetBinLowEdge(hist_map.GetNbinsY()+1)
    MapCenter_x = (MapEdge_left+MapEdge_right)/2.
    MapCenter_y = (MapEdge_lower+MapEdge_upper)/2.

    hpx, header = hp.read_map(map_file, field=0, h=True)
    #hpx, header = hp.read_map(map_file, field=1, h=True)
    npix = len(hpx)
    nside = hp.npix2nside(npix)
    for ipix in range(0,npix):
        theta, phi = hp.pix2ang(nside, ipix)
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5 * np.pi - theta)
        fits_data = hpx[ipix]
        binx = hist_map.GetXaxis().FindBin(ra)
        biny = hist_map.GetYaxis().FindBin(dec)
        if binx<1: continue
        if biny<1: continue
        if binx>hist_map.GetNbinsX(): continue
        if biny>hist_map.GetNbinsY(): continue
        hist_map.SetBinContent(binx,biny,fits_data)

    return hist_map

def MatplotlibHist2D(hist_map,fig,label_x,label_y,label_z,plotname,zmax=0,zmin=0):

    top = cm.get_cmap('Blues_r', 128) # r means reversed version
    bottom = cm.get_cmap('Oranges', 128)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')
    colormap = 'coolwarm'
    #colormap = 'viridis'

    map_nbins = hist_map.GetNbinsX()
    MapEdge_left = hist_map.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_map.GetXaxis().GetBinLowEdge(hist_map.GetNbinsX()+1)
    MapEdge_lower = hist_map.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_map.GetYaxis().GetBinLowEdge(hist_map.GetNbinsY()+1)

    deg_per_bin = (MapEdge_right-MapEdge_left)/map_nbins
    nbins_per_deg = map_nbins/(MapEdge_right-MapEdge_left)
    x_axis = np.linspace(MapEdge_left,MapEdge_right,map_nbins)
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,map_nbins)
    grid_z = np.zeros((map_nbins, map_nbins))
    max_z = 0.
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

    fig.clf()
    figsize_x = 6.
    figsize_y = 6.
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    if zmax==0 and zmin==0:
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=0)
    else:
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=0,vmin=zmin,vmax=zmax)

    #MSCW_lower_blind = -0.6
    #MSCL_lower_blind = -0.7
    #MSCW_upper_blind = 0.6
    #MSCL_upper_blind = 0.3
    #if 'Matrix' in plotname:
    #    x_upper_cut = MSCL_upper_blind
    #    y_upper_cut = MSCW_upper_blind
    #    x_lower_cut = MSCL_lower_blind
    #    y_lower_cut = MSCW_lower_blind
    #    x1, y1 = [ x_lower_cut,  x_lower_cut],  [ y_lower_cut,  y_upper_cut]
    #    x2, y2 = [ x_lower_cut,  x_upper_cut],  [ y_upper_cut,  y_upper_cut]
    #    x3, y3 = [ x_upper_cut,  x_upper_cut],  [ y_upper_cut,  y_lower_cut]
    #    x4, y4 = [ x_upper_cut,  x_lower_cut],  [ y_lower_cut,  y_lower_cut]
    #    plt.plot(x1, y1, x2, y2, x3, y3, x4, y4, color='k')

    divider = make_axes_locatable(axbig)
    cax = divider.append_axes("bottom", size="5%", pad=0.7)
    cbar = fig.colorbar(im,orientation="horizontal",cax=cax)
    cbar.set_label(label_z)
    
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()

def ConvertGalacticToRaDec(l, b):
    my_sky = SkyCoord(l*my_unit.deg, b*my_unit.deg, frame='galactic')
    return my_sky.icrs.ra.deg, my_sky.icrs.dec.deg

def GetVelocitySpectrum(map_file, roi_lon, roi_lat, roi_inner_ring, roi_outer_ring):

    roi_ra, roi_dec = ConvertGalacticToRaDec(roi_lon,roi_lat)
    print ('GetVelocitySpectrum, roi_ra = %0.1f, roi_dec = %0.1f'%(roi_ra,roi_dec))

    hdu = fits.open(map_file)[0]
    wcs = WCS(hdu.header)
    image_data = hdu.data
    dimensions = image_data.shape
    vel_dim = dimensions[2]
    lon_dim = dimensions[1]
    lat_dim = dimensions[0]

    lon_max = roi_lon-2.0
    lon_min = roi_lon+2.0
    lat_max = roi_lat+2.0
    lat_min = roi_lat-2.0
    pixs_min = wcs.all_world2pix(0.,lon_min,lat_min,1)
    pixs_max = wcs.all_world2pix(0.,lon_max,lat_max,1)
    pixs_lon_min = int(pixs_min[1])
    pixs_lat_min = int(pixs_min[2])
    pixs_lon_max = int(pixs_max[1])
    pixs_lat_max = int(pixs_max[2])

    vel_axis = []
    column_density = []
    for vel_pix in range(0,vel_dim):
        world_coord = wcs.all_pix2world(vel_pix,0,0,1) 
        velocity = world_coord[0]
        vel_axis += [velocity]
        total_pix = 0.
        avg_density = 0.
        for lon_pix in range(pixs_lon_min,pixs_lon_max):
            for lat_pix in range(pixs_lat_min,pixs_lat_max):
                world_coord = wcs.all_pix2world(vel_pix,lon_pix,lat_pix,1) 
                velocity = world_coord[0]
                lon = world_coord[1]
                lat = world_coord[2]
                distance = pow(pow(lon-roi_lon,2)+pow(lat-roi_lat,2),0.5)
                if distance<roi_inner_ring: continue
                if distance>roi_outer_ring: continue
                avg_density += image_data[lat_pix,lon_pix,vel_pix]
                total_pix += 1.
        avg_density = avg_density/total_pix
        column_density += [avg_density]

    return vel_axis, column_density

def GetGalfaHIVelocitySpectrum(map_file, roi_lon, roi_lat, roi_inner_ring, roi_outer_ring):

    vel_min = 0.
    vel_max = 50.*1e3

    roi_ra, roi_dec = ConvertGalacticToRaDec(roi_lon,roi_lat)
    print ('GetVelocitySpectrum, roi_ra = %0.1f, roi_dec = %0.1f'%(roi_ra,roi_dec))

    hdu = fits.open(map_file)[0]
    wcs = WCS(hdu.header)
    image_data = hdu.data
    dimensions = image_data.shape
    vel_dim = dimensions[0]
    dec_dim = dimensions[1]
    ra_dim = dimensions[2]

    ra_max = roi_ra-roi_outer_ring
    ra_min = roi_ra+roi_outer_ring
    dec_max = roi_dec+roi_outer_ring
    dec_min = roi_dec-roi_outer_ring
    pixs_min = wcs.all_world2pix(ra_min,dec_min, vel_min,1)
    pixs_max = wcs.all_world2pix(ra_max,dec_max, vel_max,1)
    pixs_ra_min = int(pixs_min[0])
    pixs_dec_min = int(pixs_min[1])
    pixs_vel_min = int(pixs_min[2])
    pixs_ra_max = int(pixs_max[0])
    pixs_dec_max = int(pixs_max[1])
    pixs_vel_max = int(pixs_max[2])

    vel_axis = []
    column_density = []
    for vel_pix in range(pixs_vel_min,pixs_vel_max):
        world_coord = wcs.all_pix2world(0,0,vel_pix,1) 
        velocity = world_coord[2]
        vel_axis += [velocity]
        total_pix = 0.
        avg_density = 0.
        for ra_pix in range(pixs_ra_min,pixs_ra_max):
            for dec_pix in range(pixs_dec_min,pixs_dec_max):
                if ra_pix<0: continue
                if dec_pix<0: continue
                if ra_pix>=image_data[:,:,:].shape[2]: continue
                if dec_pix>=image_data[:,:,:].shape[1]: continue
                world_coord = wcs.all_pix2world(ra_pix,dec_pix,vel_pix,1) 
                velocity = world_coord[2]
                ra = world_coord[0]
                dec = world_coord[1]
                distance = pow(pow(ra-roi_ra,2)+pow(dec-roi_dec,2),0.5)
                if distance<roi_inner_ring: continue
                if distance>roi_outer_ring: continue
                avg_density += image_data[vel_pix,dec_pix,ra_pix]
                total_pix += 1.
        avg_density = avg_density/total_pix
        column_density += [avg_density]

    return vel_axis, column_density

