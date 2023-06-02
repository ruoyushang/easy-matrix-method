
import sys,ROOT
import array
import math
from array import *
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt

runnumber = 97646
filename = '/gamma_raid/userspace/rshang/analysis/Results/v487/%s.anasum.root'%(runnumber)
InputFile = ROOT.TFile(filename)
DL3EventTree = InputFile.Get("run_%s/stereo/DL3EventTree"%(runnumber))
total_entries = DL3EventTree.GetEntries()
print ('DL3EventTree.GetEntries() = %s'%(total_entries))

# Variables in the tree:
# Erec            
# Erec_Err        
# XCore           
# YCore           
# Xderot          
# Yderot          
# NImages         
# ImgSel          
# MeanPedvar      
# MSCW            
# MSCL            
# RA              
# DEC             
# Az              
# El              
# EmissionHeight  
# Xoff            
# Yoff            
# Acceptance      

sig_height = []
bkg_height = []
sig_mscw = []
bkg_mscw = []
sig_mscl = []
bkg_mscl = []
sig_eng = []
bkg_eng = []
sig_eng_err = []
bkg_eng_err = []
sig_core = []
bkg_core = []
for entry in range(0,total_entries):

    DL3EventTree.GetEntry(entry)

    NImages = DL3EventTree.NImages
    Erec = DL3EventTree.Erec
    Erec_Err = DL3EventTree.Erec_Err
    RA = DL3EventTree.RA
    DEC = DL3EventTree.DEC
    EmissionHeight = DL3EventTree.EmissionHeight
    MSCW = DL3EventTree.MSCW
    MSCL = DL3EventTree.MSCL
    XCore = DL3EventTree.XCore
    YCore = DL3EventTree.YCore

    if Erec<0.25: continue
    if NImages<3: continue
    #if MSCW>0.5: continue
    #if MSCL>0.7: continue
    #if MSCW>1.0: continue
    #if MSCW<-1.0: continue
    #if MSCL>1.4: continue
    #if MSCL<-1.0: continue

    dist_to_crab = pow(pow(RA-83.633,2)+pow(DEC-22.014,2),0.5)
    if dist_to_crab<0.10:
        sig_height += [EmissionHeight]
        sig_mscw += [MSCW]
        sig_mscl += [MSCL]
        sig_eng += [Erec]
        sig_eng_err += [Erec_Err]
        sig_core += [pow(XCore*XCore+YCore*YCore,0.5)]
    else:
        bkg_height += [EmissionHeight]
        bkg_mscw += [MSCW]
        bkg_mscl += [MSCL]
        bkg_eng += [Erec]
        bkg_eng_err += [Erec_Err]
        bkg_core += [pow(XCore*XCore+YCore*YCore,0.5)]


fig, ax = plt.subplots()
figsize_x = 8
figsize_y = 6
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

fig.clf()
fig.set_figheight(8)
fig.set_figwidth(8)
axbig = fig.add_subplot()
axbig.scatter(bkg_core,bkg_height,color='r',alpha=0.2)
axbig.scatter(sig_core,sig_height,color='b',alpha=0.2)
axbig.set_xlabel('R core')
axbig.set_ylabel('Height')
fig.savefig("output_plots/simple_analysis_core_height_plot.png")
axbig.remove()

fig.clf()
fig.set_figheight(8)
fig.set_figwidth(8)
axbig = fig.add_subplot()
axbig.scatter(bkg_mscw,bkg_mscl,color='r',alpha=0.2)
axbig.scatter(sig_mscw,sig_mscl,color='b',alpha=0.2)
axbig.set_xlabel('MSCW')
axbig.set_ylabel('MSCL')
fig.savefig("output_plots/simple_analysis_mscw_mscl_plot.png")
axbig.remove()

fig.clf()
fig.set_figheight(8)
fig.set_figwidth(8)
axbig = fig.add_subplot()
axbig.scatter(bkg_eng,bkg_eng_err,color='r',alpha=0.2)
axbig.scatter(sig_eng,sig_eng_err,color='b',alpha=0.2)
axbig.set_xlabel('E')
axbig.set_ylabel('E err')
axbig.set_xscale('log')
fig.savefig("output_plots/simple_analysis_eng_err_plot.png")
axbig.remove()

