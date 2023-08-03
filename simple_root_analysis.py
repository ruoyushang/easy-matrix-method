
import sys,ROOT
import array
import math
from array import *
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt

def do_1dplots(histos,name):
    c = ROOT.TCanvas('c1', 'c1') 
    c.cd()
    histos[0].SetMinimum(1.)
    histos[0].Draw()
    histos[1].SetLineColor(kRed)
    histos[1].Draw("same")
    c.Modified()
    c.Update()
    c.SetLogy()
    c.SaveAs('output_plots/%s.png'%(name))

def do_2dplots(histos,name):
    c = ROOT.TCanvas('c1', 'c1') 
    c.cd()
    histos[0].Draw("COLZ")
    c.Modified()
    c.Update()
    c.SaveAs('output_plots/%s.png'%(name))

runnumber = 65312 
filename = '/gamma_raid/userspace/rshang/analysis/Results/v490/%s.anasum.root'%(runnumber)
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

hist_sci_rcore = ROOT.TH1D("hist_sci_rcore","",20,0.,1000.)
hist_bkg_rcore = ROOT.TH1D("hist_bkg_rcore","",20,0.,1000.)
hist_sci_height = ROOT.TH1D("hist_sci_height","",20,0.,30.)
hist_bkg_height = ROOT.TH1D("hist_bkg_height","",20,0.,30.)
hist_sci_mscl = ROOT.TH1D("hist_sci_mscl","",20,-1.0,2.0)
hist_bkg_mscl = ROOT.TH1D("hist_bkg_mscl","",20,-1.0,2.0)
hist_sci_mscw = ROOT.TH1D("hist_sci_mscw","",20,-1.0,2.0)
hist_bkg_mscw = ROOT.TH1D("hist_bkg_mscw","",20,-1.0,2.0)
hist_sci_nimage = ROOT.TH1D("hist_sci_nimage","",4,1,5)
hist_bkg_nimage = ROOT.TH1D("hist_bkg_nimage","",4,1,5)
hist_sci_mva = ROOT.TH1D("hist_sci_mva","",20,-1.,1.)
hist_bkg_mva = ROOT.TH1D("hist_bkg_mva","",20,-1.,1.)
hist_2d_mva_mscw = ROOT.TH2D("hist_2d_mva_mscw","",20,-1.,1.,20,-1.,2.)

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
    Xoff = DL3EventTree.Xoff
    Yoff = DL3EventTree.Yoff
    MVA = DL3EventTree.MVA
    RCore = pow(XCore*XCore+YCore*YCore,0.5)

    if Erec<0.3: continue
    if pow(Xoff*Xoff+Yoff*Yoff,0.5)>1.5: continue
    if EmissionHeight<5.: continue
    if EmissionHeight>15.: continue
    #if Erec<1.0: continue
    #if NImages<4: continue

    #if MSCW>0.5: continue
    #if MSCL>0.7: continue
    #if MSCW>1.0: continue
    #if MSCW<-1.0: continue
    #if MSCL>1.4: continue
    #if MSCL<-1.0: continue

    min_images = 3
    dist_to_crab = pow(pow(RA-83.633,2)+pow(DEC-22.014,2),0.5)
    if dist_to_crab<0.10:
        hist_sci_nimage.Fill(NImages)
        hist_sci_height.Fill(EmissionHeight)
        hist_sci_rcore.Fill(RCore)
        if NImages>=min_images:
            hist_sci_mva.Fill(MVA)
            hist_sci_mscw.Fill(MSCW)
            if MSCW<0.6:
                hist_sci_mscl.Fill(MSCL)
    else:
        hist_bkg_nimage.Fill(NImages)
        hist_bkg_height.Fill(EmissionHeight)
        hist_bkg_rcore.Fill(RCore)
        if NImages>=min_images:
            hist_bkg_mva.Fill(MVA)
            hist_bkg_mscw.Fill(MSCW)
            hist_2d_mva_mscw.Fill(MVA,MSCW)
            if MSCW<0.6:
                hist_bkg_mscl.Fill(MSCL)


histos = [hist_bkg_rcore,hist_sci_rcore]
do_1dplots(histos,'Rcore')
histos = [hist_bkg_height,hist_sci_height]
do_1dplots(histos,'Height')
histos = [hist_bkg_nimage,hist_sci_nimage]
do_1dplots(histos,'NImage')
histos = [hist_bkg_mscl,hist_sci_mscl]
do_1dplots(histos,'MSCL')
histos = [hist_bkg_mscw,hist_sci_mscw]
do_1dplots(histos,'MSCW')
histos = [hist_bkg_mva,hist_sci_mva]
do_1dplots(histos,'MVA')

histos = [hist_2d_mva_mscw]
do_2dplots(histos,'MVA_MSCW')
