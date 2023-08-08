
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

hist_sci_erecerr = ROOT.TH1D("hist_sci_erecerr","",20,0.,2.)
hist_bkg_erecerr = ROOT.TH1D("hist_bkg_erecerr","",20,0.,2.)
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
hist_2d_sci_mva_mscw = ROOT.TH2D("hist_2d_sci_mva_mscw","",20,-1.,1.,20,-1.,2.)
hist_2d_sci_mva_mscl = ROOT.TH2D("hist_2d_sci_mva_mscl","",20,-1.,1.,20,-1.,2.)
hist_2d_bkg_mva_mscw = ROOT.TH2D("hist_2d_bkg_mva_mscw","",20,-1.,1.,20,-1.,2.)
hist_2d_bkg_mva_mscl = ROOT.TH2D("hist_2d_bkg_mva_mscl","",20,-1.,1.,20,-1.,2.)

runnumber = []
#runnumber += [66730]
#runnumber += [66731]
#runnumber += [66732]
#runnumber += [66733]
#runnumber += [66786]
#runnumber += [66787]
#runnumber += [66788]
#runnumber += [66789]
#runnumber += [70489]
#runnumber += [70490]
#runnumber += [70532]
#runnumber += [70533]
#runnumber += [75888]
#runnumber += [75904]
#runnumber += [76212]
#runnumber += [76213]
#runnumber += [79157]
#runnumber += [79158]
#runnumber += [83655]
#runnumber += [83656]
#runnumber += [83702]
#runnumber += [83703]
#runnumber += [84065]
#runnumber += [85000]
#runnumber += [85455]
#runnumber += [91106]
#runnumber += [91127]
#runnumber += [91128]
#runnumber += [91152]
#runnumber += [91153]
#runnumber += [96011]
#runnumber += [96706]
#runnumber += [96786]
#runnumber += [97646]
#runnumber += [97663]
#runnumber += [97730]
runnumber += [101692]
runnumber += [101717]
runnumber += [101718]
runnumber += [101743]
runnumber += [101744]
runnumber += [101777]
runnumber += [101923]
runnumber += [101963]
runnumber += [103089]
runnumber += [103136]
runnumber += [103149]
runnumber += [103171]
runnumber += [103172]
runnumber += [103491]
runnumber += [103492]
runnumber += [103542]
runnumber += [103543]
runnumber += [103544]
runnumber += [103569]
runnumber += [103591]
runnumber += [103792]
runnumber += [103793]
runnumber += [103794]

for run in range(0,len(runnumber)):
    filename = '/gamma_raid/userspace/rshang/analysis/Results/v490/%s.anasum.root'%(runnumber[run])
    InputFile = ROOT.TFile(filename)
    DL3EventTree = InputFile.Get("run_%s/stereo/DL3EventTree"%(runnumber[run]))
    total_entries = DL3EventTree.GetEntries()
    print ('DL3EventTree.GetEntries() = %s'%(total_entries))
    
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
    
        min_images = 3
        if NImages<min_images: continue
        if Erec<0.2: continue
        if Erec>0.5: continue
        if pow(Xoff*Xoff+Yoff*Yoff,0.5)>1.5: continue
        if EmissionHeight<5.: continue
        if EmissionHeight>15.: continue
        #if Erec<1.0: continue
        #if NImages<4: continue
        #if MVA<0.: continue
    
        #if MSCW>0.5: continue
        #if MSCL>0.7: continue
        #if MSCW>1.0: continue
        #if MSCW<-1.0: continue
        #if MSCL>1.4: continue
        #if MSCL<-1.0: continue
    
        dist_to_crab = pow(pow(RA-83.633,2)+pow(DEC-22.014,2),0.5)
        if dist_to_crab<0.15:
            hist_sci_nimage.Fill(NImages)
            hist_sci_height.Fill(EmissionHeight)
            hist_sci_rcore.Fill(RCore)
            hist_sci_erecerr.Fill(Erec_Err)
            hist_sci_mva.Fill(MVA)
            hist_sci_mscw.Fill(MSCW)
            hist_2d_sci_mva_mscw.Fill(MVA,MSCW)
            hist_2d_sci_mva_mscl.Fill(MVA,MSCL)
            if MSCW<0.6:
                hist_sci_mscl.Fill(MSCL)
        else:
            hist_bkg_nimage.Fill(NImages)
            hist_bkg_height.Fill(EmissionHeight)
            hist_bkg_rcore.Fill(RCore)
            hist_bkg_erecerr.Fill(Erec_Err)
            hist_bkg_mva.Fill(MVA)
            hist_bkg_mscw.Fill(MSCW)
            hist_2d_bkg_mva_mscw.Fill(MVA,MSCW)
            hist_2d_bkg_mva_mscl.Fill(MVA,MSCL)
            if MSCW<0.6:
                hist_bkg_mscl.Fill(MSCL)


histos = [hist_bkg_erecerr,hist_sci_erecerr]
do_1dplots(histos,'Eerr')
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

histos = [hist_2d_bkg_mva_mscw]
do_2dplots(histos,'MVA_MSCW_BKG')
histos = [hist_2d_bkg_mva_mscl]
do_2dplots(histos,'MVA_MSCL_BKG')
histos = [hist_2d_sci_mva_mscw]
do_2dplots(histos,'MVA_MSCW_SCI')
histos = [hist_2d_sci_mva_mscl]
do_2dplots(histos,'MVA_MSCL_SCI')
