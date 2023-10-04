
import os
import sys

SMI_DIR = os.environ['SMI_DIR']
SMI_BIN = os.environ['SMI_BIN']
#CONDA_DIR = os.environ['CONDA_DIR']

print (SMI_BIN)

#isTraining = True
isTraining = False

source = []
for_syst = []
for_imposter = []

folder = '%s'%(SMI_BIN)

epoch_list = []
epoch_list += ['V4']
epoch_list += ['V5']
epoch_list += ['V6']

def PrepareSource(source_name, imposter_idx=0):
    global source
    global for_syst
    global for_imposter
    if imposter_idx>0:
        source += [source_name]
        for_syst += [False]
        for_imposter += [imposter_idx]
    else:
        source += [source_name]
        if '_ON' in source_name:
            for_syst += [False]
        else:
            for_syst += [True]
        for_imposter += [0]


imposter_analyses = []

#imposter_analyses += ['UrsaMajorII'] 
#imposter_analyses += ['UrsaMinor'] 
#imposter_analyses += ['RGB_J0710_p591'] 
#imposter_analyses += ['1ES0229'] 
#imposter_analyses += ['PKS1424'] 
#imposter_analyses += ['PG1553'] 
#imposter_analyses += ['3C273'] 
#imposter_analyses += ['Segue1'] 
#imposter_analyses += ['NGC1275'] 
#imposter_analyses += ['OJ287'] 
#imposter_analyses += ['Draco'] 
#imposter_analyses += ['BLLac'] 
#imposter_analyses += ['3C264'] 
#imposter_analyses += ['1ES0502'] 
#imposter_analyses += ['M82'] 
#imposter_analyses += ['1ES0647'] 
#imposter_analyses += ['1ES1011'] 
#imposter_analyses += ['1ES0414'] 
#imposter_analyses += ['H1426'] 
#
#imposter_analyses += ['CrabNebula_elev_80_90']
#imposter_analyses += ['CrabNebula_elev_70_80']
#imposter_analyses += ['CrabNebula_elev_60_70']
#imposter_analyses += ['CrabNebula_elev_50_60']
#imposter_analyses += ['CrabNebula_elev_40_50']
#imposter_analyses += ['CrabNebula_1p0wobble']
#imposter_analyses += ['CrabNebula_1p5wobble']

#imposter_analyses += ['PSR_J1856_p0245']
#imposter_analyses += ['PSR_J1907_p0602']
#imposter_analyses += ['PSR_J1928_p1746']
#
#imposter_analyses += ['SS433']
#imposter_analyses += ['SNR_G189_p03'] # IC 443
#imposter_analyses += ['PSR_J2021_p4026'] # Gamma Cygni
#imposter_analyses += ['Geminga'] 
#imposter_analyses += ['SNR_G150_p4']
imposter_analyses += ['PSR_J2021_p3651'] # Dragonfly
#imposter_analyses += ['PSR_J2032_p4127']

#imposter_analyses += ['Cas_A']
#imposter_analyses += ['CTA1']
#imposter_analyses += ['Tycho']


for analysis in range(0,len(imposter_analyses)):
    for epoch in range(0,len(epoch_list)):
        if isTraining:
            PrepareSource('%s_%s_OFF'%(imposter_analyses[analysis],epoch_list[epoch])) 
        else:
            PrepareSource('%s_%s_ON'%(imposter_analyses[analysis],epoch_list[epoch])) 
            for imposter in range(0,6):
                PrepareSource('%s_%s_Imposter%s'%(imposter_analyses[analysis],epoch_list[epoch],imposter+1),imposter_idx=imposter+1)



job_counts = 0
for s in range(0,len(source)):
    job_counts += 1
    file = open("run/ana_%s.sh"%(source[s]),"w") 
    file.write('cd %s\n'%(SMI_DIR))
    file.write('mkdir %s\n'%(folder))
    file.write('rm -r %s/%s\n'%(folder,source[s]))
    file.write('mkdir %s/%s\n'%(folder,source[s]))
    file.write('cp GetRunList.h %s/%s\n'%(folder,source[s]))
    file.write('cp ResetPublicVariables.C %s/%s\n'%(folder,source[s]))
    file.write('cp NetflixParameters.h %s/%s\n'%(folder,source[s]))
    file.write('cp FillHistograms.C %s/%s\n'%(folder,source[s]))
    file.write('cd %s/%s\n'%(folder,source[s]))
    file.write('rm *_C*\n')
    if for_syst[s]:
        file.write("root -b -l -q 'FillHistograms.C+(\"%s\",false,0)'\n"%(source[s])) 
    else:
        file.write("root -b -l -q 'FillHistograms.C+(\"%s\",true,%s)'\n"%(source[s],for_imposter[s])) 
    file.close() 

job_counts = 0
qfile = open("run/condor_ana.sh","w") 
for s in range(0,len(source)):
    job_counts += 1
    qfile.write('universe = vanilla \n')
    qfile.write('getenv = true \n')
    qfile.write('executable = /bin/bash \n')
    qfile.write('arguments = ana_%s.sh\n'%(source[s]))
    qfile.write('request_cpus = 1 \n')
    qfile.write('request_memory = 1024M \n')
    qfile.write('request_disk = 1024M \n')
    qfile.write('output = condor_ana_%s.out\n'%(source[s]))
    qfile.write('error = condor_ana_%s.err\n'%(source[s]))
    qfile.write('log = condor_ana_%s.log\n'%(source[s]))
    qfile.write('queue\n')
qfile.close() 

job_counts = 0
qfile = open("run/local_ana.sh","w") 
for s in range(0,len(source)):
    job_counts += 1
    qfile.write('sh ana_%s.sh\n'%(source[s]))
    qfile.write('sleep 2s\n')
qfile.close() 

