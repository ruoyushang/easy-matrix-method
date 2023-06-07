
import os
import sys

SMI_DIR = os.environ['SMI_DIR']
SMI_BIN = os.environ['SMI_BIN']
#CONDA_DIR = os.environ['CONDA_DIR']

print (SMI_BIN)

isGamma2 = True

source = []
for_syst = []
for_imposter = []

folder = '%s'%(SMI_BIN)
rank = 3
MSCW_cut = 0.3
MSCL_cut = 0.3

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

PrepareSource('H1426_V6_OFF') 
PrepareSource('H1426_V5_OFF') 
PrepareSource('UrsaMajorII_V6_OFF') 
PrepareSource('UrsaMajorII_V5_OFF') 
PrepareSource('UrsaMinor_V6_OFF') 
PrepareSource('UrsaMinor_V5_OFF') 
PrepareSource('RGB_J0710_p591_V6_OFF') 
PrepareSource('RGB_J0710_p591_V5_OFF') 
PrepareSource('1ES0229_V6_OFF') 
PrepareSource('1ES0229_V5_OFF') 
PrepareSource('PKS1424_V6_OFF') 
PrepareSource('PKS1424_V5_OFF') 
PrepareSource('PG1553_V6_OFF') 
PrepareSource('PG1553_V5_OFF') 
PrepareSource('3C273_V6_OFF') 
PrepareSource('3C273_V5_OFF') 
PrepareSource('Segue1_V6_OFF') 
PrepareSource('Segue1_V5_OFF') 
PrepareSource('NGC1275_V6_OFF') 
PrepareSource('NGC1275_V5_OFF') 
PrepareSource('OJ287_V6_OFF') 
PrepareSource('OJ287_V5_OFF') 
PrepareSource('Draco_V6_OFF') 
PrepareSource('Draco_V5_OFF') 
PrepareSource('BLLac_V6_OFF') 
PrepareSource('BLLac_V5_OFF') 
PrepareSource('3C264_V6_OFF') 
PrepareSource('3C264_V5_OFF') 
PrepareSource('1ES0502_V6_OFF') 
PrepareSource('1ES0502_V5_OFF') 
PrepareSource('M82_V6_OFF') 
PrepareSource('M82_V5_OFF') 
PrepareSource('1ES0414_V6_OFF') 
PrepareSource('1ES0414_V5_OFF') 
PrepareSource('1ES1011_V6_OFF') 
PrepareSource('1ES1011_V5_OFF') 
PrepareSource('1ES0647_V6_OFF') 
PrepareSource('1ES0647_V5_OFF') 

PrepareSource('CrabNebula_elev_80_90_V6_OFF')
PrepareSource('CrabNebula_elev_70_80_V6_OFF')
PrepareSource('CrabNebula_elev_60_70_V6_OFF')
PrepareSource('CrabNebula_elev_50_60_V6_OFF')
PrepareSource('CrabNebula_elev_40_50_V6_OFF')
PrepareSource('CrabNebula_elev_80_90_V5_OFF')
PrepareSource('CrabNebula_elev_70_80_V5_OFF')
PrepareSource('CrabNebula_elev_60_70_V5_OFF')
PrepareSource('CrabNebula_elev_50_60_V5_OFF')
PrepareSource('CrabNebula_elev_40_50_V5_OFF')

PrepareSource('CrabNebula_elev_80_90_V6_ON')
PrepareSource('CrabNebula_elev_70_80_V6_ON')
PrepareSource('CrabNebula_elev_60_70_V6_ON')
PrepareSource('CrabNebula_elev_50_60_V6_ON')
PrepareSource('CrabNebula_elev_40_50_V6_ON')
PrepareSource('CrabNebula_elev_80_90_V5_ON')
PrepareSource('CrabNebula_elev_70_80_V5_ON')
PrepareSource('CrabNebula_elev_60_70_V5_ON')
PrepareSource('CrabNebula_elev_50_60_V5_ON')
PrepareSource('CrabNebula_elev_40_50_V5_ON')

PrepareSource('CrabNebula_1p0wobble_V6_ON')
PrepareSource('CrabNebula_1p5wobble_V6_ON')

imposter_analyses = []

imposter_analyses += ['Sky_RA55Dec53_V6']
imposter_analyses += ['1ES0229_V6']
imposter_analyses += ['1ES0229_V5']
imposter_analyses += ['PSR_J1907_p0602_V6']
imposter_analyses += ['PSR_J1907_p0602_V5']
imposter_analyses += ['PSR_J2021_p4026_V6'] # Gamma Cygni
imposter_analyses += ['PSR_J2021_p4026_V5']
imposter_analyses += ['SNR_G189_p03_V6'] # IC 443
imposter_analyses += ['SNR_G189_p03_V5']
imposter_analyses += ['Geminga_V6'] 
imposter_analyses += ['Geminga_V5']
imposter_analyses += ['PSR_J2021_p3651_V6'] # Dragonfly
imposter_analyses += ['PSR_J2021_p3651_V5']
imposter_analyses += ['SNR_G067p6_p0p9_V6']
imposter_analyses += ['CTB109_V6']
imposter_analyses += ['CTB109_V5']

#imposter_analyses += ['SgrA_V6']


for analysis in range(0,len(imposter_analyses)):
    PrepareSource('%s_ON'%(imposter_analyses[analysis])) 
    for imposter in range(0,5):
        PrepareSource('%s_Imposter%s'%(imposter_analyses[analysis],imposter+1),imposter_idx=imposter+1)



job_counts = 0
for s in range(0,len(source)):
    job_counts += 1
    file = open("run/run_Netflix1_%s.sh"%(source[s]),"w") 
    if not isGamma2:
        file.write('#### submit_job.sh START ####\n')
        file.write('#!/bin/bash\n')
        file.write('#$ -cwd\n')
        file.write('# error = Merged with joblog\n')
        file.write('#$ -o joblog.$JOB_ID\n')
        file.write('#$ -j y\n')
        file.write('## Edit the line below as needed:\n')
        file.write('#$ -l h_rt=2:00:00,h_data=4G\n')
        file.write('## Modify the parallel environment\n')
        file.write('## and the number of cores as needed:\n')
        file.write('#$ -pe shared 1\n')
        file.write('# Email address to notify\n')
        file.write('#$ -M $USER@mail\n')
        file.write('# Notify when\n')
        #file.write('#$ -m bea\n')
        file.write('#$ -m e\n')
        file.write('. /u/local/Modules/default/init/modules.sh\n')
        file.write('module load cern_root/6.12.06\n')
        #file.write('module load anaconda3\n')
        #file.write('source %s/etc/profile.d/conda.sh\n'%(CONDA_DIR))
        #file.write('conda activate root_env\n')
    file.write('cd %s\n'%(SMI_DIR))
    #file.write('source %s/setup_env.sh\n'%(SMI_DIR))
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
qfile = open("run/qsub_Netflix1.sh","w") 
for s in range(0,len(source)):
    job_counts += 1
    qfile.write('qsub -V -N job_ana_%s run_Netflix1_%s.sh\n'%(source[s],source[s]))
    qfile.write('sleep 15s\n')
qfile.close() 


