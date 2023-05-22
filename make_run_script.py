
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

def PrepareSource(source_name):
    global source
    global for_syst
    global for_imposter
    if 'Imposter' in source_name:
        source += [source_name]
        for_syst += [False]
        for_imposter += [True]
    else:
        source += [source_name]
        if '_ON' in source_name:
            for_syst += [False]
        else:
            for_syst += [True]
        for_imposter += [False]

#PrepareSource('UrsaMajorIIV6_OFF') 
#PrepareSource('UrsaMinorV6_OFF') 
#PrepareSource('UrsaMinorV5_OFF') 
#PrepareSource('RGB_J0710_p591_V6_OFF') 
#PrepareSource('RGB_J0710_p591_V5_OFF') 
#PrepareSource('Sky_RA38Dec20_V6_OFF') 
#PrepareSource('Sky_RA38Dec20_V5_OFF') 
#PrepareSource('PKS1424V6_OFF') 
#PrepareSource('PKS1424V5_OFF') 
#PrepareSource('PG1553V6_OFF') 
#PrepareSource('PG1553V5_OFF') 
#PrepareSource('3C273V6_OFF') 
#PrepareSource('3C273V5_OFF') 
#PrepareSource('Segue1V6_OFF') 
#PrepareSource('Segue1V5_OFF') 
#PrepareSource('NGC1275V6_OFF') 
#PrepareSource('H1426V6_OFF') 
#PrepareSource('OJ287V6_OFF') 
#PrepareSource('DracoV6_OFF') 
#PrepareSource('DracoV5_OFF') 
#PrepareSource('BLLacV6_OFF') 
#PrepareSource('BLLacV5_OFF') 
#PrepareSource('3C264V6_OFF') 
#PrepareSource('1ES0502V6_OFF') 
#PrepareSource('1ES0502V5_OFF') 
#PrepareSource('M82V6_OFF') 
#PrepareSource('M82V5_OFF') 
#PrepareSource('1ES0414V5_OFF') 
PrepareSource('1ES1011_V6_OFF') 
#PrepareSource('1ES0647_V6_OFF') 

#PrepareSource('CrabNebula_elev_70_90_V6_OFF')
#PrepareSource('CrabNebula_elev_60_70_V6_OFF')
#PrepareSource('CrabNebula_elev_50_60_V6_OFF')
#PrepareSource('CrabNebula_elev_40_50_V6_OFF')
#PrepareSource('CrabNebula_elev_70_90_V6_ON')
#PrepareSource('CrabNebula_elev_60_70_V6_ON')
#PrepareSource('CrabNebula_elev_50_60_V6_ON')
#PrepareSource('CrabNebula_elev_40_50_V6_ON')

imposter_analyses = []

#imposter_analyses += ['CrabRHV_V6']
#imposter_analyses += ['LHAASO_J0621_p3755_RHV_V6']

#imposter_analyses += ['1ES0647V6']
#imposter_analyses += ['1ES1011V6']
#imposter_analyses += ['1ES0414V5']
#imposter_analyses += ['UrsaMinorV6']
#imposter_analyses += ['UrsaMinorV5']
#imposter_analyses += ['RGB_J0710_p591_V6']
#imposter_analyses += ['RGB_J0710_p591_V5']
#imposter_analyses += ['Sky_RA38Dec20_V6'] #1ES0229
#imposter_analyses += ['Sky_RA38Dec20_V5'] #1ES0229
#imposter_analyses += ['3C273V6']
#imposter_analyses += ['3C273V5']
#imposter_analyses += ['Segue1V6']
#imposter_analyses += ['Segue1V5']
#imposter_analyses += ['H1426V6']
#imposter_analyses += ['NGC1275V6']
#imposter_analyses += ['DracoV6']
#imposter_analyses += ['DracoV5']
#imposter_analyses += ['UrsaMajorIIV6']
#imposter_analyses += ['PKS1424V6']
#imposter_analyses += ['PKS1424V5']
#imposter_analyses += ['PG1553V6']
#imposter_analyses += ['PG1553V5']
#imposter_analyses += ['OJ287V6']
#imposter_analyses += ['BLLacV6']
#imposter_analyses += ['BLLacV5']
#imposter_analyses += ['3C264V6']
#imposter_analyses += ['1ES0502V6']
#imposter_analyses += ['1ES0502V5']
#imposter_analyses += ['M82V6']
#imposter_analyses += ['M82V5']

#imposter_analyses += ['Crab_Offset_1p0_V6']
#imposter_analyses += ['Crab_Offset_1p5_V6']
#imposter_analyses += ['CrabNebula_elev_70_90_V6']
#imposter_analyses += ['CrabNebula_elev_60_70_V6']
#imposter_analyses += ['CrabNebula_elev_50_60_V6']
#imposter_analyses += ['CrabNebula_elev_40_50_V6']

#imposter_analyses += ['IC443HotSpotV6']
#imposter_analyses += ['IC443HotSpotV5']
#imposter_analyses += ['MGRO_J2019_V6']
#imposter_analyses += ['MGRO_J2019_V5']
#imposter_analyses += ['BoomerangV6']
#imposter_analyses += ['BoomerangV5']
#imposter_analyses += ['TychoV6']
#imposter_analyses += ['TychoV5']
#imposter_analyses += ['LHAASO_J2032_V5']
#imposter_analyses += ['LHAASO_J2032_V6']
#imposter_analyses += ['LHAASO_J2032_Baseline_V5']
#imposter_analyses += ['LHAASO_J2032_Baseline_V6']
#imposter_analyses += ['LHAASO_J2032_Fall2017_V6']

#imposter_analyses += ['MGRO_J1908_V6']
#imposter_analyses += ['MGRO_J1908_V5']
#imposter_analyses += ['SS433_V6']
#imposter_analyses += ['SS433_V5']
#imposter_analyses += ['SS433Half1_V6']
#imposter_analyses += ['SS433Half1_V5']
#imposter_analyses += ['SS433Half2_V6']
#imposter_analyses += ['SS433Half2_V5']
#imposter_analyses += ['PSR_J1907_p0602_V6']
#imposter_analyses += ['PSR_J1907_p0602_V5']
#imposter_analyses += ['PSR_J2032_p4127_Fall2017_V6']
#imposter_analyses += ['PSR_J2032_p4127_Baseline_V6']
#imposter_analyses += ['PSR_J2032_p4127_Baseline_V5']
#imposter_analyses += ['PSR_J2032_p4127_All_V6']
#imposter_analyses += ['PSR_J2032_p4127_All_V5']
#imposter_analyses += ['PSR_J2229_p6114_V6'] # boomerang
#imposter_analyses += ['PSR_J2229_p6114_V5']
#imposter_analyses += ['PSR_J2021_p3651_V6'] # Dragonfly
#imposter_analyses += ['PSR_J2021_p3651_V5']
#imposter_analyses += ['PSR_J2021_p4026_V6'] # Gamma Cygni
#imposter_analyses += ['PSR_J2021_p4026_V5']
#imposter_analyses += ['PSR_J1856_p0245_V6']
#imposter_analyses += ['PSR_J1856_p0245_V5']
#imposter_analyses += ['PSR_J1930_p1852_V6']
#imposter_analyses += ['PSR_J1928_p1746_V6']
#imposter_analyses += ['GemingaV6']
#imposter_analyses += ['GemingaV5']
#imposter_analyses += ['SNR_G189_p03_V6'] # IC 443
#imposter_analyses += ['SNR_G189_p03_V5']
#imposter_analyses += ['PSR_J1954_p2836_V6'] # only detected at ~200 GeV, age = 6.94e+04 yr, Edot = 1.0e+36 erg/s
#imposter_analyses += ['PSR_J2238_p5903_V6']
#imposter_analyses += ['PSR_J0205_p6449_V6']

#imposter_analyses += ['PSR_B2255_p58_V6']
#imposter_analyses += ['PSR_B2255_p58_V5']
#imposter_analyses += ['PSR_B1823_m13_V6']
#imposter_analyses += ['PSR_J0633_p0632_V6']
#imposter_analyses += ['PSR_J0633_p0632_V5']
#imposter_analyses += ['PSR_J1954_p2836_V6']
#imposter_analyses += ['SNR_G109_m01_V6']
#imposter_analyses += ['SNR_G109_m01_V5']
#imposter_analyses += ['SNR_G111_m02_V6']
#imposter_analyses += ['PSR_J2022_p3842_V6']
#imposter_analyses += ['PSR_J2022_p3842_V5']
#imposter_analyses += ['PSR_B2127_p11_V5']
#imposter_analyses += ['SNR_G073_p00_V6']
#imposter_analyses += ['SNR_G073_p00_V5']
#imposter_analyses += ['PSR_J1846_m0258_V6']
#imposter_analyses += ['PSR_B1937_p21_V6']
#imposter_analyses += ['PSR_J1946_p2052_V6']
#imposter_analyses += ['PSR_J0517_p2212_V6']
#imposter_analyses += ['PSR_J0517_p2212_V5']

#imposter_analyses += ['V_V725_Tau_V6']

#imposter_analyses += ['Galactic_RA304Dec36_V6']
#imposter_analyses += ['Galactic_RA345Dec59_V6']
#imposter_analyses += ['Galactic_RA09Dec60_V6']
#imposter_analyses += ['Sky_RA55Dec53_V6']

for analysis in range(0,len(imposter_analyses)):
    PrepareSource('%s_ON'%(imposter_analyses[analysis])) 
    for imposter in range(0,5):
        PrepareSource('%s_Imposter%d'%(imposter_analyses[analysis],imposter+1))



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
    file.write('cp NetflixParameters.h %s/%s\n'%(folder,source[s]))
    file.write('cp FillHistograms.C %s/%s\n'%(folder,source[s]))
    file.write('cd %s/%s\n'%(folder,source[s]))
    file.write('rm *_C*\n')
    if for_syst[s]:
        file.write("root -b -l -q 'FillHistograms.C+(\"%s\",false,false)'\n"%(source[s])) 
    else:
        if for_imposter[s]:
            file.write("root -b -l -q 'FillHistograms.C+(\"%s\",true,true)'\n"%(source[s])) 
        else:
            file.write("root -b -l -q 'FillHistograms.C+(\"%s\",true,false)'\n"%(source[s])) 
    file.close() 

job_counts = 0
qfile = open("run/qsub_Netflix1.sh","w") 
for s in range(0,len(source)):
    job_counts += 1
    qfile.write('qsub -V -N job_%s run_Netflix1_%s.sh\n'%(source[s],source[s]))
    qfile.write('sleep 30s\n')
qfile.close() 


