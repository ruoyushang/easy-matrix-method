
import os
import sys

SMI_DIR = os.environ['SMI_DIR']

input_params = []

input_params += [ ['1ES0647'               ,102.694 ,25.050 ,30 ,90] ]
input_params += [ ['1ES1011'               ,153.767 ,49.434 ,30 ,90] ]
input_params += [ ['1ES0414'               ,64.220  ,1.089  ,30 ,90] ]
input_params += [ ['M82'                   ,148.970 ,69.679 ,30 ,90] ]
input_params += [ ['1ES0502'               ,76.983  ,67.623 ,30 ,90] ]
input_params += [ ['3C264'                 ,176.271 ,19.606 ,30 ,90] ]
input_params += [ ['BLLac'                 ,330.680 ,42.277 ,30 ,90] ]
input_params += [ ['Draco'                 ,260.059 ,57.921 ,30 ,90] ]
input_params += [ ['OJ287'                 ,133.705 ,20.100 ,30 ,90] ]
input_params += [ ['H1426'                 ,217.136  ,42.673 ,30 ,90] ]
input_params += [ ['NGC1275'               ,49.950  ,41.512 ,30 ,90] ]
input_params += [ ['Segue1'                ,151.767 ,16.082 ,30 ,90] ]
input_params += [ ['3C273'                 ,187.277 ,2.05   ,30 ,90] ]
input_params += [ ['PG1553'                ,238.936 ,11.195 ,30 ,90] ]
input_params += [ ['PKS1424'               ,216.750 ,23.783 ,30 ,90] ]
input_params += [ ['1ES0229'               ,38.222  ,20.273 ,30 ,90] ]
input_params += [ ['RGB_J0710_p591'        ,107.61  ,59.15  ,30 ,90] ]
input_params += [ ['UrsaMinor'             ,227.285 ,67.222 ,30 ,90] ]
input_params += [ ['UrsaMajorII'           ,132.875 ,63.13  ,30 ,90] ]
input_params += [ ['CrabNebula_elev_80_90' ,83.633  ,22.014 ,80 ,90] ]
input_params += [ ['CrabNebula_elev_70_80' ,83.633  ,22.014 ,70 ,80] ]
input_params += [ ['CrabNebula_elev_60_70' ,83.633  ,22.014 ,60 ,70] ]
input_params += [ ['CrabNebula_elev_50_60' ,83.633  ,22.014 ,50 ,60] ]
input_params += [ ['CrabNebula_elev_40_50' ,83.633  ,22.014 ,40 ,50] ]
input_params += [ ['CrabNebula_elev_30_40' ,83.633  ,22.014 ,30 ,40] ]
input_params += [ ['CrabNebula_1p0wobble' ,83.633  ,22.014 ,30 ,90] ]
input_params += [ ['CrabNebula_1p5wobble' ,83.633  ,22.014 ,30 ,90] ]

input_params += [ ['PSR_J1907_p0602'       ,286.975 ,6.337  ,30 ,90] ]
input_params += [ ['PSR_J2021_p4026'       ,305.37  ,40.45  ,30 ,90] ] # gamma cygni
input_params += [ ['SNR_G189_p03'          ,94.213  ,22.503 ,30 ,90] ] # ic 443
input_params += [ ['Geminga'               ,98.476  ,17.770 ,30 ,90] ]
input_params += [ ['PSR_J2021_p3651'       ,305.27  ,36.85  ,30 ,90] ] # Dragonfly

input_params += [ ['LHAASO_J1956_p2845', 299.05  , 28.75  ,30 ,90] ]

#input_params += [ ['Sky_RA55Dec53'         ,55.34   ,52.97  ,30 ,90] ]
#input_params += [ ['SgrA'                  ,266.41  ,-29.00 ,20 ,90] ] # Galactic center

# galactic observations with > 50 runs
#input_params += [ ['SNR_G067p6_p0p9', 299.44, 30.88  ,30 ,90] ]
#input_params += [ ['2HWC_J1953_p294', 298.26  , 29.48  ,30 ,90] ]
#input_params += [ ['Cisne_HS_2013', 304.08  , 36.13  ,30 ,90] ]
#input_params += [ ['CTB109', 345.28  , 58.88  ,30 ,90] ]  # SNR CTB 109 at 3.1 kpc
#input_params += [ ['G079_p00', 308.12  , 40.33  ,30 ,90] ]
#input_params += [ ['LHAASO_J0341_p5258', 55.34  , 52.97  ,30 ,90] ]
#input_params += [ ['LHAASO_J2108_p5157', 317.15  , 51.95  ,30 ,90] ]
#input_params += [ ['HESS_J0632_p057', 98.24  , 5.81  ,30 ,90] ]
#input_params += [ ['LS_V_p4417', 70.25  , 44.53  ,30 ,90] ]  # High Mass X-ray Binary
#input_params += [ ['LS_I_p61_303', 40.13  , 61.23  ,30 ,90] ] # microquasar, a binary system
#input_params += [ ['PSR_B1937_HESS_J1943', 295.45  , 21.44  ,30 ,90] ]
#input_params += [ ['PSR_J2238_p5903', 339.50  , 59.05  ,30 ,90] ]
#input_params += [ ['TeV2032_Cyg_X3', 308.  , 41.23  ,30 ,90] ] # high-mass X-ray binary
#input_params += [ ['Tycho', 6.28  , 64.17  ,30 ,90] ]

#input_params += [ ['Galactic_Plane', 0., 0.  ,30 ,90] ]

input_params += [ ['AUX_files', 0., 0.  ,30 ,90] ]

job_counts = 0
for s in range(0,len(input_params)):
    job_counts += 1
    file = open("run/vts_db_%s.sh"%(input_params[s][0]),"w") 
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
    file.write('#$ -m e\n')
    file.write('. /u/local/Modules/default/init/modules.sh\n')
    file.write('cd %s\n'%(SMI_DIR))
    file.write('python3 veritas_db_query.py "%s" %s %s %s %s "V5"\n'%(input_params[s][0],input_params[s][1],input_params[s][2],input_params[s][3],input_params[s][4]))
    file.write('python3 veritas_db_query.py "%s" %s %s %s %s "V6"\n'%(input_params[s][0],input_params[s][1],input_params[s][2],input_params[s][3],input_params[s][4]))
    file.close() 

job_counts = 0
qfile = open("run/qsub_vts_db.sh","w") 
for s in range(0,len(input_params)):
    job_counts += 1
    qfile.write('qsub -V -N job_vts_%s vts_db_%s.sh\n'%(input_params[s][0],input_params[s][0]))
    qfile.write('sleep 10s\n')
qfile.close() 


