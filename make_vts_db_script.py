
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
input_params += [ ['H1426'                 ,38.222  ,20.273 ,30 ,90] ]
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
input_params += [ ['Sky_RA55Dec53'         ,55.34   ,52.97  ,30 ,90] ]

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
    file.write('python3 veritas_db_query.py "%s" %s %s %s %s\n'%(input_params[s][0],input_params[s][1],input_params[s][2],input_params[s][3],input_params[s][4]))
    file.close() 

job_counts = 0
qfile = open("run/qsub_vts_db.sh","w") 
for s in range(0,len(input_params)):
    job_counts += 1
    qfile.write('qsub -V -N job_vts_%s vts_db_%s.sh\n'%(input_params[s][0],input_params[s][0]))
    qfile.write('sleep 10s\n')
qfile.close() 


