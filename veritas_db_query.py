
import os.path
import sys,ROOT
from ROOT import *
import sys
import pymysql
from datetime import datetime, timedelta
import math
from astropy.coordinates import SkyCoord
from astropy import units as my_unit
import numpy as np

# https://veritas.sao.arizona.edu/wiki/VOFFLINE_Database_Tables
# https://veritas.sao.arizona.edu/wiki/Ryan_Dickherber%27s_Wiki
# Log Gen script: http://veritash.sao.arizona.edu:8081/OfflineAnalysis-WG/230319_221625/query_night

all_runs_info = []

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

def ReadHAWCTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    source_flux = []
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
        if 'flux:' in line:
            target_flux = line.lstrip('          flux: ')
        if 'index systematic uncertainty down:' in line:
            source_name += [target_name]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            source_flux += [float(target_flux)/2.34204e-13]
            target_name = ''
            target_ra = ''
            target_dec = ''
            target_flux = ''
    return source_name, source_ra, source_dec, source_flux

def ConvertRaDecToGalactic(ra, dec):
    my_sky = SkyCoord(ra*my_unit.deg, dec*my_unit.deg, frame='icrs')
    return my_sky.galactic.l.deg, my_sky.galactic.b.deg

def get_sec(time_str):
    """Get Seconds from time."""
    if time_str=='None': return 0
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + int(s)

def print_all_runs_nsb():

    print ('Connect to DB...')
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    
    print ('Read tblRun_Info...')
    all_runs_src = {}
    all_runs_weather = {}
    all_runs_type = {}
    query = 'SELECT run_id,source_id,weather,run_type FROM tblRun_Info'
    crs.execute(query)
    res_run_info = crs.fetchall()
    for x in res_run_info:
        all_runs_src[x['run_id']] = x['source_id']
        all_runs_weather[x['run_id']] = x['weather']
        all_runs_type[x['run_id']] = x['run_type']

    print ('Read tblProcess_RunTime...')
    query = 'SELECT run_id,db_start_time,db_end_time FROM tblProcess_RunTime'
    crs.execute(query)
    res_runtime = crs.fetchall()

    print ('Calculate NSB...')
    with open('/gamma_raid/userspace/rshang/SMI_AUX/NSB_allruns.txt', 'w') as file:
        previous_run_id = 0
        for x in res_runtime:

            if x['run_id']==previous_run_id: continue
            previous_run_id = x['run_id']

            if x['run_id']<46642: continue
            if all_runs_weather[x['run_id']]==None: continue
            if 'C' in all_runs_weather[x['run_id']]: continue
            if 'D' in all_runs_weather[x['run_id']]: continue
            if 'F' in all_runs_weather[x['run_id']]: continue

            if all_runs_src[x['run_id']]==None: continue
            if all_runs_src[x['run_id']]=='engineering': continue
            if all_runs_src[x['run_id']]=='other': continue
            if all_runs_src[x['run_id']]=='none': continue
            if all_runs_src[x['run_id']]=='qi': continue
            if all_runs_src[x['run_id']]=='NOSOURCE': continue

            if all_runs_type[x['run_id']]=='flasher': continue

            run_start_time = x['db_start_time']
            run_end_time = x['db_end_time']
            current_avg_run = 0.
            total_channels = 0.
            for tel in range(0,1):
                for ch in range(1,10):
                    current_avg_ch = 0.
                    total_entries = 0.
                    query = "SELECT current_meas FROM tblHV_Telescope%s_Status WHERE db_start_time>'%s' AND db_start_time<'%s' AND channel=%s"%(tel,run_start_time,run_end_time,ch*40)
                    crs.execute(query)
                    res2 = crs.fetchall()
                    for y in res2:
                        current_avg_ch += y['current_meas']
                        total_entries += 1.
                    if total_entries==0.: continue
                    current_avg_ch = current_avg_ch/total_entries
                    current_avg_run += current_avg_ch
                    total_channels += 1.
            if total_channels==0.: continue
            current_avg_run = current_avg_run/total_channels
            print ('%s %0.2f'%(x['run_id'],current_avg_run))
            file.write('%s %0.2f\n'%(x['run_id'],current_avg_run))


def get_run_nsb(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    
    query = 'SELECT run_id,db_start_time,db_end_time FROM tblProcess_RunTime WHERE run_id=%s'%(run_id)
    crs.execute(query)
    res = crs.fetchall()
    run_start_time = res[0]['db_start_time']
    run_end_time = res[0]['db_end_time']
    current_avg_run = 0.
    total_channels = 0.
    for ch in range(1,500):
        current_avg_ch = 0.
        total_entries = 0.
        query = "SELECT current_meas FROM tblHV_Telescope1_Status WHERE db_start_time>'%s' AND db_start_time<'%s' AND channel=%s"%(run_start_time,run_end_time,ch)
        crs.execute(query)
        res = crs.fetchall()
        for x in res:
            current_avg_ch += x['current_meas']
            total_entries += 1.
        current_avg_ch = current_avg_ch/total_entries
        current_avg_run += current_avg_ch
        total_channels += 1.
    current_avg_run = current_avg_run/total_channels
    print ('run_id = %s, current_avg_run = %s'%(run_id,current_avg_run))

def get_all_runs_info(epoch,obs_type):

    global all_runs_info

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    all_runs_comments = {}
    query = 'SELECT run_id,usable_duration FROM tblRun_Analysis_Comments'
    crs.execute(query)
    # fetch from cursor
    res_comment = crs.fetchall()
    for x in res_comment:
        all_runs_comments[x['run_id']] = get_sec(str(x['usable_duration']))

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    
    query = 'SELECT run_id,source_id,data_start_time,data_end_time,weather,run_type FROM tblRun_Info'
    crs.execute(query)
    res_run_info = crs.fetchall()

    all_runs_info = []

    for x in res_run_info:

        if x['run_type']!=obs_type: continue
        if x['weather']==None: continue
        if 'C' in x['weather']: continue
        if 'D' in x['weather']: continue
        if 'F' in x['weather']: continue

        if x['run_id']<46642: continue
        if x['run_id']<46642:
            if not epoch=='V4': continue
        if x['run_id']>=46642 and x['run_id']<63373:
            if not epoch=='V5': continue
        if x['run_id']>=63373:
            if not epoch=='V6': continue

        if x['source_id']==None: continue
        if x['source_id']=='engineering': continue
        if x['source_id']=='other': continue
        if x['source_id']=='none': continue
        if x['source_id']=='qi': continue
        if x['source_id']=='NOSOURCE': continue

        if x['run_id'] in all_runs_comments:
            run_usable_time = all_runs_comments[x['run_id']]
            if run_usable_time<5.*60.: continue
        else:
            continue

        el_avg_run = 0.
        az_avg_run = 0.
        total_entries = 0.

        timestamp_start = '%s'%(x['data_start_time'])
        timestamp_end = '%s'%(x['data_end_time'])
        if timestamp_start=='None': continue
        if timestamp_end=='None': continue
        timestamp_start = timestamp_start.replace(' ','').replace('-','').replace(':','')
        timestamp_end = timestamp_end.replace(' ','').replace('-','').replace(':','')
        timestamp_start += '000'
        timestamp_end += '000'
        query = "SELECT elevation_target,azimuth_target FROM tblPositioner_Telescope1_Status WHERE timestamp>%s AND timestamp<%s"%(timestamp_start,timestamp_end)
        crs.execute(query)
        res2 = crs.fetchall()
        if len(res2)==0: continue
        el_avg_run += res2[0]['elevation_target']
        el_avg_run += res2[len(res2)-1]['elevation_target']
        az_avg_run += res2[0]['azimuth_target']
        az_avg_run += res2[len(res2)-1]['azimuth_target']
        total_entries = 2.
        el_avg_run = el_avg_run/total_entries*180./math.pi
        az_avg_run = az_avg_run/total_entries*180./math.pi

        all_runs_info += [[x['run_id'],x['source_id'],el_avg_run,az_avg_run,run_usable_time]]


def print_all_runs_l3rate():

    print ('Connect to DB...')
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    
    print ('Read tblRun_Info...')
    all_runs_src = {}
    all_runs_weather = {}
    all_runs_type = {}
    query = 'SELECT run_id,source_id,weather,run_type,data_start_time,data_end_time FROM tblRun_Info'
    crs.execute(query)
    res_run_info = crs.fetchall()
    for x in res_run_info:
        all_runs_src[x['run_id']] = x['source_id']
        all_runs_weather[x['run_id']] = x['weather']
        all_runs_type[x['run_id']] = x['run_type']

    print('Calculate Run L3 rate...')
    with open('/gamma_raid/userspace/rshang/SMI_AUX/l3rate_allruns.txt', 'w') as file:
        previous_run_id = 0
        for x in res_run_info:

            if x['run_id']==previous_run_id: continue
            previous_run_id = x['run_id']

            if x['run_id']<46642: continue
            if all_runs_weather[x['run_id']]==None: continue
            if 'C' in all_runs_weather[x['run_id']]: continue
            if 'D' in all_runs_weather[x['run_id']]: continue
            if 'F' in all_runs_weather[x['run_id']]: continue

            if all_runs_src[x['run_id']]==None: continue
            if all_runs_src[x['run_id']]=='engineering': continue
            if all_runs_src[x['run_id']]=='other': continue
            if all_runs_src[x['run_id']]=='none': continue
            if all_runs_src[x['run_id']]=='qi': continue
            if all_runs_src[x['run_id']]=='NOSOURCE': continue

            if all_runs_type[x['run_id']]=='flasher': continue

            timestamp_start = '%s'%(x['data_start_time'])
            timestamp_end = '%s'%(x['data_end_time'])
            if timestamp_start=='None': continue
            if timestamp_end=='None': continue
            timestamp_start = timestamp_start.replace(' ','').replace('-','').replace(':','')
            timestamp_end = timestamp_end.replace(' ','').replace('-','').replace(':','')
            timestamp_start += '000'
            timestamp_end += '000'
            l3_avg_run = 0.
            total_entries = 0.
            query = "SELECT L3 FROM tblL3_Array_TriggerInfo WHERE timestamp>%s AND timestamp<%s AND run_id=%s"%(timestamp_start,timestamp_end,x['run_id'])
            crs.execute(query)
            res_pos = crs.fetchall()
            for y in res_pos:
                l3_avg_run += float(y['L3'])
                total_entries += 1.
            if total_entries==0.: continue
            l3_avg_run = l3_avg_run/total_entries
            print ('%s %0.1f'%(x['run_id'],l3_avg_run))
            file.write('%s %0.1f\n'%(x['run_id'],l3_avg_run))

def print_all_runs_el_az():

    print ('Connect to DB...')
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    
    print ('Read tblRun_Info...')
    all_runs_src = {}
    all_runs_weather = {}
    all_runs_type = {}
    query = 'SELECT run_id,source_id,weather,run_type,data_start_time,data_end_time FROM tblRun_Info'
    crs.execute(query)
    res_run_info = crs.fetchall()
    for x in res_run_info:
        all_runs_src[x['run_id']] = x['source_id']
        all_runs_weather[x['run_id']] = x['weather']
        all_runs_type[x['run_id']] = x['run_type']

    print('Calculate Run Elev Azim...')
    with open('/gamma_raid/userspace/rshang/SMI_AUX/elaz_allruns.txt', 'w') as file:
        previous_run_id = 0
        for x in res_run_info:

            if x['run_id']==previous_run_id: continue
            previous_run_id = x['run_id']

            if x['run_id']<46642: continue
            if all_runs_weather[x['run_id']]==None: continue
            if 'C' in all_runs_weather[x['run_id']]: continue
            if 'D' in all_runs_weather[x['run_id']]: continue
            if 'F' in all_runs_weather[x['run_id']]: continue

            if all_runs_src[x['run_id']]==None: continue
            if all_runs_src[x['run_id']]=='engineering': continue
            if all_runs_src[x['run_id']]=='other': continue
            if all_runs_src[x['run_id']]=='none': continue
            if all_runs_src[x['run_id']]=='qi': continue
            if all_runs_src[x['run_id']]=='NOSOURCE': continue

            if all_runs_type[x['run_id']]=='flasher': continue

            timestamp_start = '%s'%(x['data_start_time'])
            timestamp_end = '%s'%(x['data_end_time'])
            if timestamp_start=='None': continue
            if timestamp_end=='None': continue
            timestamp_start = timestamp_start.replace(' ','').replace('-','').replace(':','')
            timestamp_end = timestamp_end.replace(' ','').replace('-','').replace(':','')
            timestamp_start += '000'
            timestamp_end += '000'
            el_avg_run = 0.
            az_avg_run = 0.
            total_entries = 0.
            query = "SELECT elevation_target,azimuth_target FROM tblPositioner_Telescope1_Status WHERE timestamp>%s AND timestamp<%s"%(timestamp_start,timestamp_end)
            crs.execute(query)
            res_pos = crs.fetchall()
            for y in res_pos:
                el_avg_run += y['elevation_target']
                az_avg_run += y['azimuth_target']
                total_entries += 1.
            if total_entries==0.: continue
            el_avg_run = el_avg_run/total_entries*180./math.pi
            az_avg_run = az_avg_run/total_entries*180./math.pi
            print ('%s %0.2f %0.2f'%(x['run_id'],el_avg_run,az_avg_run))
            file.write('%s %0.2f %0.2f\n'%(x['run_id'],el_avg_run,az_avg_run))

def get_run_el_az(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    
    query = 'SELECT run_id,data_start_time,data_end_time FROM tblRun_Info WHERE run_id=%s'%(run_id)
    crs.execute(query)
    res = crs.fetchall()
    timestamp_start = '%s'%(res[0]['data_start_time'])
    timestamp_end = '%s'%(res[0]['data_end_time'])
    if timestamp_start=='None': return 0, 0
    if timestamp_end=='None': return 0, 0
    timestamp_start = timestamp_start.replace(' ','').replace('-','').replace(':','')
    timestamp_end = timestamp_end.replace(' ','').replace('-','').replace(':','')
    timestamp_start += '000'
    timestamp_end += '000'
    #print ('timestamp_start = %s'%(timestamp_start))
    el_avg_run = 0.
    az_avg_run = 0.
    total_entries = 0.
    query = "SELECT elevation_target,azimuth_target FROM tblPositioner_Telescope1_Status WHERE timestamp>%s AND timestamp<%s"%(timestamp_start,timestamp_end)
    crs.execute(query)
    res = crs.fetchall()
    for x in res:
        el_avg_run += x['elevation_target']
        az_avg_run += x['azimuth_target']
        total_entries += 1.
    if total_entries==0.:
        return 0., 0.
    el_avg_run = el_avg_run/total_entries*180./math.pi
    az_avg_run = az_avg_run/total_entries*180./math.pi
    #print ('run_id = %s, el_avg_run = %s, az_avg_run = %s'%(run_id,el_avg_run,az_avg_run))
    return el_avg_run, az_avg_run

def get_run_category(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT run_id,data_category FROM tblRun_Analysis_Comments WHERE run_id=%s'%(run_id)
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    print(res[0]['run_id'],res[0]['data_category'])

def print_all_runs_timecut():

    print ('Connect to DB...')
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    print ('Read tblRun_Analysis_Comments...')
    query = 'SELECT run_id,time_cut_mask FROM tblRun_Analysis_Comments'
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()

    print ('Get timecut...')
    with open('/gamma_raid/userspace/rshang/SMI_AUX/timecuts_allruns.txt', 'w') as file:
        for x in res:
            print('%s %s'%(x['run_id'],x['time_cut_mask']))
            file.write('%s %s\n'%(x['run_id'],x['time_cut_mask']))

def get_run_timecut(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT run_id,time_cut_mask FROM tblRun_Analysis_Comments WHERE run_id=%s'%(run_id)
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    print('run_id = %s, time_cut_mask = %s'%(res[0]['run_id'],res[0]['time_cut_mask']))

def print_all_runs_usable_duration():

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT run_id,usable_duration FROM tblRun_Analysis_Comments'
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()

    with open('/gamma_raid/userspace/rshang/SMI_AUX/usable_time_allruns.txt', 'w') as file:
        for x in res:
            duration = get_sec(str(x['usable_duration']))
            print ('%s %s'%(x['run_id'],duration))
            file.write('%s %s\n'%(x['run_id'],duration))

def get_run_usable_duration(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT run_id,usable_duration FROM tblRun_Analysis_Comments WHERE run_id=%s'%(run_id)
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    #print('run_id = %s, usable_duration = %s'%(res[0]['run_id'],get_sec(str(res[0]['usable_duration']))))
    if len(res)==0:
        return 0.
    return get_sec(str(res[0]['usable_duration']))

def get_run_ra_dec(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT run_id,source_id,offsetRA,offsetDEC,offset_angle FROM tblRun_Info WHERE run_id=%s'%(run_id)
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    source_name = res[0]['source_id']
    query2 = "SELECT source_id,ra,decl FROM tblObserving_Sources WHERE source_id='%s'"%(source_name)
    crs.execute(query2)
    res2 = crs.fetchall()
    source_ra = res2[0]['ra']
    source_dec = res2[0]['decl']
    run_ra = (source_ra + res[0]['offsetRA'])*180./math.pi
    run_dec = (source_dec + res[0]['offsetDEC'])*180./math.pi
    #print ('source_name = %s, source_ra = %s, source_dec = %s'%(source_name,source_ra*180./math.pi,source_dec*180./math.pi))
    #print ('run_id = %s, run_ra = %s, run_dec = %s'%(run_id,run_ra,run_dec))
    return run_ra, run_dec

def get_run_weather(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT run_id,weather FROM tblRun_Info WHERE run_id=%s'%(run_id)
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    #print ('run_id = %s, weather = %s'%(run_id,res[0]['weather']))
    return res[0]['weather']

def print_all_runs_type():

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT run_id,run_type FROM tblRun_Info'
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()

    with open('/gamma_raid/userspace/rshang/SMI_AUX/runtype_allruns.txt', 'w') as file:
        for x in res:
            print ('%s %s'%(x['run_id'],x['run_type']))
            file.write('%s %s\n'%(x['run_id'],x['run_type']))

def get_run_type(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT run_id,run_type FROM tblRun_Info WHERE run_id=%s'%(run_id)
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    #print ('run_id = %s, run_type = %s'%(run_id,res[0]['run_type']))
    return res[0]['run_type']

def find_runs_near_galactic_plane(obs_name,epoch,obs_type,gal_b_low,gal_b_up):

    out_file = open('output_vts_hours/%s.txt'%(obs_name),"w")

    list_on_run_ids = []
    list_on_sources = []

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    all_src_ra = {}
    all_src_dec = {}
    runs_per_src = {}
    query = 'SELECT source_id,ra,decl FROM tblObserving_Sources'
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    for x in res:
        source_name = x['source_id']
        #print ('source_name = %s'%(source_name))
        source_ra = x['ra']*180./math.pi
        source_dec = x['decl']*180./math.pi
        all_src_ra[x['source_id']] = source_ra
        all_src_dec[x['source_id']] = source_dec
        runs_per_src[x['source_id']] = 0
        source_gal_l, source_gal_b = ConvertRaDecToGalactic(source_ra,source_dec)
        if abs(source_gal_b)>gal_b_low and abs(source_gal_b)<gal_b_up:
            list_on_sources += [source_name]
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++')
    for src in range(0,len(list_on_sources)):
        print (list_on_sources[src])

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    all_runs_duration = {}
    query = 'SELECT run_id,usable_duration FROM tblRun_Analysis_Comments'
    crs.execute(query)
    # fetch from cursor
    res_comment = crs.fetchall()
    for x in res_comment:
        all_runs_duration[x['run_id']] = get_sec(str(x['usable_duration']))

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    print ('Read tblRun_Info...')
    all_runs_src = {}
    all_runs_weather = {}
    all_runs_type = {}
    query = 'SELECT run_id,source_id,weather,run_type FROM tblRun_Info'
    crs.execute(query)
    res_run_info = crs.fetchall()
    for x in res_run_info:
        all_runs_src[x['run_id']] = x['source_id']
        all_runs_weather[x['run_id']] = x['weather']
        all_runs_type[x['run_id']] = x['run_type']

    query = 'SELECT run_id,source_id FROM tblRun_Info'
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    for x in res:

        if x['run_id']<46642: continue
        if x['run_id']<46642:
            if not epoch=='V4': continue
        if x['run_id']>=46642 and x['run_id']<63373:
            if not epoch=='V5': continue
        if x['run_id']>=63373:
            if not epoch=='V6': continue

        if x['source_id']==None: continue
        if x['source_id']=='engineering': continue
        if x['source_id']=='other': continue
        if x['source_id']=='none': continue
        if x['source_id']=='qi': continue
        if x['source_id']=='NOSOURCE': continue

        source_name = x['source_id']
        is_good_src = False
        for src in range(0,len(list_on_sources)):
            if source_name==list_on_sources[src]:
                is_good_src = True
        if not is_good_src: continue

        if x['run_id'] in all_runs_type:
            run_type = all_runs_type[x['run_id']]
            if run_type!=obs_type: continue
        else:
            continue

        if x['run_id'] in all_runs_weather:
            run_weather = all_runs_weather[x['run_id']]
            if run_weather==None: continue
            if 'C' in run_weather: continue
            if 'D' in run_weather: continue
            if 'F' in run_weather: continue
        else:
            continue

        if x['run_id'] in all_runs_duration:
            run_duration = all_runs_duration[x['run_id']]
            if run_duration<5.*60.: continue
        else:
            continue

        on_run_el, on_run_az = get_run_el_az(x['run_id'])
        if on_run_el<45.: continue

        print ('run_id = %s, source_name = %s, RA = %0.2f, Dec = %0.2f'%(x['run_id'],x['source_id'],all_src_ra[x['source_id']],all_src_dec[x['source_id']]))
        out_file.write('run_id = %s, source_name = %s, RA = %0.2f, Dec = %0.2f \n'%(x['run_id'],x['source_id'],all_src_ra[x['source_id']],all_src_dec[x['source_id']]))
        list_on_run_ids += [x['run_id']]
        runs_per_src[x['source_id']] += 1

    out_file.write('++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    out_file.write('Source list\n')
    for src in range(0,len(list_on_sources)):
        src_name = list_on_sources[src]
        #if runs_per_src[src_name]<50: continue
        if runs_per_src[src_name]<10: continue
        print ('%s, runs = %s, RA = %0.2f, Dec = %0.2f'%(src_name,runs_per_src[src_name],all_src_ra[src_name],all_src_dec[src_name]))
        out_file.write('%s, runs = %s, RA = %0.2f, Dec = %0.2f \n'%(src_name,runs_per_src[src_name],all_src_ra[src_name],all_src_dec[src_name]))

    out_file.close()

def find_on_runs_around_source(obs_name,obs_ra,obs_dec,epoch,obs_type,elev_range,search_radius):

    global all_runs_info

    out_file = open('output_vts_hours/%s.txt'%(obs_name),"w")
    out_run_file = open('output_vts_hours/RunList_%s.txt'%(obs_name),"w")

    list_on_run_ids = []
    list_on_run_elev = []
    list_on_run_azim = []
    list_on_sources = []

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT source_id,ra,decl FROM tblObserving_Sources'
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    for x in res:
        source_name = x['source_id']
        #print ('source_name = %s'%(source_name))
        source_ra = x['ra']*180./math.pi
        source_dec = x['decl']*180./math.pi
        source_gal_l, source_gal_b = ConvertRaDecToGalactic(source_ra,source_dec)
        distance = pow(pow(obs_ra-source_ra,2)+pow(obs_dec-source_dec,2),0.5)
        if distance<search_radius:
            list_on_sources += [source_name]
        #if distance>10.:
        #    if 'HWC' in source_name: continue
        #    #if abs(source_ra-83.633)<2. and abs(source_dec-22.014)<2.: continue # Crab
        #    #if abs(source_ra-166.079)<2. and abs(source_dec-38.195)<2.: continue # Mrk 421
        #    #if abs(source_ra-253.467)<2. and abs(source_dec-39.76)<2.: continue # Mrk 501
        #    if abs(source_ra-98.117)<3. and abs(source_dec-17.367)<3.: continue # Geminga
        #    if abs(source_gal_b)>10.:
        #        list_off_sources += [source_name]
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++')
    for src in range(0,len(list_on_sources)):
        print (list_on_sources[src])

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    all_runs_duration = {}
    query = 'SELECT run_id,usable_duration FROM tblRun_Analysis_Comments'
    crs.execute(query)
    # fetch from cursor
    res_comment = crs.fetchall()
    for x in res_comment:
        all_runs_duration[x['run_id']] = get_sec(str(x['usable_duration']))

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    print ('Read tblRun_Info...')
    all_runs_src = {}
    all_runs_weather = {}
    all_runs_type = {}
    query = 'SELECT run_id,source_id,weather,run_type FROM tblRun_Info'
    crs.execute(query)
    res_run_info = crs.fetchall()
    for x in res_run_info:
        all_runs_src[x['run_id']] = x['source_id']
        all_runs_weather[x['run_id']] = x['weather']
        all_runs_type[x['run_id']] = x['run_type']

    query = 'SELECT run_id,source_id FROM tblRun_Info'
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    for x in res:

        if x['run_id']<46642: continue
        if x['run_id']<46642:
            if not epoch=='V4': continue
        if x['run_id']>=46642 and x['run_id']<63373:
            if not epoch=='V5': continue
        if x['run_id']>=63373:
            if not epoch=='V6': continue

        if x['source_id']==None: continue
        if x['source_id']=='engineering': continue
        if x['source_id']=='other': continue
        if x['source_id']=='none': continue
        if x['source_id']=='qi': continue
        if x['source_id']=='NOSOURCE': continue

        source_name = x['source_id']
        is_good_src = False
        for src in range(0,len(list_on_sources)):
            if source_name==list_on_sources[src]:
                is_good_src = True
        if not is_good_src: continue

        if x['run_id'] in all_runs_type:
            run_type = all_runs_type[x['run_id']]
            if run_type!=obs_type: continue
        else:
            continue

        if x['run_id'] in all_runs_weather:
            run_weather = all_runs_weather[x['run_id']]
            if run_weather==None: continue
            if 'C' in run_weather: continue
            if 'D' in run_weather: continue
            if 'F' in run_weather: continue
        else:
            continue

        if x['run_id'] in all_runs_duration:
            run_duration = all_runs_duration[x['run_id']]
            if run_duration<5.*60.: continue
        else:
            continue


        on_run_el, on_run_az = get_run_el_az(x['run_id'])
        if on_run_el<elev_range[0]: continue
        if on_run_el>elev_range[1]: continue

        on_run_ra, on_run_dec = get_run_ra_dec(x['run_id'])
        run_offset = pow(pow(obs_ra-on_run_ra,2)+pow(obs_dec-on_run_dec,2),0.5)
        if '1p0wobble' in obs_name:
            if run_offset>1.2 or run_offset<0.6: continue
        if '1p5wobble' in obs_name:
            if run_offset>1.8 or run_offset<1.2: continue

        if use_local_data:
            file_path = '/gamma_raid/userspace/rshang/analysis/Results/v487/%s.anasum.root'%(int(x['run_id']))
            root_folder = "run_%s/stereo/DL3EventTree"%(int(x['run_id']))
            if os.path.exists(file_path):
                InputFile = ROOT.TFile(file_path)
                mytree = InputFile.Get(root_folder)
                if InputFile.IsZombie():
                    #print ('%s, something very wrong, cannot use this file'%(int(x['run_id']))
                    continue
                elif InputFile.TestBit(ROOT.TFile.kRecovered):
                    #print ('%s, the Recover procedure has been run when opening the file'%(int(x['run_id']))
                    continue
                InputFile.Close()
            else:
                continue

        print ('run_id = %s, source_name = %s'%(x['run_id'],x['source_id']))
        out_file.write('run_id = %s, source_name = %s \n'%(x['run_id'],x['source_id']))

        list_on_run_ids += [x['run_id']]
        list_on_run_elev += [on_run_el]
        list_on_run_azim += [on_run_az]

    #list_pairs = zip(list_on_run_ids, list_on_run_elev)
    #sorted_pairs = sorted(list_pairs, key=lambda x: x[1])
    #list_on_run_ids = [x[0] for x in sorted_pairs]

    out_file.write('++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    out_file.write('ON run list\n')
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('ON run list')
    for run in range(0,len(list_on_run_ids)):
        out_file.write('%s\n'%(list_on_run_ids[run]))
        out_run_file.write('%s\n'%(list_on_run_ids[run]))
        print (list_on_run_ids[run])

    out_file.close()

    return list_on_run_ids

def find_off_runs_around_source(obs_name,obs_ra,obs_dec,epoch,obs_type,elev_range,list_on_run_ids,is_imposter,file_name):

    global all_runs_info

    list_on_run_nmatch = []
    list_off_run_ids = []
    list_no_repeat_off_run_ids = []
    list_off_sources = []

    out_pair_file = open('output_vts_hours/%s_%s.txt'%(file_name,obs_name),"w")
    out_off_file = open('output_vts_hours/%s_OFFRuns_%s.txt'%(file_name,obs_name),"w")

    require_nmatch = 5
    if not is_imposter:
        require_nmatch = 3
        #require_nmatch = 16

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()

    query = 'SELECT source_id,ra,decl FROM tblObserving_Sources'
    crs.execute(query)
    # fetch from cursor
    res = crs.fetchall()
    for x in res:
        source_name = x['source_id']
        #print ('source_name = %s'%(source_name))
        source_ra = x['ra']*180./math.pi
        source_dec = x['decl']*180./math.pi
        source_gal_l, source_gal_b = ConvertRaDecToGalactic(source_ra,source_dec)
        distance = pow(pow(obs_ra-source_ra,2)+pow(obs_dec-source_dec,2),0.5)
        if distance>10.:
            if 'HWC' in source_name: continue
            #if abs(source_ra-83.633)<2. and abs(source_dec-22.014)<2.: continue # Crab
            #if abs(source_ra-166.079)<2. and abs(source_dec-38.195)<2.: continue # Mrk 421
            #if abs(source_ra-253.467)<2. and abs(source_dec-39.76)<2.: continue # Mrk 501
            if abs(source_ra-98.117)<3. and abs(source_dec-17.367)<3.: continue # Geminga
            if abs(source_gal_b)>10.:
                list_off_sources += [source_name]
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++')
    for src in range(0,len(list_off_sources)):
        print (list_off_sources[src])


    for on_run in range(0,len(list_on_run_ids)):

        on_run_el, on_run_az = get_run_el_az(list_on_run_ids[on_run])
        number_off_runs = 0

        sum_on_run_elev = 0.
        sum_off_run_elev = 0.
        for run in range(0,len(all_runs_info)):

            if number_off_runs>=require_nmatch: continue

            if all_runs_info[run][0]==list_on_run_ids[on_run]: continue

            if use_local_data:
                file_path = '/gamma_raid/userspace/rshang/analysis/Results/v487/%s.anasum.root'%(all_runs_info[run][0])
                root_folder = "run_%s/stereo/DL3EventTree"%(all_runs_info[run][0])
                if os.path.exists(file_path):
                    InputFile = ROOT.TFile(file_path)
                    mytree = InputFile.Get(root_folder)
                    if InputFile.IsZombie():
                        print ('%s, something very wrong, cannot use this file'%(all_runs_info[run][0]))
                        continue
                    elif InputFile.TestBit(ROOT.TFile.kRecovered):
                        print ('%s, the Recover procedure has been run when opening the file'%(all_runs_info[run][0]))
                        continue
                    InputFile.Close()
                else:
                    continue

            #if abs(all_runs_info[run][0]-list_on_run_ids[on_run])>40000: continue
            already_used = False
            if is_imposter:
                for off_run in range(0,len(list_off_run_ids)):
                    if all_runs_info[run][0]==list_off_run_ids[off_run][1]: already_used = True
                if already_used: continue

            if all_runs_info[run][0]<46642: continue
            if all_runs_info[run][0]<46642:
                if not epoch=='V4': continue
            if all_runs_info[run][0]>=46642 and all_runs_info[run][0]<63373:
                if not epoch=='V5': continue
            if all_runs_info[run][0]>=63373:
                if not epoch=='V6': continue

            source_name = all_runs_info[run][1]
            is_good_src = False
            for src in range(0,len(list_off_sources)):
                if source_name==list_off_sources[src]:
                    is_good_src = True
            if not is_good_src: continue

            off_run_duration = all_runs_info[run][4]
            if off_run_duration<20.*60.: continue

            off_run_el = all_runs_info[run][2]
            off_run_az = all_runs_info[run][3]
            if on_run_el==0.: continue
            if off_run_el==0.: continue
            
            #if off_run_az>180.: off_run_az = 360.-off_run_az 
            #if on_run_az>180.: on_run_az = 360.-on_run_az 
            #delta_azim = off_run_az-on_run_az
            delta_azim = abs(off_run_az-on_run_az)
            if delta_azim>180.: delta_azim = 360.-delta_azim

            delta_airmass = (1./math.sin(on_run_el*math.pi/180.)-1./math.sin(off_run_el*math.pi/180.));
            delta_elev = off_run_el-on_run_el

            if is_imposter:
                if abs(delta_elev)>5.: continue
            else:
                if (sum_off_run_elev-sum_on_run_elev)>0.:
                    if (off_run_el-on_run_el)>0.: continue
                else:
                    if (off_run_el-on_run_el)<0.: continue
                #if abs(delta_airmass)>0.1: continue
                if abs(delta_elev)>5.: continue
                if abs(delta_azim)>10.: continue

            list_off_run_ids += [[int(list_on_run_ids[on_run]),int(all_runs_info[run][0]),on_run_el,off_run_el]]
            number_off_runs += 1
            sum_on_run_elev += on_run_el
            sum_off_run_elev += off_run_el

            already_used = False
            for off_run in range(0,len(list_no_repeat_off_run_ids)):
                if all_runs_info[run][0]==list_no_repeat_off_run_ids[off_run]: already_used = True
            if not already_used:
                list_no_repeat_off_run_ids += [all_runs_info[run][0]]

        list_on_run_nmatch += [number_off_runs]

    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('%s list'%(file_name))
    for run in range(0,len(list_off_run_ids)):
        print ('%s %s'%(list_off_run_ids[run][0],list_off_run_ids[run][1]))
        out_pair_file.write('%s %s\n'%(list_off_run_ids[run][0],list_off_run_ids[run][1]))
    out_pair_file.close()

    out_file = open('output_vts_hours/%s.txt'%(obs_name),"a")
    out_file.write('++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    out_file.write('%s ON runs do not have enought matches\n'%(file_name))
    for run in range(0,len(list_on_run_ids)):
        if list_on_run_nmatch[run]<require_nmatch:
            out_file.write('%s has %s OFF runs.\n'%(list_on_run_ids[run],list_on_run_nmatch[run]))


    out_file.write('++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    out_file.write('ON run list\n')
    out_file.write('%s OFF/ON ratio = %0.2f\n'%( file_name, float(len(list_off_run_ids))/float(len(list_on_run_ids)) ))
    out_file.close()

    for run in range(0,len(list_no_repeat_off_run_ids)):
        out_off_file.write('%s\n'%(list_no_repeat_off_run_ids[run]))
    out_off_file.close()
    
    list_off_run_ids = np.array(list_off_run_ids)
    return list_off_run_ids[:,1]

use_local_data = False
#run_obs_type = 'obsLowHV' # RHV
run_obs_type = 'observing'
find_off = False

#target_hwc_name, target_hwc_ra, target_hwc_dec, target_hwc_flux = ReadHAWCTargetListFromFile('Cat_3HWC.txt')
#for hwc in range(3,len(target_hwc_name)):
#    print ('Search VTS data for %s'%(target_hwc_name[hwc]))
#    obs_ra = target_hwc_ra[hwc]
#    obs_dec = target_hwc_dec[hwc]
#    obs_name = target_hwc_name[hwc].replace(' ','_').replace('+','_p').replace('-','_m')
#    obs_name += '_%s'%(run_epoch)
#    find_runs_around_source(obs_name,obs_ra,obs_dec,run_epoch,run_obs_type,find_off)

input_name = sys.argv[1]
input_ra = float(sys.argv[2])
input_dec = float(sys.argv[3])
input_elev_low = float(sys.argv[4])
input_elev_up = float(sys.argv[5])
input_epoch = sys.argv[6]

run_epoch = input_epoch # 'V5' or 'V6'

if input_name=='AUX_files':

    print_all_runs_usable_duration()
    print_all_runs_type()
    print_all_runs_timecut()
    print_all_runs_el_az()
    print_all_runs_nsb()
    #print_all_runs_l3rate()  # do not use

elif input_name=='Galactic_Plane':

    obs_name = 'Galactic_Plane_%s'%(run_epoch)
    find_runs_near_galactic_plane(obs_name,run_epoch,run_obs_type,0.0,5.0)

elif input_name=='LHAASO_Catalog':

    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('Get all runs El Az...')
    get_all_runs_info(run_epoch,run_obs_type)
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++')

    target_lhs_name, target_lhs_ra, target_lhs_dec = ReadLhaasoListFromFile()
    for lhs in range(0,len(target_lhs_name)):
        print ('Search VTS data for %s'%(target_lhs_name[lhs]))
        obs_ra = target_lhs_ra[lhs]
        obs_dec = target_lhs_dec[lhs]
        obs_name = target_lhs_name[lhs].replace(' ','_').replace('+','_p').replace('-','_m')
        obs_name += '_%s'%(run_epoch)
        run_elev_range = [input_elev_low,input_elev_up]
        my_list_on_run_ids = find_on_runs_around_source(obs_name,obs_ra,obs_dec,run_epoch,run_obs_type,run_elev_range,2.0)

else:

    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('Get all runs El Az...')
    get_all_runs_info(run_epoch,run_obs_type)
    print ('++++++++++++++++++++++++++++++++++++++++++++++++++++')
    
    search_radius = 2.0
    job_tag = '2p0deg'

    find_off = True
    find_imposter = True
    obs_name = '%s_%s_%s'%(input_name,job_tag,run_epoch)
    obs_ra = input_ra
    obs_dec = input_dec
    run_elev_range = [input_elev_low,input_elev_up]
    
    use_local_data = False
    my_list_on_run_ids = find_on_runs_around_source(obs_name,obs_ra,obs_dec,run_epoch,run_obs_type,run_elev_range,search_radius)
    use_local_data = True
    my_list_off_run_ids = find_off_runs_around_source(obs_name,obs_ra,obs_dec,run_epoch,run_obs_type,run_elev_range,my_list_on_run_ids,False,'PairList')
    use_local_data = True
    my_list_imposter_run_ids = find_off_runs_around_source(obs_name,obs_ra,obs_dec,run_epoch,run_obs_type,run_elev_range,my_list_on_run_ids,True,'ImposterList')
    my_list_imposter_off_run_ids = find_off_runs_around_source(obs_name,obs_ra,obs_dec,run_epoch,run_obs_type,run_elev_range,my_list_imposter_run_ids,False,'ImposterPairList')


run_id = 103322
#run_id = 104633
#get_run_category(run_id)
#get_run_timecut(run_id)
#get_run_usable_duration(run_id)
#get_run_ra_dec(run_id)
#get_run_weather(run_id)
#get_run_type(run_id)
#get_run_nsb(run_id)
#get_run_el_az(run_id)

