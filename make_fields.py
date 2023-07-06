from datetime import datetime
from awlofar.database.Context import context
from awlofar.main.aweimports import CorrelatedDataProduct, \
    FileObject,  Observation, SubArrayPointing, AveragingPipeline
from awlofar.toolbox.LtaStager import LtaStager, LtaStagerError
import os,sys,time,argparse,numpy as np,multiprocessing as mp,json
import surveys_db; from surveys_db import SurveysDB

## all LoTSS project IDs and co-observing project IDs
projects = [ 'LC2_038', 'LC3_008', 'LC4_034', 'LT5_007', 'LC6_015', 'LC7_024', 'LC8_022', 'LC9_030', 'LT10_010', 'LC8_014', 'LC8_030', 'DDT9_001', 'LC9_011', 'LC9_012', 'LC9_019', 'LC9_020', 'COM10_001', 'LC10_001', 'LC10_010', 'LC10_014', 'LT10_012', 'LC11_013', 'LC11_016', 'LC11_019', 'LC11_020', 'LC12_014' ]

## possible calibrator names
calnames = ['3C196','3C295','3C048','3C48','3C147','3C380','3C287','4C59.16','3C244.1','Calibrator']

calsolutions = None
try:
    calsolutions = np.loadtxt('calibrator_solutions.dat',dtype='str')
except:
    print ('**********  NO CALSOLUTION POINTERS ************')
    print ('** To include these provide calsolutions.dat **')
    
fields = []

def fillobs (this, ob):
    this['lotssfield']=ob['field'] if ob is not None else ''
    this['lotssstatus']=ob['status'] if ob is not None else ''
    this['date']=ob['date'].ctime() if ob is not None else ''
    this['integration']=ob['integration'] if ob is not None else ''
    this['nchan']=ob['nchan'] if ob is not None else ''
    this['nsb']=ob['nsb'] if ob is not None else ''
    this['location']=ob['location'] if ob is not None else ''
    this['calibrator_id']=ob['calibrator_id'] if ob is not None else ''
    this['calibrator_nsb']=ob['calibrator_nsb'] if ob is not None else ''
    this['calibrator_name']=ob['calibrator_name'] if ob is not None else ''
    this['calibrator_date']=ob['calibrator_date'].ctime() if ob is not None else ''
    this['nr_international']=ob['nr_international'] if ob is not None else ''
    this['nr_core']=ob['nr_core'] if ob is not None else ''
    this['nr_remote']=ob['nr_remote'] if ob is not None else ''
    this['int_processed']=ob['int_processed'] if ob is not None else ''

def findsas(product, npoint):
    if npoint not in [2,3]:
        print ('Not 2 or 3 pointings!!!')
        return []
    else:
        if npoint == 2:
            sys.stdout.write ('Warning 2 pointings: ')
        nsub = len(product)
        sas = product[0].dataProductIdentifierName.split('_')[0].split('L')[1]
        ra = product[0].subArrayPointing.pointing.rightAscension
        dec = product[0].subArrayPointing.pointing.declination
        for i in range(nsub-1,0,-1):
            sys.stdout.write('%d.'%i); sys.stdout.flush()
            sas1 = product[i].dataProductIdentifierName.split('_')[0].split('L')[1]
            ra1 = product[i].subArrayPointing.pointing.rightAscension
            dec1 = product[i].subArrayPointing.pointing.declination
            if sas1!=sas:
                sys.stdout.write('\n')
                return [int(sas), float(ra), float(dec)], [int(sas1), float(ra1), float(dec1)]
        print ('**** WARNING returning only one SAS-ID ****')
        return [[int(sas), float(ra), float(dec)]]
        
                                                                       

print('Looking up Surveys observations and fields databases...')
sdb = SurveysDB(readonly=True)
sdb.execute('select * from observations')
obs_results = sdb.cur.fetchall()
sdb.close()
sdb = SurveysDB(readonly=True)
sdb.execute('select * from fields')
fld_results = sdb.cur.fetchall()
sdb.close()
os.system ('rm temp.txt; rm lc_surveysdb.txt')
obsdict = dict()

#for project in projects:
project = 'LC3_008'
if 1:
    print('Beginning project: %s, reading LTA table'%project)
    fo = open('lc_surveysdb.txt','a')
    for i in range(len(obs_results)):
        obs_result = obs_results[i]
        if obs_result['project_code']==None:
            continue
        if obs_result['project_code']==project:
            n=fo.write('%s %s %s %s %s %s\n'%\
                (obs_result['project_code'],obs_result['id'],\
                 obs_result['field'],obs_result['status'],\
                 obs_result['calibrator_id'],\
                 obs_result['calibrator_name']))
    fo.close()
    cls = CorrelatedDataProduct     # observation from awlofar
    query_observations = Observation.select_all().project_only(project)
    query_pipelines = AveragingPipeline.select_all().project_only(project)
    pnames, psas = np.array([],dtype='str'),np.array([],dtype='str')
    for pipeline in query_pipelines:
        pnames = np.append(pnames, pipeline.pipelineName)
        psas = np.append(psas, pipeline.observationId)
    maxfreq = 168.5
    fo = open('temp.txt','a')
#    observation = query_observations[0]
    for observation in query_observations:
        obsname_long = observation.observationDescription
        obsname = obsname_long.split('/')[0]
        if obsname in calnames:
            continue
        sys.stdout.write('Processing '+project+':'+obsname+': '); sys.stdout.flush()
        dataproduct_query = cls.observations.contains(observation)
        dataproduct_query &= cls.isValid==1
        dataproduct_query &= cls.maximumFrequency < maxfreq
        dataproduct = dataproduct_query[0]
        sas_id = findsas(dataproduct_query,\
                    observation.numberOfSubArrayPointings)
        n=fo.write('%s %s %6.2f %5.2f %s\n'%\
        (observation.get_project(),sas_id[0][0],sas_id[0][1],sas_id[0][2],obsname))
        if len(sas_id)>1:
            n=fo.write('%s %s %6.2f %5.2f %s\n'%\
            (observation.get_project(),sas_id[1][0],sas_id[1][1],sas_id[1][2],obsname))
        fo.flush()
        # do not iterate over the dataproduct_query to find sasId, it takes ages
        for sas in sas_id:
            this_obs = dict()
            this_obs['sas'] = sas[0]
            this_obs['ra'] = sas[1]
            this_obs['dec'] = sas[2]
            this_obs['project'] = project
            this_obs['description'] = observation.observationDescription
            if calsolutions is not None:
                calsas = np.asarray(calsolutions[:,0],dtype='int')
                try:
                    calsas_idx = np.argwhere(calsas==sas[0])[0][0]
                    this_obs['linccal']=calsolutions[calsas_idx,1]
                except:
                    this_obs['linccal']=None
            fillobs (this_obs, None)
            for obsr in obs_results:
                if obsr['id']==sas[0]:
                    fillobs (this_obs,obsr)
            found_existing_field = False
            for field in fields:
                decsep = sas[2]-field['dec']
                rasep = sas[1]-field['ra']
                sep = np.hypot(decsep,rasep*np.cos(np.rad2deg(field['dec'])))
                if sep<0.1:  # new observation in old field
                    field['obs'].append(this_obs)
                    field['nobs']+=1
                    found_existing_field = True
            if not found_existing_field:
                field=dict()
                field['ra'], field['dec'] = sas[1],sas[2]
                field['obs']=[]
                field['obs'].append(this_obs)
                field['nobs']=1
                fields.append(field)
    
    fo.close()
    os.system('sort -k 2 -n temp.txt >lc_lta.txt')
    os.system('rm temp.txt;rm fields.dat')
    json.dump (fields, open('fields.dat','w'))
