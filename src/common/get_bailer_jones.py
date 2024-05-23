#function to upload a set of Gaia IDs to the gea archive and query the Gaia DR3 Bailer Jones Distances
#created by Zack Reeves February 2024

import numpy as np
from astropy.table import Table, join
from astroquery.gaia import Gaia
from astroquery.utils.tap.core import TapPlus
import astropy.units as u

#when calling the function, source_id should be set as whatever the Gaia EDR3, DR3, or DR2 ID is called in the table
#the IDs should be the raw id string, no prefix
#context should be DR2, EDR3, or DR3
def get_bj_distances(data:Table, source_id='source_id', columns=None, get_motion=False, context='DR3'):
    
    #log in to Gaia Server - Can change to different credentials file for a different user
    Gaia.login(credentials_file='../common/gaia_credentials.txt')
    
    #get username from credentials file for query
    with open('../common/gaia_credentials.txt', 'r') as file:
        username = file.readline()
        
    #set proper motion and radial velocity columns if needed
    if(get_motion):
        motion_columns='gs.pmra, gs.pmdec, gs.radial_velocity, '
    else:
        motion_columns=''
     
    #construct gaia source columns string
    gaia_source_columns = "gs.parallax, gs.parallax_error, "+motion_columns
    
    #choose the Gaia source table, Bailer-Jones table and desired columns based on data release context
    if((context=='DR3')or(context=='EDR3')):
        bj_table = 'external.gaiaedr3_distance'
        bj_columns = 'bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo'
        gaia_source_table='gaiadr3.gaia_source'
    elif(context=='DR2'):
        bj_table = 'external.gaiadr2_geometric_distance'
        bj_columns = 'bj.r_est, bj.r_hi, bj.r_lo'
        gaia_source_table='gaiadr2.gaia_source'
    else:
        raise Exception('Context must be DR2, EDR3, or DR3')
        
    #construct query
    query = "select a.*, "+gaia_source_columns+bj_columns+" "+"from user_"+username+".gaia_ids a left join "+gaia_source_table+" gs on a."+source_id+" = gs.source_id left join "+bj_table+" bj on a."+source_id+" = bj.source_id"

    
    #create a table with the columns specified
    #by default, this is just the Gaia source ID.  If other columns are desired, specify with the columns argument
    if(columns==None):
        id_table=data[[source_id]]
    else:
        id_table=data[columns]

    #Upload table (table name will be forced to lowercase)
    job = Gaia.upload_table(upload_resource=id_table,
                            table_name="gaia_ids", format="csv")

    #Query Gaia DR3 source for parallaxes
    job = Gaia.launch_job_async(query, dump_to_file=False)

    #Put the resulting table into a Table
    distances = job.get_results()
    
    #dropping the original index (artifact from the query)
    distances.remove_column('gaia_ids_oid')
    
    #Deleting table and job from Gaia ESA server so we don't clog the memory
    Gaia.delete_user_table('gaia_ids')
    Gaia.remove_jobs(job.jobid)

    Gaia.logout()
    
    #Choose distance we want (for EDR3 and DR3: dist>500 -> use photogeo if it exists, otherwise, geo)
    #also set dcalc: 1 for geometric distance, 2 for photogeometric
    if(context=='DR2'):
        distances['bj_distance'] = distances['r_est']
        distances['dcalc'] = [1]*len(distances)
    else:
        distances['bj_distance'] = [distances['r_med_photogeo'][i]*u.pc if((not(np.ma.is_masked(distances['r_med_photogeo'][i])))and(distances['r_med_geo'][i]>500)) else distances['r_med_geo'][i]*u.pc for i in range(len(distances))]
        distances['dcalc'] = [1 if((not(np.ma.is_masked(distances['r_med_photogeo'][i])))and(distances['r_med_geo'][i]>500)) else 2 for i in range(len(distances))]
    
    #Calculate distance error, same criteria as distance
    if(context=='DR2'):
        distances['e_bj_dist'] = [(distances['r_hi'][i]-distances['r_lo'][i])/2 for i in range(len(distances))]
    else:
        distances['e_bj_dist'] = [((distances['r_hi_photogeo'][i]-distances['r_lo_photogeo'][i])/2)*u.pc if((not(np.ma.is_masked(distances['r_med_photogeo'][i])))and(distances['r_med_geo'][i]>500)) else ((distances['r_hi_geo'][i]-distances['r_lo_geo'][i])/2)*u.pc for i in range(len(distances))]
    
    
    #Calculate distance error percentage
    distances['bj_error_over_distance'] = [distances['e_bj_dist'][i]/distances['bj_distance'][i] for i in range(len(distances))]
    
    return distances
    #return join(data, distances, keys=source_id, join_type='left')