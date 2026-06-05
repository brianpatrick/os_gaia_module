# PROCESS THE GAIA CATALOG OF STARS with CHEMISTRY from APOGEE and GALAH:
#currently just has apogee

# Zack Reeves
# Created: 2024
# Edited by Cade Mohrhardt: 2026

# Versions:
#  1.1  JUN 2024 CREATE JUPYTER NOTEBOOK
#  1.02  MAY 2026 MODERNIZE CODE TO MATCH THE DIGITAL UNIVERSE

#  Python 3.12.12

import pandas as pd
import numpy as np
import sys
import os
import collections
import astropy.units as u
import astropy.coordinates
from astropy.table import Table, join, unique, vstack
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import SciServer
from SciServer import CasJobs
from matplotlib import pyplot as plt, colors

sys.path.insert(0, '..')
from common import file_functions, calculations, gaia_functions, asset_creation


GENERATE_SPECK = True
GENERATE_ASSET_FILE = True
READ_LOCAL_CATALOG = False


# Define the metadata for the data set. 
def generate_metadata():
    metadata = {}

    metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
    metadata['sub_project'] = 'Gaia Catalog of Nearby Stars'

    metadata['catalog'] = 'The Gaia Catalogue of Nearby Stars (Gaia Collaboration, 2021)'
    metadata['catalog_author'] = 'Gaia Collaboration'
    metadata['catalog_year'] = '2021'
    metadata['catalog_doi'] = 'doi:10.1051/0004-6361/202039498'
    metadata['catalog_bibcode'] = '2021A&A...649A...6G'

    metadata['prepared_by'] = 'Brian Abbott, Zack Reeves, Cade Mohrhardt'
    metadata['version'] = '1.1'

    metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
    metadata['raw_data_dir'] = ''

    metadata['data_group_title'] = 'Chemistry'
    metadata['data_group_desc'] = 'Chemistry of stars in the Milky Way mapped by Gaia and measured by APOGEE/GALAH'
    metadata['data_group_desc_long'] = 'Have you ever wondered what is out there in space? Now, thanks to Gaia EDR3, the solar neighbourhood has been mapped with great precision out to 100 pc (326 light years)'
    metadata['fileroot'] = 'chem'

    metadata['OS_identifier'] = metadata['sub_project']
    metadata['OS_gui_name'] = metadata['data_group_title']
    metadata['OS_gui_path'] = '/Milky Way/Stars'
    metadata['OS_gui_description'] = metadata['data_group_desc_long']

    return metadata


metadata = generate_metadata()

file_functions.generate_license_file(metadata)
#file_functions.generate_asset_file(metadata)



if READ_LOCAL_CATALOG == False:
    #logging into CasJobs
    #to make a new account: https://apps.sciserver.org/login-portal/Account/Login
    from SciServer import Authentication

    #SciServer_credentials.txt should have your username on line 0 and your password on line 1
    file = open('../common/SciServer_credentials.txt', 'r')
    lines = file.readlines()
    file.close()
    Authentication_loginName = lines[0].strip()
    Authentication_loginPassword = lines[1].strip()

    manualtoken = Authentication.login(Authentication_loginName, Authentication_loginPassword)
    manualtokenvalue = Authentication.token.value

    #Querying the aspcapStar table from SDSS to get chemistry info

    #we select the apogee_id (2MASS style ID) and any chemistry info we want
    #we may want to thresh the quality of the chemistry data
    query = 'select apogee_id as apogee_id, fe_h, fe_h_flag, si_fe, alpha_m from dr18.aspcapStar' ###change top 100 limiter
    apogee = Table.from_pandas(CasJobs.executeQuery(query, context='DR18', format='pandas'))


    #creating a new column of apogee IDS without the 2M string in front
    apogee['twomass_id'] = [apogee['apogee_id'][i][2:] if apogee['apogee_id'][i][:2]=='2M' else apogee['apogee_id'][i] for i in range(len(apogee))]


    apogee['survey'] = ['apogee']*len(apogee)

    # #retrieve the GALAH DR3 dataset from Vizier
    # #reading in the catalogue
    catalog = Vizier(catalog='J/MNRAS/506/150/stars', columns=['**'], row_limit=-1).query_constraints()
    ###
    print(catalog[0])
    ###
    # #choosing our GALAH columns
    # #we take the 2MASS id and any chemistry data we want
    galah = catalog[0][['2MASS', '[C/Fe]']] ### check that [C/Fe] is same as '__c_fe_', which it was previously
    galah.rename_column('2MASS', 'twomass_id')

    galah['survey'] = ['galah']*len(galah)

    #threshing galah data
    data = vstack([apogee, galah])


    #Query Gaia ESA ADQL server using apogee_ids to match to Gaia stars and obtain proper motion to calculate uvw as well as photometric data
    #log in to Gaia Server - Can change to different credentials file for a different user
    Gaia.login(credentials_file='../common/gaia_credentials.txt')

    #grab username from file
    file = open('../common/gaia_credentials.txt', 'r')
    username = file.readline().strip()

    #Delete a table that may have been created on a previous attempt
    #Gaia.delete_user_table('chemistry_stars')
    # #Upload table (table name will be forced to lowercase)
    job = Gaia.upload_table(upload_resource=apogee[['twomass_id']], table_name="chemistry_stars", format="csv")

    #Query Gaia DR3 source for parallaxes
    #Potentially want Bailer Jones distances pending figuring out the parallax error issue
    job = Gaia.launch_job_async("select a.twomass_id, "
                                "gaia_match.source_id, "
                                "bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo, "
                                "c.ra, c.dec, c.pmra, c.pmdec, c.radial_velocity, c.phot_g_mean_mag, c.bp_g, c.teff_gspphot "
                                "from user_"+username+".chemistry_stars a inner join gaiadr3.tmass_psc_xsc_best_neighbour gaia_match on a.twomass_id = gaia_match.original_ext_source_id "
                                "inner join external.gaiaedr3_distance bj on gaia_match.source_id = bj.source_id "
                                "inner join gaiadr3.gaia_source c on gaia_match.source_id = c.source_id "
                                "where c.phot_g_mean_mag > 0",
                                dump_to_file=False)

    #put the resulting table into a dataframe and fix the dtype of the apogee_id column
    query_data = job.get_results()
    query_data['twomass_id'] = query_data.Column(data=[str(query_data['twomass_id'][i]) for i in range(len(query_data))],
                                                    meta=collections.OrderedDict([('ucd', 'meta.id.tmass')]), 
                                                    description='2MASS ID')

    #Deleting table and job from Gaia ESA server so we don't clog the memory
    Gaia.delete_user_table(table_name="user_"+username+".chemistry_stars")
    Gaia.remove_jobs(job.jobid)

    Gaia.logout()

    #join table onto data
    data = unique(join(apogee, query_data, keys='twomass_id', join_type='inner'), keys=['twomass_id', 'source_id'])

    data.remove_rows(np.where(data['fe_h']<-999)[0])
    data.remove_rows(np.where(data['alpha_m']<-999)[0])
    data.remove_rows(np.where(data['si_fe']<-999)[0])

    # Uncomment to download the query results to a csv
    download = data.to_pandas()
    download.to_csv('raw_data/chemistryquery.csv', index=False)

if READ_LOCAL_CATALOG == True:
    data = pd.read_csv('raw_data/chemistryquery.csv')
    data = Table.from_pandas(data)



len(data[data['si_fe']<-999])

#setting dcalc based on r_med_geo (if>500pc and photogeo exists, we choose photogeo and set dcalc to 1, else geo and dcalc to 2)
data['dcalc'] = [1 if((not(np.ma.is_masked(data['r_med_photogeo'][i])))and(data['r_med_geo'][i]>500)) else 2 for i in range(len(data))]

#setting metadata for dcalc
data['dcalc'] = data.Column(data['dcalc'],
                            meta=collections.OrderedDict([('ucd', 'meta.dcalc')]),
                            description='Distance Indicator: 1 indicates a Bailer-Jones photogeometric distance; 2 indicates a Bailer-Jones geometric distance')

#Choosing distance based on dcalc
data['bj_distance'] = [data['r_med_photogeo'][i] if data['dcalc'][i]==1 else data['r_med_geo'][i] for i in range(len(data))]
data['bj_distance'].unit=u.pc

#Choosing and calculating distance error based on the distance we chose
data['e_bj_dist'] = [((data['r_hi_photogeo'][i]-data['r_lo_photogeo'][i])/2)*u.pc if((not(np.ma.is_masked(data['r_med_photogeo'][i])))and(data['r_med_geo'][i]>500)) else ((data['r_hi_geo'][i]-data['r_lo_geo'][i])/2)*u.pc for i in range(len(data))]


#calculating distance in light years and parsecs
calculations.get_distance(data, dist='bj_distance', use='distance')

#calculating cartesian coordinates
calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='radial_velocity', frame='icrs')

gaia_functions.get_magnitudes(data)


#calculate luminosity based on absolute magnitude
data['lum'] = [10**(1.89 - 0.4*data['absmag'][i]) for i in range(len(data))]
small_luminosities = np.where((data['lum']>0.0) & (data['lum']<0.001))[0]
data['lum'][small_luminosities] = [0.001]*len(small_luminosities)

data['lum'] = data.MaskedColumn(data=data['lum'],
                             unit=u.solLum,
                             meta=collections.OrderedDict([('ucd', 'phys.luminosity')]),
                             format='{:.6f}',
                             description='Stellar Luminosity')

#setting color and visualizing
data['color'] = data.MaskedColumn(data=data['bp_g'],
                             unit=u.solLum,
                             meta=collections.OrderedDict([('ucd', 'phys.color')]),
                             format='{:.2f}',
                             description='Gaia BP-G color')


# plt.hist(data['color'], bins=250);

# plt.hist(data['fe_h'], bins = 250);

# #2D Visualization
# fig, ax = plt.subplots(1, 2)

# #XY Plane
# ax[0].scatter(data['x'], data['y'])
# ax[0].set_title('XY Plane')

# #XZ Plane
# ax[1].scatter(data['x'], data['z'])
# ax[1].set_title('XZ Plane')

# #set good spacing
# fig.tight_layout()
# fig.set_size_inches(10, 4, forward=True)
# plt.show

# #2D Density Visualization
# fig, ax = plt.subplots(1, 2)

# #XY Plane
# ax[0].hist2d(data['x'], data['y'], 
#            bins = 200,  
#            norm = colors.LogNorm(),  
#            cmap = "RdYlGn_r",) 
# ax[0].set_title('XY Plane')

# #XZ Plane
# ax[1].hist2d(data['x'], data['z'], 
#            bins = 200,  
#            norm = colors.LogNorm(),  
#            cmap = "RdYlGn_r",) 
# ax[1].set_title('XZ Plane')

# #set good spacing
# fig.tight_layout()
# fig.set_size_inches(10, 4, forward=True)
# #plt.show


#construct a speck comment column
data['speck_label'] = data.Column(data=['#  '+str(name) for name in data['source_id']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                  description='Gaia DR3 Source ID')

#construct a label column
data['label'] = ['GaiaDR3_'+ str(source) for source in data['source_id']]  #leaving for now in case we want to add other labels


#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')


data['fe_h'] = data.Column(data=data['fe_h'], 
                           unit=u.dex,
                           meta=collections.OrderedDict([('ucd', 'spectroscopy.metallicity')]),
                           description='APOGEE Metallicity')

data['alpha_m'] = data.Column(data=data['alpha_m'], 
                           unit=u.dex,
                           meta=collections.OrderedDict([('ucd', 'spectroscopy.metallicity')]),
                           description='APOGEE [alpha/Fe]')

data['si_fe'] = data.Column(data=data['si_fe'], 
                           unit=u.dex,
                           meta=collections.OrderedDict([('ucd', 'spectroscopy.metallicity')]),
                           description='APOGEE [Si/Fe]')


#Getting the column metadata
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'fe_h', 'alpha_m', 'si_fe', 'speck_label'])



if GENERATE_SPECK:
    # Print the speck file using the to_speck function in file_functions
    file_functions.to_speck(metadata, Table.to_pandas(data), columns)


# Print the csv file using the to_csv function in file_functions
file_functions.to_csv(metadata, Table.to_pandas(data), columns)

# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(data))


df = Table.to_pandas(data)
file_functions.generate_plot_pdf(df[columns['name']], metadata)



def asset_main():
    """Generate the asset file for stars"""

    metadata = generate_metadata()
    datainfo = {
        "renderable": "RenderableStars",
        "filename": metadata['fileroot'],
		"asset_dir": "",
        "local_modules": True,
        "data": {
            "File": metadata['fileroot']+".speck",
            "Name": metadata['data_group_title']+" Speck Files",
            "Identifier": "gaia_"+metadata['fileroot']+"_speck",
            "Version": 6
            },
        "Texture": {
            "Glare": "halo.png",
            "Core": "glare.png",
            "Name": "Stars Textures",
            "Identifier": "stars_textures",
            "Version": 1
            },
        "ColorMap":{
            "ColorMap": "colorbv.cmap",
            "OtherDataColorMap": "viridis.cmap",
            "Name": "Stars Color Table",
            "Identifier": "stars_colormap",
            "Version": 3
            },
        "Identifier": metadata['fileroot'],
        "Bv_column": "color",
        "Luminance_column": "lum",
        "AbsoluteMagnitude_column": "absmag",
        "ApparentMagnitude_column": "appmag",
        "Vx_column": "u",
        "Vy_column": "v",
        "Vz_column": "w",
        "Speed_column": "speed",
        "GUI": {
            "Name": metadata['OS_gui_name'],
            "Path": metadata['OS_gui_path'],
            "Description": metadata['OS_gui_description']
            },
        "meta_name": metadata['sub_project'],
		"author": metadata['prepared_by']
    }
    asset_creation.write_asset(datainfo)


if __name__ == "__main__" and GENERATE_ASSET_FILE:
    asset_main()
    print("Asset file for "+metadata['data_group_title']+" generated successfully.")