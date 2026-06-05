# PROCESS THE BAILER-JONES NEAR ENCOUNTERS CATALOG:
# https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/935/L9

# Zack Reeves
# Created: 2024
# Edited by Cade Mohrhardt: 2026
#
# Versions:
#  1.1  MAR 2024 CREATE JUPYTER NOTEBOOK
#  1.02 MAY 2026 MODERNIZE CODE TO MATCH THE DIGITAL UNIVERSE

#  Python 3.12.12 OCT 2025

import pandas as pd
import numpy as np
import sys
import os
import collections
import astropy.units as u
import astropy.coordinates
from astropy.table import Table, join
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
from matplotlib import pyplot as plt, colors

sys.path.insert(0, '..')
from common import file_functions, calculations, asset_creation


GENERATE_SPECK = True
GENERATE_ASSET_FILE = True
READ_LOCAL_CATALOG = True


# Define the metadata for the data set. 
def generate_metadata():
    metadata = {}

    metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
    metadata['sub_project'] = 'Close encounters to the Sun'

    metadata['catalog'] = 'Close encounters to the Sun in Gaia DR3 (Bailer-Jones, 2022)'
    metadata['catalog_author'] = 'Bailer-Jones'
    metadata['catalog_year'] = '2022'
    metadata['catalog_doi'] = 'doi:10.26093/cds/vizier.19359009'
    metadata['catalog_bibcode'] = '2022yCat..19359009B'

    metadata['prepared_by'] = 'Zack Reeves (AMNH), Cade Mohrhardt (AMNH)'
    metadata['version'] = '1.1'

    metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
    metadata['raw_data_dir'] = ''

    metadata['data_group_title'] = 'Near Encounters to the Sun'
    metadata['data_group_desc'] = 'Near Encounters to the Sun'
    metadata['data_group_desc_long'] = 'This is a catalog of stars that pass closer than 1pc to the Sun within the past or future 6Myr.'
    metadata['fileroot'] = 'near_encounters'

    metadata['OS_identifier'] = metadata['sub_project']
    metadata['OS_gui_name'] = metadata['data_group_title']
    metadata['OS_gui_path'] = '/Milky Way/Stars'
    metadata['OS_gui_description'] = metadata['data_group_desc_long']

    return metadata


metadata = generate_metadata()

file_functions.generate_license_file(metadata)
#file_functions.generate_asset_file(metadata)



if READ_LOCAL_CATALOG == False:
    #reading in the catalogue
    catalog = Vizier(catalog='J/ApJ/935/L9/table12', columns=['**'], row_limit=-1).query_constraints()
    catalog[0]

    #reducing data down to the necessary columns
    data = catalog[0][['GaiaDR3', 'tphmed', 'dphmed', 'vphmed', 'Plx', 'e_Plx', 'RV', 'Gmag', 'GMAG', 'GLON', 'GLAT']]

    #Query Gaia ESA ADQL server using Gaia EDR3 IDs to obtain proper motion and RA/DEC to calculate uvw

    #log in to Gaia Server - Can change to different credentials file for a different user
    Gaia.login(credentials_file='../common/gaia_credentials.txt')

    #grab username from file
    file = open('../common/gaia_credentials.txt', 'r')
    username = file.readline().strip()

    #Upload table (table name will be forced to lowercase)
    job = Gaia.upload_table(upload_resource=data[['GaiaDR3']], table_name="near_encounters", format="csv")

    #Query Gaia DR3 source for parallaxes
    #Potentially want Bailer Jones distances pending figuring out the parallax error issue
    job = Gaia.launch_job_async("select a.GaiaDR3, b.ra, b.dec, b.pmra, b.pmdec, bp_g "
                                "from user_"+username+".near_encounters a left join gaiadr3.gaia_source b on a.GaiaDR3 = b.source_id ",
                                dump_to_file=False)

    #put the resulting table into a dataframe and drop the unnecessary index column
    data = join(data, job.get_results(), keys='GaiaDR3', join_type='left')
    
    #Download the query results to a csv
    download = data.to_pandas()
    download.to_csv('raw_data/nearencountersquery.csv', index=False)

    #Deleting table and job from Gaia ESA server so we don't clog the memory
    Gaia.delete_user_table(table_name="user_"+username+".near_encounters")
    Gaia.remove_jobs(job.jobid)

    Gaia.logout()

if READ_LOCAL_CATALOG == True:
    data = pd.read_csv('raw_data/nearencountersquery.csv')
    data = Table.from_pandas(data)



#setting units and metadata for important columns
data['ra'] = data.MaskedColumn(data=data['ra'], 
                                       unit=u.deg,
                                       meta = collections.OrderedDict([('ucd', 'pos.eq.ra')]),
                                       format='{:.6f}', 
                                       description='Right Ascension')

data['dec'] = data.MaskedColumn(data=data['dec'], 
                                       unit=u.deg,
                                       meta = collections.OrderedDict([('ucd', 'pos.eq.dec')]),
                                       format='{:.6f}', 
                                       description='Declination')

data['pmra'] = data.MaskedColumn(data=data['pmra'], 
                                        unit=u.mas/u.yr,
                                        meta = collections.OrderedDict([('ucd', 'pos.eq.ra')]),
                                        format='{:.6f}', 
                                        description='Proper Motion of RA')

data['pmdec'] = data.MaskedColumn(data=data['pmdec'], 
                                        unit=u.mas/u.yr,
                                        meta = collections.OrderedDict([('ucd', 'pos.eq.dec')]),
                                        format='{:.6f}', 
                                        description='Proper Motion of DEC')

data['RV'] = data.MaskedColumn(data=data['RV'], 
                                        unit=u.km/u.s,
                                        meta = collections.OrderedDict([('ucd', 'pos.eq.dec')]),
                                        format='{:.6f}', 
                                        description='Radial Velocity')

data['Plx'] = data.MaskedColumn(data=data['Plx'], 
                                        unit=u.mas,
                                        meta = collections.OrderedDict([('ucd', 'pos.parallax.trig')]),
                                        format='{:.6f}', 
                                        description='Median parallax estimated using the selected members')

data['tphmed'] = data.MaskedColumn(data=data['tphmed'], 
                                       unit=u.kyr,
                                       meta = collections.OrderedDict([('ucd', 'src.orbital.tp')]),
                                       format='{:.6f}', 
                                       description='Median perihelion time')

data['dphmed'] = data.MaskedColumn(data=data['dphmed'], 
                                       unit=u.pc,
                                       meta = collections.OrderedDict([('ucd', 'src.orbital.distPeriastron')]),
                                       format='{:.6f}', 
                                       description='Median perihelion distance')

data['vphmed'] = data.MaskedColumn(data=data['vphmed'], 
                                       unit=u.km/u.s,
                                       meta = collections.OrderedDict([('ucd', 'src.orbital.velocity')]),
                                       format='{:.6f}', 
                                       description='Median perihelion velocity')



#calculating distance in light years and parsecs
calculations.get_distance(data, parallax='Plx', use='parallax')

#calculating cartesian coordinates
calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='RV', frame='icrs')


#setting dcalc
#since we calculate distance only using parallax in this dataset, dcalc is always 2
data['dcalc'] = data.Column(data=[2]*len(data),
                            meta=collections.OrderedDict([('ucd', 'meta.dcalc')]),
                            description='Distance Indicator: 1 indicates a Bailer-Jones photogeometric distance; 2 indicates a Bailer-Jones geometric distance')


#calculating absolute magnitudes
#calculate absolute V mag based on apparent magnitude and distance
data['appmag'] = data.MaskedColumn(data=data['Gmag'],
                             unit=u.mag,
                             meta=collections.OrderedDict([('ucd', 'phot.mag;em.opt.G')]),
                             format='{:.6f}',
                             description='Apparent magnitude in Gaia G-band')
data['absmag'] = data.MaskedColumn(data=[data['appmag'][i]+5-5*np.log10(data['dist_pc'][i]) for i in range(len(data))],
                             unit=u.mag,
                             meta=collections.OrderedDict([('ucd', 'phot.magAbs;em.opt.G')]),
                             format='{:.6f}',
                             description='Absolute magnitude in Gaia G-band')

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
plt.hist(data['color'], bins=10)


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


data['error_over_parallax']=[data['e_Plx'][i]/data['Plx'][i] for i in range(len(data))]

len(data[data['error_over_parallax']>0.2])


#construct a speck comment column
data['speck_label'] = data.Column(data=['"'+str(name)+'"' for name in data['GaiaDR3']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                  description='Gaia DR3 Source ID')

#construct a label column
data['label'] = ['GaiaDR3_'+ str(source) for source in data['GaiaDR3']]  #leaving for now in case we want to add other labels

#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')


#Getting the column metadata
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'tphmed', 'dphmed', 'vphmed', 'speck_label'])



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