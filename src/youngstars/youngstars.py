# PROCESS THE YOUNG STARS CATAOLOG:
# https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/620/A172

# Zack Reeves
# Created: 2023
# Edited by Cade Mohrhardt: 2026

# Versions:
#  1.1  OCT 2023 CREATE JUPYTER NOTEBOOK
#  1.02  MAY 2026 MODERNIZE CODE TO MATCH THE DIGITAL UNIVERSE

#  Python 3.12.12

import pandas as pd
import numpy as np
import sys
import collections
from astropy.io import ascii
import astropy.units as u
import astropy.coordinates
from astropy.coordinates import Angle
from astropy.table import unique, vstack, Table, join
from astroquery.vizier import Vizier
from matplotlib import pyplot as plt, colors

sys.path.insert(0, '..')
from common import file_functions, calculations, get_bailer_jones, gaia_functions, asset_creation


GENERATE_SPECK = True
GENERATE_ASSET_FILE = True
READ_LOCAL_CATALOG = True


# Define the metadata for the data set. 
#https://www.aanda.org/articles/aa/full_html/2023/06/aa43964-22/aa43964-22.html
def generate_metadata():
    metadata = {}

    metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
    metadata['sub_project'] = 'Young Stars'

    metadata['catalog'] = '3D mapping of young stars in the solar neighbourhood with Gaia DR2 (Zari+, 2023)'  #need to edit
    metadata['catalog_author'] = 'Zari+'
    metadata['catalog_year'] = '2023'
    metadata['catalog_doi'] = 'doi:10.1051/0004-6361/202039498' #need to fix
    metadata['catalog_bibcode'] = '2021A&A...649A...6G' #need to fix

    metadata['prepared_by'] = 'Brian Abbott, Zack Reeves, Cade Mohrhardt'
    metadata['version'] = '1.1'

    metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
    metadata['raw_data_dir'] = ''

    metadata['data_group_title'] = 'Young Stars'
    metadata['data_group_desc'] = 'Young Stars'
    metadata['data_group_desc_long'] = 'Young Stars in the Milky Way mapped by Gaia'
    metadata['fileroot'] = 'young_stars'
    
    metadata['OS_identifier'] = metadata['sub_project']
    metadata['OS_gui_name'] = metadata['data_group_title']
    metadata['OS_gui_path'] = '/Milky Way/Stars'
    metadata['OS_gui_description'] = metadata['data_group_desc_long']

    return metadata


metadata = generate_metadata()

file_functions.generate_license_file(metadata)
#file_functions.generate_asset_file(metadata)



if READ_LOCAL_CATALOG == False:
    #Reading in the catalog with Vizier
    #We specify the row limit to make sure we get all the stars in the catalog
    #We place constraints on the Parallax as a preliminary thresh
    #We specify columns = ['**'] to get all of the columns, not just the default ones
    catalog = Vizier(catalog='J/A+A/620/A172', columns=['**'], row_limit=-1).query_constraints(Plx='> 0.0')

    #This catalog comes with 4 tables:
    # - Pre main sequence (has SIMBAD column)
    # - Upper main sequence (has SIMBAD column)
    # - Pre main sequence S=2 tangential velocity
    # - Pre main sequence S=3 tangential velocity

    #We first label each object with the table it came from
    catalog[0]['table'] = catalog[0].Column(data=['Pre-main sequence']*len(catalog[0]),
                                            meta = collections.OrderedDict([('ucd', 'meta.table')]),
                                            description='Catalog Table')

    catalog[1]['table'] = catalog[1].Column(data=['Upper main sequence']*len(catalog[1]),
                                            meta = collections.OrderedDict([('ucd', 'meta.table')]),
                                            description='Catalog Table')

    catalog[2]['table'] = catalog[2].Column(data=['Pre-main sequence S=2']*len(catalog[2]),
                                            meta = collections.OrderedDict([('ucd', 'meta.table')]),
                                            description='Catalog Table')

    catalog[3]['table'] = catalog[3].Column(data=['Pre-main sequence S=3']*len(catalog[3]),
                                            meta = collections.OrderedDict([('ucd', 'meta.table')]),
                                            description='Catalog Table')


    #We concatenate these tables into one for a full catalog
    #Some stars exist in multiple tables and present as duplicate objects
    #We remove duplicate objects using the unique function
    data = unique(vstack([catalog[0], catalog[1], catalog[2], catalog[3]], 
                metadata_conflicts='silent'), keys='Source', keep='first')
    
    #querying Gaia for bailer-jones distances
    distances = get_bailer_jones.get_bj_distances(data, source_id='Source')

    data = join(data, distances, keys='Source', join_type='inner')

    #Download the query results to a csv
    download = data.to_pandas()
    download.to_csv('raw_data/youngstarsquery.csv', index=False)

if READ_LOCAL_CATALOG == True:
    data = pd.read_csv('raw_data/youngstarsquery.csv')
    data = Table.from_pandas(data)



#setting units and metadata for important columns
data['GLON'] = data.MaskedColumn(data=data['GLON'], 
                                      unit=u.deg,
                                      meta = collections.OrderedDict([('ucd', 'pos.eq.ra')]),
                                      format='{:.6f}', 
                                      description='Right Ascension')

data['GLAT'] = data.MaskedColumn(data=data['GLAT'], 
                                      unit=u.deg,
                                      meta = collections.OrderedDict([('ucd', 'pos.eq.dec')]),
                                      format='{:.6f}', 
                                      description='Declination')

data['pmGLON'] = data.MaskedColumn(data=data['pmGLON'], 
                                       unit=u.mas/u.yr,
                                       meta = collections.OrderedDict([('ucd', 'pos.eq.ra')]),
                                       format='{:.6f}', 
                                       description='Proper Motion of RA')

data['pmGLAT'] = data.MaskedColumn(data=data['pmGLAT'], 
                                       unit=u.mas/u.yr,
                                       meta = collections.OrderedDict([('ucd', 'pos.eq.dec')]),
                                       format='{:.6f}', 
                                       description='Proper Motion of DEC')

data['RV'] = data.MaskedColumn(data=data['RV'], 
                                       unit=u.km/u.s,
                                       meta = collections.OrderedDict([('ucd', 'pos.eq.dec')]),
                                       format='{:.6f}', 
                                       description='Radial Velocity')


gaia_functions.set_bj_distance(data)

# #fixing parallax units (Vizier labels it as a magnitude, probably meant milliarcseconds (mag versus mas))
# data['Plx'].unit=u.mas

# #thresh on parallax error (cutting on >10% error removes 1614 stars)
# data['parallax_over_error'] = [data['Plx'][i] / data['e_Plx'][i] for i in range(len(data))]
# data.remove_rows(np.where(data['parallax_over_error']<10)[0])

#calculating distance in light years and parsecs
#calculations.get_distance(data, parallax='Plx', use='parallax')

calculations.get_distance(data, dist='bj_distance', use='distance')

len(data)

#threshing on distance
data.remove_rows(np.where(data['dist_pc']<0.10)[0])

#calculating cartesian coordinates
calculations.get_cartesian(data, glon='GLON', glat='GLAT', pmglon='pmGLON', pmglat='pmGLAT', 
                           radial_velocity='RV', frame='galactic')

gaia_functions.get_magnitudes(data, gmag='Gmag')
gaia_functions.get_luminosity(data)
data['bp_rp'] = [data['BPmag'][i] - data['RPmag'][i] for i in range(len(data))]
gaia_functions.get_bp_g_color(data, color='bp_rp')


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
# plt.show


#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')

#construct a speck comment column
data['speck_label'] = data.Column(data=['#  '+str(name) for name in data['Source']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                 description='Gaia DR2 Source ID')

#construct a label column
data['label'] = ['GaiaDR2_'+ str(source) for source in data['Source']]  #leaving for now in case we want to add other labels


#construct a metadata table
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'appmag', 'absmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'speck_label'])



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