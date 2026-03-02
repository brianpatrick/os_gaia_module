# PROCESS THE OPEN CLUSTER CATAOLOG:
# https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/673/A114
#
#
# ZACK REEVES
# CREATED: 2023
# CADE MOHRHARDT
#
# VERSIONS:
#  1.1  OCT 2023 CREATE JUPYTER NOTEBOOK
#  Python 3.12.12 OCT 2025

# Define the metadata for the data set. 
PRINT_CENSUS = True
GENERATE_ASSET_FILE = True

def generate_metadata():
    metadata = {}

    metadata['project'] = 'Digital Universe Atlas'
    metadata['sub_project'] = 'Open Clusters'

    metadata['catalog'] = 'Improving the open cluster census. II. An all-sky cluster catalogue with Gaia DR3 (Hunt+, 2023)'
    metadata['catalog_author'] = 'Hunt+'
    metadata['prepared_by'] = 'Zack Reeves (AMNH), Cade Mohrhardt (AMNH)'
    metadata['catalog_year'] = '2023' #Cade added this, not sure if this is what the catalog_year refers to
    metadata['version'] = '1.1'

    metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
    metadata['raw_data_dir'] = ''

    metadata['data_group_title'] = 'Open Clusters'
    metadata['data_group_desc'] = 'Open Cluster catalog'
    metadata['data_group_desc_long'] = "A catalog of stars from open clusters"
    metadata['fileroot'] = 'oc'
    return metadata

metadata = generate_metadata()

import pandas as pd
import numpy as np
import sys
import collections

from astropy.io import ascii
import astropy.units as u
import astropy.coordinates
from astropy.coordinates import Angle
from astropy.table import Table

from astroquery.vizier import Vizier

sys.path.insert(0, '..')
from common import file_functions, calculations, asset_creation

import matplotlib.pyplot as plt

#file_functions.generate_asset_file(metadata)
#Reading in the catalog with Vizier
#We specify the row limit to make sure we get all the stars in the catalog
#We place constraints on the Parallax and Probability of being a White Dwarf as a preliminary thresh
#We specify columns = ['**'] to get all of the columns, not just the default ones
catalog = Vizier(catalog='J/A+A/673/A114/clusters', columns=['**'], row_limit=-1).query_constraints(dist50='> 0.0')
data = catalog[0]


data

#calculating distance in light years and parsecs
calculations.get_distance(data, dist='dist50', use='distance')

#calculating cartesian coordinates
calculations.get_cartesian(data, ra='RA_ICRS', dec='DE_ICRS')

#playing around with threshing on distance
data.remove_rows(np.where(data['dist_pc']>20000)[0])

#2D Visualization
fig, ax = plt.subplots(1, 2)

#XY Plane
ax[0].scatter(data['x'], data['y'])
ax[0].set_title('XY Plane')

#XZ Plane
ax[1].scatter(data['x'], data['z'])
ax[1].set_title('XZ Plane')

#set good spacing
fig.tight_layout()
fig.set_size_inches(10, 4, forward=True)
plt.show

#construct a speck comment column
data['speck_label'] = data.Column(data=['#__'+name for name in data['Name']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                 description='Object ID')

#construct a label column
data['label'] = data['Name']  #leaving for now in case we want to add other labels

#construct a metadata table
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'dist_ly', 'N', 'r50', 'speck_label'])
columns

# Print the speck file using the to_speck function in file_functions
file_functions.to_speck(metadata, Table.to_pandas(data), columns)

# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(data))

# Print the csv file using the to_csv function in file_functions
file_functions.to_csv(metadata, Table.to_pandas(data), columns)


def asset_main():
    """Generate the asset file for Gaia Stellar Clusters."""

    datainfo = {
        # Core script requirements
        "renderable": "RenderablePolygonCloud",
        "filename": metadata['fileroot'],
        "Identifier": "GaiaOpenClusters",
        "local_modules": True,
        "asset_dir": "",
        "Enabled": "true",
        "name": "Gaia Open Clusters",
        "data": {
            "File": "oc.speck",
            "Name": "Open Clusters Speck Files",
            "Identifier": "gaia_openclusters_speck",
            "Version": 3
            },
        "Texture": {
            "File": "stars_textures"
        },
        
        "local colormaps": """asset.resource({
    Name = "Stars Color Table",
    Type = "HttpSynchronization",
    Identifier = "stars_colormap",
    Version = 3
})""",

        "local textures": """asset.resource({
    Name = "Stars Textures",
    Type = "HttpSynchronization",
    Identifier = "stars_textures",
    Version = 1
})""",

        # Renderable Parameters
        "opacity": 0.9,
        "PolygonSides": "12",
        "Unit": "pc",
        "FixedColor": "0.80, 0.0, 0.50",
        
        # Label Settings
        "labels": {
            "file": "oc.label",
            "color": "{ 0.0, 0.36, 0.14 }",
            "size": 15.5,
            "min_max": "{ 4, 30 }"
        },
        
        # Size Settings
        "ScaleExponent": "15.7",
        "MaxSize": "23.0",
        "enable_max_size": "true",

        # GUI & Meta
        "GUI": {
        "Path": "/Milky Way/Star Clusters",
        "Name": "Gaia Open Clusters",
        "Description": "This is a large catalog containing 'Milky Way open clusters' using Gaia DR3 data.",
        "url": "https://www.amnh.org/research/hayden-planetarium/digital-universe",
        "license": "AMNH Gaia"
        },
        "meta_name": metadata['sub_project'],
		"author": metadata['prepared_by']
    }

    # Generate the file using your internal utility
    asset_creation.write_asset(datainfo)

if __name__ == "__main__":
    asset_main()
    print("open clusters asset generated successfully.")