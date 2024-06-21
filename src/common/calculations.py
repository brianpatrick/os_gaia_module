# calculations v.2
# created by Zack Reeves
# Sept 30, 2023
# Updated Nov 19, 2023

# Calculations necessary in processing data for the Digital Universe

# functions:

# get_distance() - takes an Astropy Table and calculates columns for distance in parsecs and distance in light years

# get_cartesian() - takes an Astropy Table and adds calculated columns for XYZ and UVW if given proper motions and radial velocity

#get_photometric_distance() - takes stellar effective temperature, radius, and Kepler magnitude to calculate a distance using astronomical equations

# get_num_nearby() - takes an Astropy Table and calculates the number of objects within the table that are within a specified 
#                    distance to each object

#get_redshift_distance() - calculates the lookback and comoving distances of objects in a table given redshifts

import pandas as pd
import numpy as np

from scipy.spatial import cKDTree

from astropy.table import Table
import astropy.coordinates
import astropy.units as u
from astropy import constants as const
import astropy.cosmology.units as cu
from astropy.cosmology import WMAP9

import collections

from tqdm import tqdm

# takes a distance held in 'data' and converts it to distances in parsecs and light years
# If both a parallax and a distance exist in the data, the parallax is used by default.  If distance is preferred, change the 'use' argument to 'distance'
# data must be an Astropy Table
# if given a parallax or a distance, the given column must have a unit attached
def get_distance(data:Table, parallax = 'Plx', dist='Dist', use='parallax'):
    #
    if((parallax not in data.columns)&(dist not in data.columns)):
        raise Exception('calculations.get_distance(): parallax or distance must be in data, neither were found')
    if((parallax in data.columns)&(use=='parallax')):    
        data['dist_pc'] = data.MaskedColumn(data=data[parallax].to(u.pc, equivalencies=u.parallax()), 
                                      meta=collections.OrderedDict([('ucd', 'pos.distance')]), 
                                      format='{:.6f}',
                                      description='Distance from Sun (pc)')
    elif((dist in data.columns)&(use=='distance')):
        data['dist_pc'] = data.MaskedColumn(data=data[dist].to(u.pc), 
                                      meta=collections.OrderedDict([('ucd', 'pos.distance')]), 
                                      format='{:.6f}',
                                      description='Distance from Sun (pc)')
    else:
        raise Exception('calculations.get_distance(): use argument must exist in data')
    
    #calculate distance in light years based on dist_pc
    data['dist_ly'] = data.MaskedColumn(data=data['dist_pc'].to(u.lyr), 
                                  meta=collections.OrderedDict([('ucd', 'pos.distance')]), 
                                  format='{:.1f}', 
                                  description='Distance from Sun (lyr)')

# transforms a given set of coordinates in a pandas df (RA/DEC, L/B) to Cartesian XYZ
# if given proper motions and radial velocities, also returns UVW and speed
def get_cartesian(data:Table, frame='icrs', dist='dist_pc', ra='ra', dec='dec', glon='GLON', glat='GLAT', pmra='pmra', pmde='pmdec', pmglon='pmglon', pmglat='pmglat', radial_velocity='radial_velocity', epoch='J2000'):
    
    #Raise exception if distance is not in data
    if(dist not in data.columns):
        raise Exception('distance must be provided')
        
    #setting the units based on the distance
    dist_unit = str(data[dist].unit)
    if(dist_unit=='pc'):
        dist_unit='parsecs'
    elif(dist_unit=='Mpc'):
        dist_unit='Megaparsecs'
    
    #We use the calculate_velocities indicator to inform the code whether or not proper motions and radial velocities have been provided.  It is defined as False unless pmra, pmde, and radial_velocity have been called
    calculate_velocities = False
    
    #define frame
    #The Galactic frame is ultimately necessary. If given an ICRS frame, we transform to Galactic before getting the Cartesian representation.
    #To calculate velocities, the ICRS frame must be given.  If in the ICRS frame, we check to see if proper motions and radial velocity are in data and set calculate_velocity accordingly
    if(frame=='icrs'):
        
        #raise exception if ra and dec are not given
        if((ra not in data.columns)|(dec not in data.columns)):
            raise Exception('RA and Dec must be provided to calculate ICRS position')
            
        if((pmra in data.columns)&(pmde in data.columns)&(radial_velocity in data.columns)):
            #set velocity variable to true to indicate that we will calculate velocities
            calculate_velocities=True  
            
    elif(frame=='galactic'):
        
        if((glon not in data.columns)|(glat not in data.columns)):
            raise Exception('GLON and GLAT must be provided to calculate Galactic position')
        else:
            #set glon and glat arrays
            #This is important, as if an icrs frame is given, glon and glat need to be calculated,
            #but if a galactic frame is given, these need to be used
            glon = data[glon]
            glat = data[glat]
            
        if((pmglon in data.columns)&(pmglat in data.columns)&(radial_velocity in data.columns)):
            #set velocity variable to true to indicate that we will calculate velocities
            calculate_velocities=True 
            
    else:
        raise Exception('\'frame\' argument must be \'icrs\' or \'galactic\'')
    
        
    #set proper motions and radial velocity to be arrays of np.nan if velocities won't be calculated
    #This is to ensure that the SkyCoord argument runs without issues
    if(not calculate_velocities):
        pmra=pmde=pmglon=pmglat=[np.nan*u.mas/u.yr]*len(data)
        radial_velocity = [np.nan*u.km/u.yr]*len(data)
    else:
        #A little scuffed but creating these variables are necessary to make the transforms work smoothly
        radial_velocity = data[radial_velocity]
        if(frame=='icrs'):
            pmra=data[pmra]
            pmde=data[pmde]
        elif(frame=='galactic'):
            pmglon=data[pmglon]
            pmglat=data[pmglat]
        
    
    #If in ICRS frame, transform to Galactic
    if(frame=='icrs'):
        #create array of ICRS SkyCoord objects
        icrs_coords = astropy.coordinates.SkyCoord(
            ra=data[ra],
            dec=data[dec],
            distance=data[dist],
            pm_ra_cosdec=pmra,
            pm_dec=pmde,
            radial_velocity=radial_velocity,
            frame='icrs'
            #obstime=epoch
            )
    
        #calculate galactic positions and proper motions
        glon = icrs_coords.galactic.l
        glat = icrs_coords.galactic.b
        pmglon = icrs_coords.galactic.pm_l_cosb
        pmglat = icrs_coords.galactic.pm_b
        

    #Calculate Cartesian representation from Galactic frame
    
    #Create array of Galactic SkyCoord objects
    galactic_coords = astropy.coordinates.SkyCoord(
        l=glon,
        b=glat,
        distance=data[dist],
        pm_l_cosb=pmglon,
        pm_b=pmglat,
        radial_velocity=radial_velocity,
        frame='galactic'
        )
    
    #get Cartesian representation and set metadata
    
    data['x'] = data.MaskedColumn(data=galactic_coords.cartesian.x, 
                            meta = collections.OrderedDict([('ucd', 'pos.cartesian.x')]),
                            format='{:.6f}', 
                            description='x position (galactic cartesian coordinates) in '+dist_unit)
    
    data['y'] = data.MaskedColumn(data=galactic_coords.cartesian.y, 
                            meta = collections.OrderedDict([('ucd', 'pos.cartesian.y')]),
                            format='{:.6f}', 
                            description='Position (y coordinate) in '+dist_unit)
    
    data['z'] = data.MaskedColumn(data=galactic_coords.cartesian.z, 
                            meta = collections.OrderedDict([('ucd', 'pos.cartesian.z')]),
                            format='{:.6f}', 
                            description='Position (z coordinate) in '+dist_unit)
    
    
    #if proper motions and rv was given, calculate UVW velocities and speed and set metadata
    if(calculate_velocities):
        
        data['u'] = data.MaskedColumn(data=galactic_coords.velocity.d_x, 
                                meta = collections.OrderedDict([('ucd', 'vel.cartesian.u')]),
                                format='{:.6f}', 
                                description='Heliocentric velocity towards Galactic Center')
        
        data['v'] = data.MaskedColumn(data=galactic_coords.velocity.d_y, 
                                meta = collections.OrderedDict([('ucd', 'vel.cartesian.v')]),
                                format='{:.6f}', 
                                description='Heliocentric velocity towards Galactic Rotation')
        
        data['w'] = data.MaskedColumn(data=galactic_coords.velocity.d_z, 
                                meta = collections.OrderedDict([('ucd', 'vel.cartesian.w')]),
                                format='{:.6f}', 
                                description='Heliocentric velocity towards Galactic North Pole')
        
        data['speed'] = data.MaskedColumn(data=[np.sqrt(data['u'][i]**2 + data['v'][i]**2 + data['w'][i]**2) for i in range(len(data))], 
                                    meta = collections.OrderedDict([('ucd', 'vel.speed')]), format='{:.6f}', 
                                    description='Total heliocentric velocity')
        
        
# Calculates photometric distance given stellar teff, radius, and apparent magnitude
# magnitude is assumed to be in the Kepler band
# Rstar is assumed to be in units of R_sun
def get_photometric_distance(Teff, Rstar, KEPmag):
    
    
    #calculate stellar luminosity in Watts
    Lstar = 4.0 * np.pi * Rstar**2 * const.R_sun.value**2 * const.sigma_sb.value * Teff**4.0
    
    #calculate absolute magnitude of star
    Mstar = 4.75 - 2.5 * np.log(Lstar/const.L_sun.value)/np.log(10.0)
    
    #calculate and return distance in parsecs
    dist_pc = 10.0**((KEPmag - Mstar + 5.0) / 5.0)*u.pc
    return dist_pc

# Calculates the number of nearby objects for each object in an Astropy Table
# Assumes Astropy table data has columns x, y, and z calculated and accordingly named
# Distance factor should be in same units as XYZ and distance
def get_num_nearby(data:Table, distance_factor:float, dist='comoving_distance'):
    #Thank you ChatGPT <3

    #Create dataframe with data
    df = Table.to_pandas(data)

    # Extract the cartesian coordinates and distances
    coordinates = df[['x', 'y', 'z']].values
    distances = df['comoving_distance'].values

    # Build a KD-tree for efficient spatial queries
    kdtree = cKDTree(coordinates)

    # For each galaxy, query the KD-tree to find neighbors within the distance factor
    num_nearby_galaxies = []
    for i in tqdm(range(len(coordinates)), desc='Processing', position=0, leave=True):
        nearby_indices = kdtree.query_ball_point(coordinates[i], distance_factor)
        # Exclude the galaxy itself
        num_nearby = len(nearby_indices) - 1
        num_nearby_galaxies.append(num_nearby)

    # Add the results as a new column to your DataFrame
    data['num_nearby_galaxies'] = data.Column(data=num_nearby_galaxies,
                                              meta=collections.OrderedDict([('ucd', 'meta.number')]),
                                              description='Number of nearby galaxies') #amend to add number of parsecs considered


#calculates the lookback and comoving distances of objects in a table given redshifts
def get_redshift_distance(data:Table, redshift='z'):
    #raise exception if redshift is not in data - could also mean that redshift is named differently
    if(redshift not in data.columns):
        raise Exception('redshift must exist in data')
    
    #calculating lookback and comoving distances
    lookback = data[redshift].to(u.lyr, cu.redshift_distance(WMAP9, kind="lookback"))
    comoving = data[redshift].to(u.Mpc, cu.redshift_distance(WMAP9, kind="comoving"))
    
    #calculating lookback time in Gyrs
    lookback_time = lookback.value / 10**9
    
    #setting columns and metadata
    data['lookback_time'] = data.MaskedColumn(data=lookback_time,
                                              unit=u.Gyr,
                                              meta = collections.OrderedDict([('ucd', 'time.lookback')]),
                                              format='{:.6f}', 
                                              description='Redshift-based lookback time')
    data['comoving_distance'] = data.MaskedColumn(data=comoving,
                                                  meta = collections.OrderedDict([('ucd', 'pos.distance.comoving')]),
                                                  format='{:.6f}', 
                                                  description='Redshift-based comoving distance')



           
    