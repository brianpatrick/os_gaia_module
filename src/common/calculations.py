# calculations v.2
# created by Zack Reeves
# Sept 30, 2023

# Calculations necessary in processing data for the Digital Universe

# functions:

# get_distance() - takes an Astropy Table and calculates columns for distance in parsecs and distance in light years

# get_cartesian() - takes an Astropy Table and adds calculated columns for XYZ and UVW if given proper motions and radial velocity

# get_num_nearby() - takes an Astropy Table and calculates the number of objects within the table that are within a specified 
#                    distance to each object

import pandas as pd
from astropy.table import Table
import astropy.coordinates
import astropy.units as u
import numpy as np
import collections

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
def get_cartesian(data:Table, frame='icrs', dist='dist_pc', ra='RAdeg', dec='DEdeg', glon='GLON', glat='GLAT', pmra='pmra', pmde='pmde', pmglon='pmglon', pmglat='pmglat', radial_velocity='radial_velocity'):
    
    #Raise exception if distance is not in data
    if(dist not in data.columns):
        raise Exception('distance must be provided')
    
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
                            description='Position (x coordinate) in parsecs')
    
    data['y'] = data.MaskedColumn(data=galactic_coords.cartesian.y, 
                            meta = collections.OrderedDict([('ucd', 'pos.cartesian.y')]),
                            format='{:.6f}', 
                            description='Position (y coordinate) in parsecs')
    
    data['z'] = data.MaskedColumn(data=galactic_coords.cartesian.z, 
                            meta = collections.OrderedDict([('ucd', 'pos.cartesian.z')]),
                            format='{:.6f}', 
                            description='Position (z coordinate) in parsecs')
    
    
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

#calculates the number of nearby objects for each object in an Astropy Table
#XYZ must already be calculated and be labeled as 'x', 'y', 'z'
def get_num_nearby(data:Table, distfactor:float):
    distfactor=10 # distance factor in parsecs
    n_near = [len(data[np.sqrt((data['x']-data['x'][i])**2 + (data['y']-data['y'][i])**2 + 
                                       (data['z']-data['z'][i])**2) < distfactor]) for i in range(len(data))]
    data['N_near'] = data.MaskedColumn(data=n_near,
                                       dtype=int,
                                       meta = collections.OrderedDict([('ucd', 'meta.number')]),
                                       description='Number of objects in the table within '+str(distfactor)+' parsecs)

           
    