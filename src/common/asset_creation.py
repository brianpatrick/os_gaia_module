from pathlib import Path
import sys


# functions to write parts of the asset file
def asset_outpath(datainfo):
    """
    Generate the output path for the asset file.
    :param datainfo: Metadata about the dataset.
    :type datainfo: dict of {str : list}
    :return: The output path for the asset file.
    :rtype: str
    """
    # Get the current working directory
    cwd = Path.cwd()
    
    # Create the output path using the current working directory and the asset name
    try:
        outfile = datainfo['filename'] + '.asset'
        outpath = cwd / datainfo['asset_dir']
    except KeyError:
        raise KeyError(f'need to set asset_dir in the datainfo dictionary')
        
    # Create the directory if it doesn't exist
    outpath.mkdir(parents=True, exist_ok=True)

    # Return the output path
    return outpath / outfile


def import_local_modules(datainfo):
    # importing local modules
    print('-- Set some parameters for OpenSpace settings')
    try:
        print(f'local data_file = asset.resource("{datainfo["data"]["File"]}")')
    except KeyError:
        raise KeyError(f'need to set File in the datainfo["data"] dictionary')

    if datainfo.get("Texture") != None:
        print(f'local texture = asset.resource("{datainfo["Texture"]["File"]}")')
    
    if datainfo.get("ColorMapping") != None:
        print(f'local cmap = asset.resource("{datainfo["ColorMapping"]["File"]}")')

    if datainfo.get("Labels") != None:
        if datainfo["Labels"].get("File") != None:
            print(f'local label = asset.resource("{datainfo["Labels"]["File"]}")')
    
    if datainfo["data"].get("NamesFile") != None:
        print(f'local names = asset.resource("{datainfo["data"]["NamesFile"]}")')

    if datainfo.get("data_folder") != None: # used for star orbits to iterate through multiple speck files
        print(f'local folder = asset.resource("{datainfo["data_folder"]}")')

    print()

def import_OpenSpace_modules(datainfo):
    # importing local modules locations
    print('-- Set some parameters for OpenSpace settings')
    try:
        print('local data_folder = asset.resource({')
        print(f'  Name = "{datainfo["data"]["Name"]}",')
        print( '  Type = "HttpSynchronization",')
        print(f'  Identifier = "{datainfo["data"]["Identifier"]}",')
        print(f'  Version = {datainfo["data"]["Version"]}')
        print('})')
    except KeyError:
        raise KeyError(f'need to set Name, Version, and Identifier \
                        in the datainfo["data"] dictionary when local_modules = False')

    if datainfo.get("Texture") != None: # if texture is needed its found in this folder
        try:
            print('local texture_folder = asset.resource({')
            print(f'  Name = "{datainfo["Texture"]["Name"]}",')
            print( '  Type = "HttpSynchronization",')
            print(f'  Identifier = "{datainfo["Texture"]["Identifier"]}",')
            print(f'  Version = {datainfo["Texture"]["Version"]}')
            print('})')
        except KeyError:
            raise KeyError(f'need to set Name, Identifier, and Version \
                           in the datainfo["Texture"] dictionary when local_modules = False \
                           and texture is provided')
        
    if datainfo.get("ColorMapping") != None: # if texture is needed its found in this folder
        try:
            print('local cmap_folder = asset.resource({ -- often the same as the data_folder but not always')
            print(f'  Name = "{datainfo["ColorMapping"]["Name"]}",')
            print( '  Type = "HttpSynchronization",')
            print(f'  Identifier = "{datainfo["ColorMapping"]["Identifier"]}",')
            print(f'  Version = {datainfo["ColorMapping"]["Version"]}')
            print('})')
        except KeyError:
            raise KeyError(f'need to set Name, Identifier, and Version \
                           in the datainfo["ColorMapping"] dictionary when local_modules = False \
                           and texture is provided')
        
    # assign variables for the OpenSpace modules
    try:
        print(f'local data_file = data_folder .. "{datainfo["data"]["File"]}"')
    except KeyError:
        raise KeyError(f'need to set File in the datainfo["data"] dictionary')

    if datainfo.get("Texture") != None:
        print(f'local texture = texture_folder .. "{datainfo["Texture"]["File"]}"')
    
    if datainfo.get("ColorMapping") != None:
        print(f'local cmap = cmap_folder .. "{datainfo["ColorMapping"]["File"]}"')

    if datainfo.get("Labels") != None:
        if datainfo["Labels"].get("File") != None:
            print(f'local label = data_folder .. "{datainfo["Labels"]["File"]}"')
    
    if datainfo["data"].get("NamesFile") != None:
        print(f'local names = data_folder .. "{datainfo["data"]["NamesFile"]}"')
    
    print()

def import_modules(datainfo):
    if datainfo["local_modules"] == True:
        import_local_modules(datainfo)
    else:
        import_OpenSpace_modules(datainfo)

def import_local_star_modules(datainfo):
        # importing local modules
        print('-- Set some parameters for OpenSpace settings')
        try:
            print(f'local data_file = asset.resource("{datainfo["data"]["File"]}")')
            print(f'local core_texture = asset.resource("{datainfo["Texture"]["Core"]}")')
            print(f'local glare_texture = asset.resource("{datainfo["Texture"]["Glare"]}")')
            print(f'local cmap = asset.resource("{datainfo["ColorMap"]["ColorMap"]}")')
            print(f'local other_cmap = asset.resource("{datainfo["ColorMap"]["OtherDataColorMap"]}")')
        except KeyError:
            raise KeyError(f'need to set File in datainfo["data"], Core and Glare in datainfo["Texture"], and ColorMap and OtherDataColorMap in datainfo')

        print()

def import_OpenSpace_star_modules(datainfo):
    # importing local modules locations
    print('-- Set some parameters for OpenSpace settings')
    try:
        print('local data_folder = asset.resource({')
        print(f'  Name = "{datainfo["data"]["Name"]}",')
        print( '  Type = "HttpSynchronization",')
        print(f'  Identifier = "{datainfo["data"]["Identifier"]}",')
        print(f'  Version = {datainfo["data"]["Version"]}')
        print('})')
        print()
    except KeyError:
        raise KeyError(f'need to set Name, Version, and Identifier \
                        in the datainfo["data"] dictionary when local_modules = False')

    try:
        print('local texture_folder = asset.resource({')
        print(f'  Name = "{datainfo["Texture"]["Name"]}",')
        print( '  Type = "HttpSynchronization",')
        print(f'  Identifier = "{datainfo["Texture"]["Identifier"]}",')
        print(f'  Version = {datainfo["Texture"]["Version"]}')
        print('})')
        print()
    except KeyError:
        raise KeyError(f'need to set Name, Identifier, and Version \
                        in the datainfo["Texture"] dictionary when local_modules = False')
    
    try:
        print('local colormaps_folder = asset.resource({')
        print(f'  Name = "{datainfo["ColorMap"]["Name"]}",')
        print( '  Type = "HttpSynchronization",')
        print(f'  Identifier = "{datainfo["ColorMap"]["Identifier"]}",')
        print(f'  Version = {datainfo["ColorMap"]["Version"]}')
        print('})')
        print()
    except KeyError:
        raise KeyError(f'need to set Name, Identifier, and Version \
                        in the datainfo["ColorMap"] dictionary when local_modules = False')
        
    # assign variables for the OpenSpace modules
    try:
        print(f'local data_file = data_folder .. "{datainfo["data"]["File"]}"')
        print(f'local core_texture = texture_folder .. "{datainfo["Texture"]["Core"]}"')
        print(f'local glare_texture = texture_folder .. "{datainfo["Texture"]["Glare"]}"')
        print(f'local cmap = colormaps_folder ..  "{datainfo["ColorMap"]["ColorMap"]}"')
        print(f'local other_cmap =  colormaps_folder .. "{datainfo["ColorMap"]["OtherDataColorMap"]}"')
    except KeyError:
        raise KeyError(f'need to set File in datainfo["data"], Core and Glare in datainfo["Texture"], and ColorMap and OtherDataColorMap in the datainfo dictionary')

    
    print()

def import_star_modules(datainfo):
    if datainfo["local_modules"] == True:
        import_local_star_modules(datainfo)
    else:
        import_OpenSpace_star_modules(datainfo)

def set_color_parameters(datainfo):
    # if color map is provided 
    if datainfo.get("ColorMapping") != None:
        print('    Coloring = {')
        # provided fix color if applicable
        if datainfo.get("FixedColor") != None:
            print('      FixedColor = {' + datainfo["FixedColor"] +'},')
        print('      ColorMapping = {')
        print('        File = cmap,')
        print('        Enabled = ' + datainfo["ColorMapping"]["Enabled"] + ',')
        print('        ParameterOptions = {')
        
        ## input color map parameters
        for index, param in enumerate(datainfo["ColorMapping"]["ParameterOptions"]):
            # not the last parameter, add a comma
            if index != len(datainfo["ColorMapping"]["ParameterOptions"]) - 1:
                print('          { Key = "' + param["Key"] + '", Range = {' + param["Range"] + '} },')
            # last parameter, no comma
            else:
                print('          { Key = "' + param["Key"] + '", Range = {' + param["Range"] + '} }')
        print('        }')
        print('      }')
        print('    },')
    # if color map is not provided
    elif datainfo.get("FixedColor") != None:
        print('    Coloring = {')
        print('      FixedColor = {' + datainfo["FixedColor"] +'}')
        print('    },')
    
    # set opacity if provided
    if datainfo.get("Opacity") != None:
        print('    Opacity = ' + datainfo["Opacity"] + ',') 

def set_size_parameters(datainfo):
    try:
        print('    SizeSettings = {') 
        if datainfo.get("SizeMapping_ParameterOptions") != None:
            print('      SizeMapping = {')
            print('        ParameterOptions = { "' + datainfo["SizeMapping_ParameterOptions"] + '" }')
            print('        },')
        print('      ScaleExponent = ' + datainfo["ScaleExponent"] + ',') 
        print('      MaxSize = ' + datainfo["MaxSize"] + ',')
        print('      EnableMaxSizeControl = true')
        print('     },')

    except KeyError:
        raise KeyError(f'need to set ScaleExponent and MaxSize in the datainfo dictionary')

def set_fading_parameters(datainfo):
    if datainfo.get("FadeInDistances") != None:
        print('    Fading = {')
        print('      FadeInDistances = {' + datainfo["FadeInDistances"] + '}')
        print('    },')

def set_label_parameters(datainfo):
    
    if datainfo.get("Labels") != None:
        ## input label settings if csv_labels = True
        if datainfo["Labels"].get("csv") == True:

            # check if label column name is provided
            if datainfo["Labels"].get("column") != None:
                label_column = datainfo["Labels"]["column"]
            else:
                label_column = "label"
            
            print('    DataMapping = {')
            print('      Name = "' + label_column + '" },')

        if datainfo["Labels"].get("csv") == True or datainfo["Labels"].get("File") != None:
            # check used in label file or csv
            try:
                print('    Labels = {')
            
                if datainfo["Labels"].get("File") != None: # if label file is provided
                    print('      File = label,')

                if datainfo["Labels"].get("Enabled") == True:
                    print('      Enabled = true,')

                print('      Color = { ' + datainfo["Labels"]["Color"] + ' },')
                print('      Size = ' + datainfo["Labels"]["Size"] + ',')
                print('      MinMaxSize = { ' + datainfo["Labels"]["MinMaxSize"] + ' },')
                print('      Unit = "' + datainfo["Unit"] + '"')
                print('    }')
            except KeyError:
                raise KeyError(f'need to set Color, Size, and MinMaxSize in the datainfo["Labels"] dictionary')
            
    

def define_GUI(datainfo):
    """
    Define the GUI settings for the asset file.
    :param datainfo: Metadata about the dataset.
    :type datainfo: dict of {str : list}
    """
    # GUI settings
    try:
        print('  GUI = {')
        print('    Name = "' + datainfo['GUI']['Name'] + '",')
        print('    Path = "' + datainfo['GUI']['Path'] + '",')
        print('    Description = [[' + datainfo['GUI']['Description'] + ']]')
        print('  }') # end of GUI
    except KeyError:
        raise KeyError(f'need to set Name, Path, and Description in the datainfo["GUI"] dictionary')
    
def asset_metadata(datainfo):

    try:
        print('asset.meta = {')
        print('  Name = "' + datainfo["meta_name"] + '",')
        print('  Description = Object.GUI.Description,')
        print('  Author = "' + datainfo["author"] + '",')
        print('  URL = "https://www.amnh.org/research/hayden-planetarium/digital-universe",')
        print('  License = "AMNH Digital Universe"')
        print('}')
    
    except KeyError:
        raise KeyError(f'need to set meta_name and author in the datainfo dictionary')

def initialize_asset_functions(Object = "Object"):
    """
    Initialize the asset functions for the asset file.
    :param datainfo: Metadata about the dataset.
    :type datainfo: dict of {str : list}
    """
    print('asset.onInitialize(function()')
    print('    openspace.addSceneGraphNode(' + Object + ')')
    print('end)')
    print()
    print('asset.onDeinitialize(function()')
    print('    openspace.removeSceneGraphNode(' + Object + ')')
    print('end)')
    print()
    print('asset.export(' + Object + ')')
    print()
    print()

def RenderableConstellationLines_Object(datainfo):
    # START OF OBJECT 
    print('local Object = {')
    print('  Identifier = "' + datainfo["Identifier"] + '",')
    # START OF RENDERABLE
    print('  Renderable = {')
    print('    Type = "RenderableConstellationLines",')
    # input file settings
    print('    File = data_file,')
    print('    NamesFile = names,')
    print('    Unit = "' + datainfo["Unit"] + '",')
    print('    Enabled = false,')
    # color, opacity, texture, size settings
    if datainfo.get("LineWidth") != None:
        print('    LineWidth= ' + datainfo["LineWidth"] + ',')
    print('    Colors = { ' + datainfo["Colors"] + ' },')
    if datainfo.get("Opacity") != None:
        print('    Opacity = ' + datainfo["Opacity"] + ',') 
    print('    DimInAtmosphere = true,')

    set_label_parameters(datainfo)
    print('  },') # END OF RENDERABLE

    print('  Tag = { "daytime_hidden" },')
    # Print GUI settings
    define_GUI(datainfo)
    print('}') # END OF OBJECT

def RenderableDUMeshes_Object(datainfo, object):

    object_dict = datainfo["objects"][object]
    
    # START OF OBJECT 
    print('local ' + object + ' = {')    
    
    # variables from object dictionary
    try:
        print('  Identifier = "' + object_dict["Identifier"] + '",')
        # START OF RENDERABLE
        print('  Renderable = {')
        print('    Type = "RenderableDUMeshes",')
        print('    File = folder .. "' + object_dict['speck_file'] + '",')
        print('    MeshColor = { {' + object_dict['MeshColor'] + '} },')
    except KeyError:
        raise KeyError(f'need to set identifier, speck_file, and MeshColor in the datainfo["objects"]["{object}"] dictionary')

    # variables in main dictionary
    try:
        print('    Opacity = ' + datainfo['Opacity'] + ',')
        print('    Enabled = ' + datainfo['Enabled'] + ',')
        print('    Unit = "' + datainfo["Unit"] + '",')
    except KeyError:
        raise KeyError(f'need to set Opacity, Enabled, and Unit in the datainfo dictionary')
    
    print('  },') # END OF RENDERABLE

    try:
        define_GUI(object_dict)
    except KeyError:
        raise KeyError(f'need to set Name, Path, and Description in the datainfo["objects"]["{object}"]["GUI"] dictionary')
    
    print('}') # END OF OBJECT
    print()

# function to create the asset file for various renderables
def make_RenderableConstellationLines_asset(datainfo):
        # We shift the stdout to our filehandle so that we don't have to keep putting
    # the filehandle in every print statement.
    ## Save the original stdout so we can switch back later
    original_stdout = sys.stdout
    ## Open the file to write to
    outpath = asset_outpath(datainfo)

    with open(outpath, 'wt') as asset:
        # Switch stdout to the file
        sys.stdout = asset
        
        # Print the header
        print("-- This file is auto-generated in the " + make_RenderablePointCloud_asset.__name__ + "() function inside " + Path(__file__).name)
        print()

        import_modules(datainfo)
        
        #  create zodiac objects if applicable
        if datainfo.get("zodiacs") == True:
            print('local constellations_helper = asset.require("util/constellations_helper") ')
            
            print('''local zodiacs = {
"CNC", "TAU", "PSC", "ARI", "LIB", "AQR", "CAP", "SCO", "VIR", "SGR", "GEM", "LEO"
}

local function zodiacsString(zodiacsList)
local zodiacsString = "{"
local isFirst = true

for k, zodiac in pairs(zodiacsList) do
    local fullName = constellations_helper.findFullName(zodiac)
    if fullName ~= nil then
    if isFirst then
        isFirst = false
    else
        zodiacsString = zodiacsString .. ", "
    end

    zodiacsString = zodiacsString .. [["]] .. fullName .. [["]]
    end
end

zodiacsString = zodiacsString .. "}"
return zodiacsString
end''') 
            print()
            
        # START OF OBJECT
        RenderableConstellationLines_Object(datainfo)
        
        # initialize asset functions 
        initialize_asset_functions()
        
        if datainfo.get("zodiacs") == True:
            print('local zodiacsString = zodiacsString(zodiacs)')
            print('''local ShowZodiacs = {
  Identifier = "os.constellation.ShowZodiacsAlly",
  Name = "Show zodiac",
  Command = [[
  openspace.setPropertyValueSingle("Scene.Constellations.Renderable.ConstellationSelection", ]] .. zodiacsString .. [[)
  openspace.fadeIn("Scene.Constellations.Renderable")
  ]],
  Documentation = "Shows the zodiac constellations lines",
  GuiPath = "/Constellations/Lines",
  IsLocal = false
  }''')
            # initialize asset functions
            print('asset.onInitialize(function()')
            print('  openspace.action.registerAction(ShowZodiacs)')
            print('end)')

            print('asset.onDeinitialize(function()')
            print('  openspace.action.removeAction(ShowZodiacs)')
            print('end)')
        
        if datainfo.get("constellation_actions") == True:
            print('''-- Actions
local ShowConstellations = {
  Identifier = "os.constellations.ShowConstellationsAlly",
  Name = "Show all",
  Command = [[
      openspace.setPropertyValueSingle("Scene.Constellations.Renderable.ConstellationSelection", {})
      openspace.fadeIn("Scene.Constellations.Renderable")
  ]],
  Documentation = "Shows all the constellations lines",
  GuiPath = "/Constellations/Lines",
  IsLocal = false
  }

local HideConstellations = {
  Identifier = "os.constellations.HideConstellationsAlly",
  Name = "Hide all",
  Command = [[
    openspace.fadeOut("Scene.Constellations.Renderable", nil, "openspace.setPropertyValueSingle('Scene.Constellations.Renderable.Enabled', false); openspace.setPropertyValueSingle('Scene.Constellations.Renderable.ConstellationSelection', {})")
  ]],
  Documentation = "Hides all the constellations lines",
  GuiPath = "/Constellations/Lines",
  IsLocal = false
}''')
            # initialize asset functions
            print('asset.onInitialize(function()')
            print('  openspace.action.registerAction(ShowConstellations)')
            print('  openspace.action.registerAction(HideConstellations)')
            print('end)')
            print()

            print('asset.onDeinitialize(function()')
            print('  openspace.action.removeAction(ShowConstellations)')
            print('  openspace.action.removeAction(HideConstellations)')
            print('end)')
            print()
        
        
        # Print the metadata for the asset
        asset_metadata(datainfo)
    
    # Close the file and switch back to original stdout
    sys.stdout = original_stdout

    # Report to stdout
    print()

def make_RenderablePointCloud_asset(datainfo):
    """
    Generate the asset file for the deep sky survey data.
    :param datainfo: Metadata about the dataset.
    :type datainfo: dict of {str : list}
    Output files:
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    :file:`[{order}]/[{version}]/consensus_species.asset`
        The asset file containing the OpenSpace configurations for the consensus species.
    """
    # We shift the stdout to our filehandle so that we don't have to keep putting
    # the filehandle in every print statement.
    ## Save the original stdout so we can switch back later
    original_stdout = sys.stdout
    ## Open the file to write to
    outpath = asset_outpath(datainfo)

    with open(outpath, 'wt') as asset:
        # Switch stdout to the file
        sys.stdout = asset
        
        # Print the header
        print("-- This file is auto-generated in the " + make_RenderablePointCloud_asset.__name__ + "() function inside " + Path(__file__).name)
        print()
        
        import_modules(datainfo)

        # START OF OBJECT 
        print('local Object = {')
        print('  Identifier = "' + datainfo["Identifier"] + '",')
        # START OF RENDERABLE
        print('  Renderable = {')
        print('    Type = "RenderablePointCloud",')
        
        # input file settings
        print('    File = data_file,')
        print('    Unit = "' + datainfo["Unit"] + '",')
        print('    Enabled = ' + datainfo["Enabled"] + ',')
        
        # color, opacity, texture, size settings
        set_color_parameters(datainfo)
        ## input texture settings
        print('    Texture = {')
        print('      File = texture },')
        ## size settings
        set_size_parameters(datainfo)
        ## fading settings
        set_fading_parameters(datainfo)
        ## input label settings if csv_labels = True
        set_label_parameters(datainfo)
        
        print('  },') # END OF RENDERABLE
        
        # input GUI settings
        define_GUI(datainfo)

        print('}') # END OF OBJECT
        print()

        # initialize asset functions 
        initialize_asset_functions()

        # Print the metadata for the asset
        asset_metadata(datainfo)

    # Close the file and switch back to original stdout
    sys.stdout = original_stdout

    # Report to stdout
    print()

def make_star_labels_asset(datainfo):
    """
    Generate the asset file for the deep sky survey data.
    :param datainfo: Metadata about the dataset.
    :type datainfo: dict of {str : list}
    Output files:
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    :file:`[{order}]/[{version}]/consensus_species.asset`
        The asset file containing the OpenSpace configurations for the consensus species.
    """
    # We shift the stdout to our filehandle so that we don't have to keep putting
    # the filehandle in every print statement.
    ## Save the original stdout so we can switch back later
    original_stdout = sys.stdout
    ## Open the file to write to
    outpath = asset_outpath(datainfo)

    # set label enabled to true
    datainfo["Labels"]["Enabled"] = True

    with open(outpath, 'wt') as asset:
        # Switch stdout to the file
        sys.stdout = asset

        # Print the header
        print("-- This file is auto-generated in the " + make_star_labels_asset.__name__ + "() function inside " + Path(__file__).name)
        print()
        
        import_modules(datainfo)

        # START OF OBJECT 
        print('local Object = {')
        print('  Identifier = "' + datainfo["Identifier"] + '",')
        # START OF RENDERABLE
        print('  Renderable = {')
        print('    File = data_file,')
        print('    Type = "RenderablePointCloud",')
        print('    Enabled = ' + datainfo["Enabled"] + ',')
        
        # opacity, units, size settings
        print('    Opacity = ' + datainfo["Opacity"] + ',')
        print('    Unit = "' + datainfo["Unit"] + '",')
        print('    SizeSettings = {') 
        print('      ScaleExponent = 0,') 
        print('      MaxSize = 0,')
        print('      EnableMaxSizeControl = false')
        print('     },')

        ## input label settings if csv_labels = True
        set_label_parameters(datainfo)
        
        print('  },') # END OF RENDERABLE
        
        # input GUI settings
        define_GUI(datainfo)

        print('}') # END OF OBJECT
        print()

        # initialize asset functions 
        initialize_asset_functions()

        # Print the metadata for the asset
        asset_metadata(datainfo)
    
    # Close the file and switch back to original stdout
    sys.stdout = original_stdout

    # Report to stdout
    print()

def make_RenderablePolygonCloud_asset(datainfo):
    """
    Generate the asset file for the deep sky survey data.
    :param datainfo: Metadata about the dataset.
    :type datainfo: dict of {str : list}
    Output files:
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    :file:`[{order}]/[{version}]/consensus_species.asset`
        The asset file containing the OpenSpace configurations for the consensus species.
    """
    # We shift the stdout to our filehandle so that we don't have to keep putting
    # the filehandle in every print statement.
    ## Save the original stdout so we can switch back later
    original_stdout = sys.stdout
    ## Open the file to write to
    outpath = asset_outpath(datainfo)

    with open(outpath, 'wt') as asset:
        # Switch stdout to the file
        sys.stdout = asset
        
        # Print the header
        print("-- This file is auto-generated in the " + make_RenderablePolygonCloud_asset.__name__ + "() function inside " + Path(__file__).name)
        print()
        
        import_modules(datainfo)

        # START OF OBJECT 
        print('local Object = {')
        print('  Identifier = "' + datainfo["Identifier"] + '",')
        # START OF RENDERABLE
        print('  Renderable = {')
        print('    Type = "RenderablePolygonCloud",')
        print('    Enabled = ' + datainfo["Enabled"] + ',')
        
        # input file settings
        print('    File = data_file,')
        print('    Unit = "' + datainfo["Unit"] + '",')
        print('    PolygonSides = ' + datainfo["PolygonSides"] + ',')
        
        # color, opacity, texture, size settings
        set_color_parameters(datainfo)

        ## size settings
        set_size_parameters(datainfo)
        ## fading settings
        set_fading_parameters(datainfo)
        ## input label settings if csv_labels = True
        set_label_parameters(datainfo)
        
        print('  },') # END OF RENDERABLE
        
        # input GUI settings
        define_GUI(datainfo)

        print('}') # END OF OBJECT
        print()

        # initialize asset functions 
        initialize_asset_functions()

        # Print the metadata for the asset
        asset_metadata(datainfo)
    
    # Close the file and switch back to original stdout
    sys.stdout = original_stdout

    # Report to stdout
    print()

def make_RenderableStars_asset(datainfo):
    """
    Generate the asset file for the deep sky survey data.
    :param datainfo: Metadata about the dataset.
    :type datainfo: dict of {str : list}
    Output files:
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    :file:`[{order}]/[{version}]/consensus_species.asset`
        The asset file containing the OpenSpace configurations for the consensus species.
    """
    # We shift the stdout to our filehandle so that we don't have to keep putting
    # the filehandle in every print statement.
    ## Save the original stdout so we can switch back later
    original_stdout = sys.stdout
    ## Open the file to write to
    outpath = asset_outpath(datainfo)

    with open(outpath, 'wt') as asset:
        # Switch stdout to the file
        sys.stdout = asset
        
        # Print the header
        print("-- This file is auto-generated in the " + make_RenderableStars_asset.__name__ + "() function inside " + Path(__file__).name)
        print()
        
        import_star_modules(datainfo)

        # START OF OBJECT 
        print('local Object = {')
        print('  Identifier = "' + datainfo["Identifier"] + '",')
        # START OF RENDERABLE
        print('  Renderable = {')
        print('    Type = "RenderableStars",')

        # input file settings
        print('    File = data_file,')
        
        # input glare and halo settings
        print('    Core = {')
        print('      Texture = core_texture,')
        print('      Multiplier = 15.0,')
        print('      Gamma = 1.66,')
        print('      Scale = 0.18')
        print('      },')
        print('    Glare = {')
        print('      Texture = glare_texture,')
        print('      Multiplier = 0.65')
        print('      },')
        
        ## misc color size settings
        print('    MagnitudeExponent = 6.325,')
        print('    ColorMap = cmap,')
        print('    OtherDataColorMap = other_cmap,')
        print('    SizeComposition = "Distance Modulus",')
        
        print('    DataMapping = {')
        print('      Bv = "' + datainfo['Bv_column']+'",')
        print('      Luminance = "' + datainfo['Luminance_column']+'",')
        print('      AbsoluteMagnitude = "' + datainfo['AbsoluteMagnitude_column']+'",')
        print('      ApparentMagnitude = "' + datainfo['ApparentMagnitude_column']+'",')
        print('      Vx = "' + datainfo['Vx_column']+'",')
        print('      Vy = "' + datainfo['Vy_column']+'",')
        print('      Vz = "' + datainfo['Vz_column']+'",')
        print('      Speed = "' + datainfo['Speed_column']+'",')
        print('    },')
        
        print('    DimInAtmosphere = true')
        
        print('  },') # END OF RENDERABLE
        
        print('  Tag = { "daytime_hidden" },')

        # input GUI settings
        define_GUI(datainfo)

        print('}') # END OF OBJECT
        print()

        # initialize asset functions 
        initialize_asset_functions()

        # Print the metadata for the asset
        asset_metadata(datainfo)
    
    # Close the file and switch back to original stdout
    sys.stdout = original_stdout

    # Report to stdout
    print()

def make_RenderableDUMeshes_asset(datainfo):
     # We shift the stdout to our filehandle so that we don't have to keep putting
    # the filehandle in every print statement.
    ## Save the original stdout so we can switch back later

    ## this function is used to create an asset for star orbits which creates multiple objects 

    original_stdout = sys.stdout
    ## Open the file to write to
    outpath = asset_outpath(datainfo)

    with open(outpath, 'wt') as asset:
        # Switch stdout to the file
        sys.stdout = asset
        
        # Print the header
        print("-- This file is auto-generated in the " + make_RenderableDUMeshes_asset.__name__ + "() function inside " + Path(__file__).name)
        print()
        
        if datainfo["local_modules"] == True:
            try:
                print(f'local folder = asset.resource("{datainfo["data_folder"]}")')
            except KeyError:
                raise KeyError(f'need to set data_folder in the datainfo dictionary when local_modules = True')
        else:
            try:
                print('local folder = asset.resource({')
                print(f'  Name = "{datainfo["data"]["Name"]}",')
                print( '  Type = "HttpSynchronization",')
                print(f'  Identifier = "{datainfo["data"]["Identifier"]}",')
                print(f'  Version = {datainfo["data"]["Version"]}')
                print('})')
            except KeyError:
                raise KeyError(f'need to set data_folder_Name, data_folder_Version, and data_folder_Identifier \
                                in the datainfo dictionary when local_modules = False')


        for object in datainfo["objects"].keys():
            
            RenderableDUMeshes_Object(datainfo, object)

        # initialize asset functions 
        print('asset.onInitialize(function()')
        for object in datainfo["objects"].keys():
            print('    openspace.addSceneGraphNode(' + object + ')')
        print('end)')
        print()

        # deinitialize asset functions
        print('asset.onDeinitialize(function()')
        for object in datainfo["objects"].keys():
            print('    openspace.removeSceneGraphNode(' + object + ')')
        print('end)')
        print()

        # export asset functions
        for object in datainfo["objects"].keys():
            print('asset.export(' + object + ')')
        print()
        print()

        # Print the metadata for the asset
        try:
            print('asset.meta = {')
            print('  Name = "' + datainfo["meta_name"] + '",')
            print('  Description = [[' + datainfo["meta_description"] + ']],')
            print('  Author = "' + datainfo["author"] + '",')
            print('  URL = "https://www.amnh.org/research/hayden-planetarium/digital-universe",')
            print('  License = "AMNH Digital Universe"')
            print('}')
        
        except KeyError:
            raise KeyError(f'need to set meta_name, author, and license in the datainfo dictionary')
    
    # Close the file and switch back to original stdout
    sys.stdout = original_stdout

    # Report to stdout
    print()

def make_RenderableConstellationBounds_asset(datainfo):
    # We shift the stdout to our filehandle so that we don't have to keep putting
    # the filehandle in every print statement.
    ## Save the original stdout so we can switch back later
    original_stdout = sys.stdout
    ## Open the file to write to
    outpath = asset_outpath(datainfo)

    with open(outpath, 'wt') as asset:
        # Switch stdout to the file
        sys.stdout = asset
        
        # Print the header
        print("-- This file is auto-generated in the " + make_RenderableConstellationBounds_asset.__name__ + "() function inside " + Path(__file__).name)
        print()
        
        
        import_modules(datainfo)
        print('local coreKernels = asset.require("spice/core")')
        print()

        # START OF OBJECT 
        print('local Object = {')
        print('  Identifier = "' + datainfo["Identifier"] + '",')
        # START OF RENDERABLE
        print('  Renderable = {')
        print('    Type = "RenderableConstellationBounds",')
        print('    Enabled = ' + datainfo["Enabled"] + ',')
        
        # input file settings
        print('    File = data_file,')
        print('    NamesFile = names,')
        
        print('  },') # END OF RENDERABLE
        
        # START OF TRANSFORM
        print('  Transform = {')
        print('    Rotation = {') # START OF ROTATION
        print('      Type = "SpiceRotation",')
        print('      SourceFrame = coreKernels.Frame.J2000,')
        print('      DestinationFrame = coreKernels.Frame.Galactic')
        print('    },') # END OF ROTATION
        print('    Scale = {') # START OF SCALE
        print('      Type = "StaticScale",')
        print('      Scale = 10e17')
        print('    }') # END OF SCALE
        print('  },') # END OF TRANSFORM

    
        # input GUI settings
        define_GUI(datainfo)

        print('}') # END OF OBJECT
        print()

        # initialize asset functions 
        initialize_asset_functions()

        # Print the metadata for the asset
        asset_metadata(datainfo)
    
    # Close the file and switch back to original stdout
    sys.stdout = original_stdout

    # Report to stdout
    print()

# master function to write the asset file
# This function will call the appropriate asset creation function based on the datainfo provided.
def write_asset(datainfo):
    """
    Write the asset file for the deep sky survey data.
    :param datainfo: Metadata about the dataset.
    :type datainfo: dict of {str : list}
    """

    if datainfo["renderable"] == "RenderablePointCloud":
        make_RenderablePointCloud_asset(datainfo)   
    
    elif datainfo["renderable"] == "StarLabels":
        make_star_labels_asset(datainfo)
    
    elif datainfo["renderable"] == "RenderablePolygonCloud":
        make_RenderablePolygonCloud_asset(datainfo)
    
    elif datainfo["renderable"] == "RenderableStars":
        make_RenderableStars_asset(datainfo)
    
    elif datainfo["renderable"] == "RenderableConstellationLines":
        make_RenderableConstellationLines_asset(datainfo)

    elif datainfo["renderable"] == "RenderableDUMeshes":
        make_RenderableDUMeshes_asset(datainfo)
    
    elif datainfo["renderable"] == "RenderableConstellationBounds":
        make_RenderableConstellationBounds_asset(datainfo)
    else:
        raise ValueError(f'Unknown renderable type: {datainfo["renderable"]} \n \
                         Possible values are: RenderablePointCloud, StarLabels, RenderablePolygonCloud, RenderableStars, RenderableConstellationLines, RenderableDUMeshes, RenderableConstellationBounds')
