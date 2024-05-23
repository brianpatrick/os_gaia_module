import sys
from pathlib import Path
import csv
import pandas as pd
from astropy.table import Table



# -----------------------------------------------------------------------------
def test_input_file(path):
    """
    Test for the existence of an input file, exit program if it's not there.

    :param path: A python path to the potential input file.
    :type path: pathlib.PosixPath
    :raises FileNotFoundError: Raised if the file does not exist.
    """
    if not path.is_file():
        #sys.exit('\n\tInput file does not exist:\n\t' + str(path) + '\n\tExiting.\n')
        raise FileNotFoundError('input file does not exist:\n\t' + str(path) + '\n\tExiting.\n')





# -----------------------------------------------------------------------------
def test_path(path):
    """
    Test is a directory (path) exists, and create any part of it that does not exist.

    :param path: A python path to a directory where a file will potentially be written.
    :type path: pathlib.PosixPath
    """
    # Set a relative filepath so we don't include the /home/...
    relative_filepath = path.relative_to(Path.cwd())

    # If the path doesn't exist, find a graceful way to exit.
    if not Path.exists(path):
        permission_create_dir = input('\tCreate directory ./' + str(relative_filepath) + '? (y/n/q): ')
        
        if permission_create_dir == 'y':
            Path(path).mkdir(parents=True)
            print('\t -- Created directory: ./' + str(relative_filepath))
        elif permission_create_dir == 'n':
            sys.exit('\n\t -- Cannot write output file. Rerun and create the directory.\n\tExiting.\n')
        elif permission_create_dir == 'q':
            sys.exit("\n\t -- You've chosen to quit.\n\tExiting.\n")
        else:
            sys.exit("\n\t -- Not a valid choice. Choose 'y' to create the necessary directory.\n\tExiting.\n")
    # else:   # debugging purposes
    #     print('Path exists???')





# -----------------------------------------------------------------------------
def write_header(metadata, out):
    """
    Write header for a speck or label file
    
    Write a header for a speck or label file using the metadata about the data

    :param metadata: A dataframe with metadata about the data set.
    :type metadata: DataFrame
    :param out: A file to print out to (needs to be opened before function call)
    :type df: file
    """
    
    # Print the metadata defined in the main function
    print('# ' + metadata['project'] + ': ' + metadata['sub_project'], file=out)
    print('# ' + metadata['catalog'] + ', ' + metadata['catalog_author'] + ', ' + metadata['catalog_year'], file=out)
    print('# Prepared by: ' + metadata['prepared_by'], file=out)
    print('# Version ' + metadata['version'], file=out)
    print('# ' + metadata['data_group_desc'], file=out)
    print('#', file=out)




 # -----------------------------------------------------------------------------
def to_csv(metadata, df, columns):
    """
    Write dataframe to a comma separated-formatted file.

    Write an OpenSpace-ready csv file using the main dataframe, the metadata about the data, and the columns to print.

    :param metadata: A dataframe with metadata about the data set.
    :type metadata: DataFrame
    :param df: A dataframe of the main data set.
    :type df: DataFrame
    :param columns: A dataframe which contains the columns of ``df`` to be included in the necessary order. It also contains a 'description' column which contains descriptions that are attached to each column in the csv file.
    :type columns: DataFrame
    """

    #set our dataframe with our desired columns and rename speck_label column
    df_csv = df[columns['name']].rename(columns={'speck_label':'identifiers'}, inplace=False)
    
    filename = metadata['fileroot'] + '.csv'
    out = open(filename, 'w', encoding='UTF-8')
    
    # Print the metadata defined in the main function
    write_header(metadata, out)

    # Print the column lines.
    # We look at the columns after the x,y,z (start range at 3), 
    # and up until "speck_label" column, which must appear in 
    # every speck file as named.
    column_lines = ''
    for i in range(len(df_csv.columns)):
        column_lines += '# COLUMN ' + ' ' + df_csv.columns[i] + '\t' + columns['description'][i] + '\n'
    # names
    print(column_lines+'#', file=out)

    # Print the data
    print(df_csv.to_csv(path_or_buf=None, index=False), file=out)





 # -----------------------------------------------------------------------------
def to_speck(metadata, df, columns):
    """
    Write dataframe to a speck-formatted file.

    Write an OpenSpace-ready speck file using the main dataframe, the metadata about the data, and the columns to print.

    :param metadata: A dataframe with metadata about the data set.
    :type metadata: DataFrame
    :param df: A dataframe of the main data set.
    :type df: DataFrame
    :param columns: A dataframe which contains the columns of ``df`` to be included in the necessary order. It also contains a 'description' column which contains descriptions that are attached to each column in the speck file.
    :type columns: DataFrame
    """
    # 'columns' argument is a dataframe containing metadata on the columns we want to print to the speck file

    df_speck = df[columns['name']]
    
    filename = metadata['fileroot'] + '.speck'
    out = open(filename, 'w', encoding='UTF-8')
    
    #replace any spaces in label column with double underscores for speck file formatting
    df_speck.loc[0:len(df_speck),'speck_label'].replace(' ', '__', regex=True, inplace=True)


    # Print the metadata defined in the main function
    write_header(metadata, out)
#     print('# ' + metadata['project'] + ': ' + metadata['sub_project'], file=out)
#     print('# ' + metadata['catalog'] + ', ' + metadata['catalog_author'] + ', ' + metadata['catalog_year'], file=out)
#     print('# Prepared by: ' + metadata['prepared_by'], file=out)
#     print('# Version ' + metadata['version'], file=out)
#     print('# ' + metadata['data_group_desc'], file=out)
#     print(file=out)


    # Print the datavar lines.
    # We look at the columns after the x,y,z (start range at 3), 
    # and up until "speck_label" column, which must appear in 
    # every speck file as named.
    datavar_lines = ''
    for i in range(3, list(df_speck.columns).index('speck_label')):
        datavar_lines += 'datavar ' + str(i-3) + ' ' + df_speck.columns[i] + '\t # ' + columns['description'][i] + '\n'
    # names
    print(datavar_lines, file=out)

    # Print the data
    # We replace the '__' with a space because we add the '__' for spaces in the names and 
    # speck comment so that 
    print(df_speck.to_csv(path_or_buf=None, sep=' ', na_rep='0', header=False, index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ", lineterminator='\n', float_format='%.8f').replace('__', ' '), file=out)





# -----------------------------------------------------------------------------
def to_label(metadata, df):
    """
    Write to a label formatted file.

    Write an OpenSpace-ready label file. The incoming ``df`` *must* have a column called 'label'.

    :param metadata: A dataframe with metadata about the data set.
    :type metadata: DataFrame
    :param df: A dataframe of the main data set.
    :type df: DataFrame
    """    
    # df must have a column called 'label' containing the primary label
    # should probably test for the existence of a label column and throw exception if doesn't exist.
    
    if('label' not in df.columns):
        raise Exception('DataFrame must have a label column called \'label\'')

    df['text'] = ['text']*len(df)

    df_label = df.loc[[i for i in range(len(df)) if len(df['label'][i])>0], ['x', 'y', 'z', 'text', 'label']]
    
    # #replace any spaces in label column with double underscores for label file formatting
    # df_label.loc[0:len(df_label), 'label'].replace(' ', '__', regex=True, inplace=True)

    df_label['label'] = ['__'.join(i.split()) for i in df_label['label']]
    # df_label.loc[0:len(df_label), 'label'] = ['__'.join(i.split()) for i in df_label['label']]

    filename = metadata['fileroot'] + '.label'
    out = open(filename, 'w', encoding='UTF-8')

    # Print the metadata defined in the main function
    write_header(metadata, out)
#     print('# ' + metadata['project'] + ': ' + metadata['sub_project'], file=out)
#     print('# ' + metadata['catalog'] + ', ' + metadata['catalog_author'] + ', ' + metadata['catalog_year'], file=out)
#     print('# Prepared by: ' + metadata['prepared_by'], file=out)
#     print('# Version ' + metadata['version'], file=out)
#     print('# ' + metadata['data_group_desc'], file=out)
#     print(file=out)

    # Print the data
    print(df_label.to_csv(path_or_buf=None, sep=' ', na_rep='0', header=False, index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ", lineterminator='\n', float_format='%.8f').replace('__', ' '), file=out)





# -----------------------------------------------------------------------------
def get_metadata(table:Table, columns:list):
    """
    Construct a DataFrame for column metadata.

    Given a list ``columns``, construct a DataFrame of the metadata for each of those columns.

    :param table: An astropy table of the main data set.
    :type table: Table
    :param columns: A list of column names from ``data``.
    :type columns: list of str
    :raises Exception: Raised if one column name from ``columns`` is not in the ``data`` table.
    :return: A DataFrame of metadata for each column.
    :rtype: DataFrame
    """    
    for col in columns:
        if(col not in table.columns):
            raise Exception('file_functions.get_metadata: \'' +col+'\' not found in table')
    
    speck_table = table[columns]
    return pd.DataFrame({
        'name': [x for x in speck_table.columns],
        'unit': [speck_table[x].unit if speck_table[x].unit!=None else '' for x in speck_table.columns],
        'datatype': [speck_table[x].dtype if '<U' not in speck_table[x].dtype.str else 'str' for x in speck_table.columns],
        'width': [speck_table[x].format.split('.')[0][2:] if speck_table[x].format!=None else '' for x in speck_table.columns],
        'precision': [speck_table[x].format.split('.')[1][:-2] if speck_table[x].format!=None else '' for x in speck_table.columns],
        'arraysize': ['*' if '<U' in speck_table[x].dtype.str else '' for x in speck_table.columns],
        'ucd': [list(speck_table[x].meta.items())[0][1] for x in speck_table.columns],
        'description': [speck_table[x].description for x in speck_table.columns]})



# -----------------------------------------------------------------------------
def generate_license_file(metadata):
    """
    This function creates a new file called license_*.lua and 
    prints the license for a particular dataset.

    :param metadata: A dataframe with metadata about the data set.
    :type metadata: DataFrame
    """    

    filename = 'license_' + metadata['fileroot'] + '.lua'
    out = open(filename, 'w')


    DU_license = r"""Copyright © American Museum of Natural History. All rights reserved. 

Downloading the Atlas:
AMNH offers the Atlas free of charge via our website, http://www.haydenplanetarium.org/. The User is free to download and/or duplicate verbatim copies of the Atlas provided this license and copyright information accompany the Atlas.

Modifying the Atlas:
The user is free to modify the Atlas by either adding data or altering existing data, provided it is for personal use only. Once the user modifies the Atlas, it is no longer part of AMNH's Atlas and cannot be used publicly alongside or within the Atlas without written permission from AMNH.

Distributing the Atlas:
The user is forbidden to distribute and use the Atlas for profit, as part of a software and/or production system that will subsequently be redistributed, or for public consumption (via print, electronic media, or broadcast/produced pieces) without written permission from AMNH.

Neither the names of American Museum of Natural History and Hayden Planetarium nor the names of their contributors may be used to endorse or promote products derived from this Atlas without specific, prior written permission.

The Atlas is free but is offered WITHOUT ANY WARRANTY of any kind. We provide the Atlas as is and take no responsibility for any damage resulting from the use of this Atlas. The entire risk as to the quality and performance of this product is with the user.

For more information, please visit https://www.amnh.org/research/hayden-planetarium/digital-universe."""

    # Print the metadata in the license formatted strings
    print('return {', file=out)
    print('    Name = "' + metadata['data_group_title'] + '",', file=out)
    print('    Version = "' + metadata['version'] + '",', file=out)
    print('    Description = "' + metadata['data_group_desc'] + '",', file=out)
    print('    Reference = "' + metadata['catalog'] + ', ' + metadata['catalog_author'] + ', ' + metadata['catalog_year'] + '",', file=out)
    print('    PreparedBy = "' + metadata['prepared_by'] + '",', file=out)
    print('    License = [[' + DU_license + ']]', file=out)
    print('}', file=out)
