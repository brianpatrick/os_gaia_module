import sys
from pathlib import Path
import csv
from astropy.table import Table
import pandas as pd


# test the input file for existence. Exit if it's not there.
# -----------------------------------------------------------------------------
def test_input_file(path):
    '''
    Test for the existence of an input file,
    Exit if it's not there.
    '''
    if not path.is_file():
        #print(str(path) + ' does not exist.')
        sys.exit('\n\tInput file does not exist:\n\t' + str(path) + '\n\tExiting.\n')



# test a path for potential output file.
# -----------------------------------------------------------------------------
def test_path(path):
    '''
    Test is a directory (path) exists, 
    and create any part of it that does not exist.
    '''
    # Set a relative filepath so we don't include the /home/...
    relative_filepath = path.relative_to(Path.cwd())

    #print(path, relative_filepath)

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




# Write to a speck formatted file
# -----------------------------------------------------------------------------
def to_speck(metadata, df, columns):

    # 'columns' argument is a dataframe which contains the columns of df to be included in the necessary order
    # 'columns' also contains a 'description' column which contains descriptions to be attached to each column

    df = df[columns['name']]
    
    filename = metadata['fileroot'] + '.speck'
    out = open(filename, 'w')


    # Print the metadata defined in the main function
    print('# ' + metadata['project'] + ': ' + metadata['sub_project'], file=out)
    print('# ' + metadata['catalog'] + ', ' + metadata['author'], file=out)
    print('# Prepared by: ' + metadata['prepared_by'], file=out)
    print('# Version ' + metadata['version'], file=out)
    print('# ' + metadata['data_group_desc'], file=out)
    print(file=out)


    # Print the datavar lines.
    # We look at the columns after the x,y,z (start range at 3), 
    # and up until "speck_label" column, which must appear in 
    # every speck file as named.
    datavar_lines = ''
    for i in range(3, list(df.columns).index('speck_label')):
        datavar_lines += 'datavar ' + str(i-3) + ' ' + df.columns[i] + '\t # ' + columns['description'][i] + '\n'
# names
    print(datavar_lines, file=out)

    # Print the data
    # We replace the '__' with a space because we add the '__' for spaces in the names and 
    # speck comment so that 
    print(df.to_csv(path_or_buf=None, sep=' ', na_rep='0', header=False, index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ", lineterminator='\n', float_format='%.8f').replace('__', ' '), file=out)



# Write to a label formatted file
# -----------------------------------------------------------------------------
def to_label(metadata, df):

    # df must have a column called 'label' containing the primary label

    df['text'] = ['text']*len(df)

    df_label = df[['x', 'y', 'z', 'text', 'label']]    

    filename = metadata['fileroot'] + '.label'
    out = open(filename, 'w')

    # Print the metadata defined in the main function
    print('# ' + metadata['project'] + ': ' + metadata['sub_project'], file=out)
    print('# ' + metadata['catalog'] + ', ' + metadata['author'], file=out)
    print('# Prepared by: ' + metadata['prepared_by'], file=out)
    print('# Version ' + metadata['version'], file=out)
    print('# ' + metadata['data_group_desc'], file=out)
    print(file=out)

    # Print the data
    print(df_label.to_csv(path_or_buf=None, sep=' ', na_rep='0', header=False, index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ", lineterminator='\n', float_format='%.8f').replace('__', ' '), file=out)
    
    

# Construct a DataFrame for column metadata
def get_metadata(table:Table, columns:list):
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

