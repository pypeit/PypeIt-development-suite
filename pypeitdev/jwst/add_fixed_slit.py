import numpy as np

from astropy.io import fits
from astropy.table import Table


def add_fs_target(visit_rate_file, slit_name, ra, dec, stellarity=1.0, new_source_id=None):
    # Read in a rate file to get the msa_meta_id for that visit
    rate_file_header = fits.getheader(visit_rate_file)
    msa_meta_id = rate_file_header['MSAMETID']
    program = rate_file_header['PROGRAM']

    # MOS 3 shutter nod offset size
    nod_offset = 0.529 # arcsec

    # S200A1/A2 approximate slit length
    slit_length = 3.3 # arcsec

    # For a three point dither, need to add three rows, one for each dither point
    dithers = [1, 2, 3]

    # Calculate the x/y positions of the source in the slit (assuming perfectly centered initially)
    initial_x = 0.5
    initial_y = 0.5

    # Assuming no spectral (x) offsets
    xpos = [initial_x, initial_x, initial_x]
    ypos = [initial_y, initial_y + nod_offset/slit_length, initial_y - nod_offset/slit_length]

    # If a source id is provided we don't need to make a new one
    if new_source_id is None:
        new_source_id = np.max(source_table['source_id']) + 1
        
    # Make a new slitlet_id
    new_slit_id = np.max(shutter_table['slitlet_id']) + 1

    # Set up the initial structure that we will fill for each dither position
    shutter_keys = ('dither_point_index', 'msa_metadata_id', 'background',
                    'estimated_source_in_shutter_x',
                    'estimated_source_in_shutter_y',
                    'primary_source', 'shutter_column', 'shutter_quadrant',
                    'shutter_row', 'shutter_state', 'slitlet_id',
                    'source_id', 'fixed_slit')
    shutter_values = (1, msa_meta_id, 'N', 0.5, 0.5, 'Y', 0, 0, 0, 'OPEN',
                    new_slit_id, new_source_id, slit_name)
    source_keys = ('program', 'source_id', 'source_name', 'alias',
                'ra', 'dec', 'preimage_id', 'stellarity')
    source_values = (program, new_source_id, f'{program}_{new_source_id}',
                    f'{new_source_id}', ra, dec, 'None', stellarity)


    for i, dither in enumerate(dithers):
        new_entry = dict(zip(shutter_keys, shutter_values))
        new_entry['dither_point_index'] = dither
        new_entry['estimated_source_in_shutter_x'] = xpos[i]
        new_entry['estimated_source_in_shutter_y'] = ypos[i]
        new_entry['fixed_slit'] = slit_name
        new_entry['slitlet_id'] = new_slit_id
        new_entry['source_id'] = new_source_id

        shutter_table.add_row(new_entry)

        if i == 0 and new_source_id not in source_table['source_id']:
            new_entry = dict(zip(source_keys, source_values))
            new_entry['source_id'] = new_source_id
            new_entry['source_name'] = f"{program}_{new_source_id}"
            new_entry['alias'] = f'{new_source_id}'

            source_table.add_row(new_entry)

    return new_slit_id, new_source_id, shutter_table, source_table

msa_file = 'jw01219006001_01_msa.fits'

# Provide one of the rate files for the visit that you want to add
# Each visit may have a separate msa_meta_id.  You may need to loop through
# an example rate file from each visit and the names of the fixed slits 
# that any targets are in.
visit_rate_file = 'stage1/jw01219006001_06101_00001_nrs1_rate.fits'

# Which slit is the target in?
slit_name = 'S200A1'

# Source RA/Dec
ra = 0.0
dec = 0.0

# point source = 1, extended = 0
stellarity = 1.0

msa_hdu_list = fits.open(msa_file)
shutter_table = Table(msa_hdu_list['SHUTTER_INFO'].data)
source_table = Table(msa_hdu_list['SOURCE_INFO'].data)

new_slit_id, new_source_id, shutter_table, source_table = add_fs_target(visit_rate_file, 
                                                                        slit_name, ra, 
                                                                        dec, 
                                                                        stellarity=stellarity)

msa_hdu_list['SHUTTER_INFO'] = fits.table_to_hdu(shutter_table)
msa_hdu_list['SOURCE_INFO'] = fits.table_to_hdu(source_table)

msa_hdu_list[2].name = 'SHUTTER_INFO'
msa_hdu_list[3].name = 'SOURCE_INFO'

msa_hdu_list.writeto(msa_file, overwrite=True)
msa_hdu_list.close()
