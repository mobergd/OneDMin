"""
functions for the driver
"""

def get_geometry(geom_str, fileflag=False, roundflag=False):
    """ Gets the geometry information.
    """

    if file_exists is False:

        if geom_str is smiles:
            # Convert to Inchi string; get geometry
            ich = automol.smiles.inchi(geom_str)
            geom = automol.inchi.geometry(ich)
        if geom_str is inchi:
            # Use Inchi string to get geometry
            geom = automol.inchi.geometry(geom_str)

    else:
        geom = read_geom_path

    if roundflag is True:
        geom = rot_geom(geom)

    return geom


def parse_input_file(filename):
    """ Opens the Auto1DMIN options file and reads in the keywords.
    """

    # Read the file into a string
    with open(option_filename, 'r') as optfile:
        opt_str = optfile.read()

    # Get the block of qc options
    qc_option_begin = '$QC'
    qc_option_end = '$END'
    qc_option_block(qc_option_begin,
                    qc_option_end,
                    opt_str)

    # Get the block of qc options
    prog_option_begin = '$1DMIN'
    prog_option_end = '$END'
    prog_option_block(prog_option_begin,
                      prog_option_end,
                      opt_str)

    # Search the blocks for keywords

    return stuff

