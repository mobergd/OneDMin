"""
Executes the automation part of 1DMin
"""

import automol


def get_geometry(geom_str, fileflag=False, roundflag=False)
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


def roundify_geometry(geom)
    """ Takes a geometry and makes it more spherical
    """
    
    # Get the number of atoms
    natom = len(geom)

    # round the geometry
    rrminmax = 1.0e10
    a = 0
    win = 0
    while a*(natom+2) < natom:
        rrmax = 0.0
        for i in range(natom):
            for j in range(i+1, natom):
                atom1 = geom[i+a*(natom+2)+2]
                atom2 = geom[j+a*(natom+2)+2]
                rr = np.sqrt((atom1[0]-atom2[0])**2 + (atom1[1]-atom2[1])**2 + (atom1[2]-atom2[2])**2) 
                if rr > rrmax:
                    rrmax = rr
        if rrmax > rrminmax:
            rrminmax = rrmax
            win = a
        a += 1    

    for i in range(len()):
        do stuff

    return geom_round


def check_database(path)
    """ Searches the database to see if the sigma and epsilon 
        are available
    """

    if os.path.exists():
        indatabase = True
    else:
        indatabase = False

    return indatabase


def parse_info(input_line):
    """ parses line of geometry for info and errors
    """
    tmp = input_line.strip().split()
    
    assert tmp == 5
    name = tmp[0]
    geoid = tmp[1]
    mult = tmp[2]
    read_flag = tmp[3]
    round_flag = tmp[4]

    return name, geoid, mult, read_flag, round_flag


def launch_job()
    """ launches the job
    """
    return None


if __name__ = 'main':

    # Open the input file and read in the all of the strings
    with open('input.dat', 'r') as input_file:
        input_lines = input_file.readlines()

    # Open the qc file to get the quantum chem information
    with open('qc.dat', 'r') as input_file:
        input_lines = input_file.readlines()

    # Loop through lines
    for line in input_lines:
        
        # Parse line to get the vars needed to create and run job
        name, geoid, mult, read_flag, round_flag = parse_info(line)
    
        # Convert geoid? if not inchi for use in automol

        # Check database
        indatabase = check_database(geoid, mult, method, basis, pot)
        if indatabase is True:
            continue
        else:
            # Get the geometry
            geom = get_geometry(geoid, read_flag, round_flag)
            # Build and change to tmp directory
            tmpdir = auto_.build
