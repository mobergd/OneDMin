"""
Executes the automation part of 1DMin
"""

import tempfile
import autodir


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


def roundify_geometry(geom):
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
                rr = np.sqrt((atom1[0]-atom2[0])**2 +
                             (atom1[1]-atom2[1])**2 +
                             (atom1[2]-atom2[2])**2)
                if rr > rrmax:
                    rrmax = rr
        if rrmax > rrminmax:
            rrminmax = rrmax
            win = a
        a += 1

#    for i in range(len()):
#        print(a)

    return geom_round


def check_database(path):
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


def launch_job(ncpus):
    """ launches the job
    """

    # Set the dictionary for the 1DMin input file
    fill_vals = {
        "ncpus": ncpus
    }

    # Fill the template
    tempname = 'submit.mako'
    submit_str = Template(filename=template_file_path).render(**fill_vals)

    # Write the file
    with open('batch_submit.sh', 'w') as input_file:
        input_file.write(submit_str)

    # Submit the job to Slurm
    subprocess.check_call(['sbatch', 'batch_submit.sh'])

    return None


def write_1dmin_inp(seed, nsamp, bath, smin, smax):
    """ writes the 1dmin input file for each instance
    """

    # Set the dictionary for the 1DMin input file
    fill_vals = {
        "ranseed": seed,
        "nsamples": nsamp,
        "bath": bath,
        "smin": smin,
        "smax": smax
    }

    # Fill the template
    tempname = '1dmin_inp.mako'
    input_str = Template(filename=template_file_path).render(**fill_vals)

    # Write the file
    with open('1dmin.inp', 'w') as input_file:
        input_file.write(input_str)

    return None


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

