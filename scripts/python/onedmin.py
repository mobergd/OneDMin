"""
Executes the automation part of 1DMin
"""

import os
import shutil
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')

def get_geometry(SPC_INFO,
                 THRY_LVL,
                 SPC_SAVE_PREFIX,
                 geom_dct=GEOM_DCT):
    """ get the geometry
    """

    # function incorrect; wont grab all the confs different ones for that

    # Obtain the reference geometry for the species
    spc_geo = moldr.util.reference_geometry(
        SPC_INFO,
        THRY_LVL,
        SPC_SAVE_PREFIX,
        geom_dct=GEOM_DCT)
    print('spc geo; ref geom call')
    print(spc_geo)

    # Obtain the most spherical geometry for species if desired
    if RND_RUN:
        spc_geo = onedmin.roundify_geom(spc_geo)
    
    # Format the geoms into xyz strings
    spc_geo = automol.geom.string(spc_geo)

    return geom

def roundify_geometry(output_string):
    """ Finds the smallest geometry (by volume) from a list of
        conformer geometries
        :param str output_string: string for a conformer
                                  trajectory file
        :r
    """

    # Get lines for the output string
    lines = output_string.splitlines()

    # Get the number of atoms
    natom = int(lines[0])

    # loop over the lines to find the smallest geometry
    rrminmax = 1.0e10
    ngeom = 0
    small_geo_idx = 0
    while ngeom*(natom+2) < len(lines):
        rrmax = 0.0
        for i in range(natom):
            for j in range(i+1, natom):
                # Get the line
                xyz1 = lines[i+ngeom*(natom+2)+2].split()[1:]
                xyz2 = lines[j+ngeom*(natom+2)+2].split()[1:]
                # Get the coordinates
                atom1 = [float(val) for val in xyz1]
                atom2 = [float(val) for val in xyz2]
                # Calculate the interatomic distance
                rrtest = np.sqrt((atom1[0]-atom2[0])**2 +
                                 (atom1[1]-atom2[1])**2 +
                                 (atom1[2]-atom2[2])**2)
                # Check and see if distance is more than max
                if rrtest > rrmax:
                    rrmax = rrtest
        # If max below moving threshold, set to smallest geom
        if rrmax < rrminmax:
            rrminmax = rrmax
            small_geo_idx = ngeom
        ngeom += 1

    # Set the output geometry
    geom_str = '{}\n'.format(natom)
    for i in range(natom):
        geom_str += lines[i+small_geo_idx*(natom+2)+2] + '\n'

    return geom_str


def write_submission_script(nprocs):
    """ launches the job
    """

    # Set the dictionary for the 1DMin input file
    fill_vals = {
        "nprocs": nprocs,
        "scratch": '/scratch/$USER'
    }

    # Set template name and path for the 1dmin input file
    template_file_name = '1dmin_submit.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build the 1dmin input string
    sub_str = Template(filename=template_file_path).render(**fill_vals)

    return sub_str


def write_onedmin_inp(seed, nsamp, bath, smin, smax):
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

    # Set template name and path for the 1dmin input file
    template_file_name = '1dmin_inp.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build the 1dmin input string
    input_str = Template(filename=template_file_path).render(**fill_vals)

    return input_str


def copy_exe(path):
    """ copy the auto1dmin.x file
    """

    exe_src = os.path.join(SRC_PATH, 'auto1dmin.x')
    exe_dest = os.path.join(path, 'auto1dmin.x')
    shutil.copyfile(exe_src, exe_dest)
