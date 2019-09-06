"""
reads info
"""

import os
import subprocess
import statistics
import random
import automol
import autofile
import moldr
import elstruct
import py1dmin.interface


def write_input(job_dir_path, nsamp,
                target_name='target.xyz', bath_name='bath,xyz',
                smin=2, smax=6):
    """ write the input file
    """

    ranseed = random.randrange(1E8, 1E9)
    inp_str = py1dmin.interface.writer.onedmin_input(
        ranseed, nsamp, target_name, bath_name, smin, smax)

    job_file_path = os.path.join(job_dir_path, 'input.dat')
    with open(job_file_path, 'w') as input_file:
        input_file.write(inp_str)


def write_xyz(job_dir_path, target_geo, bath_geo):
    """ write the target and bath xyz files
    """

    job_file_path = os.path.join(job_dir_path, 'target.xyz')
    with open(job_file_path, 'w') as xyz_file:
        xyz_file.write(automol.geom.string(target_geo))

    job_file_path = os.path.join(job_dir_path, 'bath.xyz')
    with open(job_file_path, 'w') as xyz_file:
        xyz_file.write(automol.geom.string(bath_geo))


def write_elstruct_inp(job_dir_path,
                       charge, mult, method, basis, thry_lvl,
                       prog, memory):
    """ writes the electronic structure input file
    """
    assert prog in ('g09', 'molpro')

    elstruct_inp_str = elstruct.writer.energy(
        geom='GEOMETRY',
        charge=charge,
        mult=mult,
        method=method,
        basis=basis,
        prog=prog,
        mol_options=('nosym', 'noorient', 'angstrom'),
        memory=memory,
        comment='SAMPLE GEOM',
        orb_restricted=moldr.util.orbital_restriction(
            ['', charge, mult], thry_lvl),
    )

    elstruct_inp_name = os.path.join(job_dir_path, 'qc.mol')
    with open(elstruct_inp_name, 'w') as elstruct_inp_file:
        elstruct_inp_file.write(elstruct_inp_str)


def write_elstruct_sub(job_dir_path, elstruct_sub_str, sub_name):
    """ writes the elstruct submission string
    """

    # Write the electronic structure sumbission script
    elstruct_sub_name = os.path.join(job_dir_path, sub_name)
    with open(elstruct_sub_name, 'w') as elstruct_sub_file:
        elstruct_sub_file.write(elstruct_sub_str)

    # Make the auto1dmin.x and elstruct.x file executable
    subprocess.check_call(['chmod', '+x', elstruct_sub_name])


def submit_job(run_path, nprocs):
    """ submit the 1DMin job
    """

    submit_str = py1dmin.interface.writer.submission_script(nprocs)

    submit_name = os.path.join(run_path, 'onedmin.batch')
    with open(submit_name, 'w') as submit_file:
        submit_file.write(submit_str)

    # moldr.run_script(1dmin_scr)


def obtain_overall_lj(run_path):
    """ get the lj params from each run and average them together
    """

    avg_sigma = None
    avg_epsilon = None

    job_dirs = os.listdir(run_path)
    sigmas, epsilons = [], []
    for job_dir in job_dirs:
        os.chdir(job_dir)
        with open(os.path.join(job_dir, 'lj.dat'), 'r') as lj_file:
            output_string = lj_file.read()
        sig, eps = py1dmin.interface.reader.lennard_jones(output_string)
        if sig is not None and eps is not None:
            sigmas.append(sig)
            epsilons.append(eps)

    avg_sigma = statistics.mean(sigmas)
    avg_epsilon = statistics.mean(epsilons)

    return avg_sigma, avg_epsilon


def read_lj_from_save(target_save_prefix, target_info, theory_level):
    """ lj
    """

    sigma = None
    epsilon = None

    print(theory_level)
    # Search save file system for LJ params
    tgt_save_fs = autofile.fs.species(target_save_prefix)
    if tgt_save_fs.leaf.exists(target_info):
        tgt_save_path = tgt_save_fs.leaf.path(target_info)
        etrans_save_fs = autofile.fs.energy_transfer(tgt_save_path)
        if etrans_save_fs.leaf.exists(theory_level):
            sigma = etrans_save_fs.leaf.file.lennard_jones_sigma.read(
                theory_level)
            epsilon = etrans_save_fs.leaf.file.lennard_jones_epsilon.read(
                theory_level)

    return sigma, epsilon


def write_lj_to_save(sigma, epsilon,
                     target_save_prefix, target_info, theory_level):
    """ lj
    """

    sigma = None
    epsilon = None

    # Search save file system for LJ params
    tgt_save_fs = autofile.fs.species(target_save_prefix)
    if tgt_save_fs.leaf.exists(target_info):
        tgt_save_path = tgt_save_fs.leaf.path(target_info)
        etrans_save_fs = autofile.fs.energy_transfer(tgt_save_path)
        if etrans_save_fs.leaf.exists(theory_level):
            sigma = etrans_save_fs.leaf.file.lennard_jones_sigma.write(
                theory_level)
            epsilon = etrans_save_fs.leaf.file.lennard_jones_epsilon.write(
                theory_level)

    return sigma, epsilon


def get_geometry(spc_info, thry_lvl, save_prefix,
                 geom_dct={}, conf='low', minmax=False):
    """ get the geometry
    """

    assert conf in ('low', 'all')

    # Obtain the reference geometry for the species
    geo = moldr.util.reference_geometry(
        spc_info,
        thry_lvl,
        save_prefix,
        geom_dct=geom_dct)

    # Obtain the desired conformer(s)
    cnf_save_fs = autofile.fs.conformer(save_prefix)
    if conf == 'low':
        min_cnf_locs = moldr.util.min_energy_conformer_locators(save_prefix)
        if min_cnf_locs:
            geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
    elif conf == 'all':
        # add a reader to the trajectory function
        cnf_save_fs.trunk.file.trajectory.read()

    # Obtain the most spherical geometry for species if desired
    if minmax:
        geo = py1dmin.interface.util.roundify_geometry(geo)

    # Format the geoms into xyz strings
    geo_str = automol.geom.string(geo)

    return geo_str
