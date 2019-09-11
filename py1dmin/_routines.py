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
        xyz_file.write(target_geo)

    job_file_path = os.path.join(job_dir_path, 'bath.xyz')
    with open(job_file_path, 'w') as xyz_file:
        xyz_file.write(bath_geo)


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


def write_elstruct_sub(job_dir_path, drive_path, prog):
    """ writes the elstruct submission string
    """

    if prog == 'molpro':
        sub_name = 'm.x'
    else:
        raise NotImplementedError

    # Read the original input script as a string
    elstruct_sub_name_in = os.path.join(drive_path, 'elstruct.x')
    with open(elstruct_sub_name_in, 'w') as elstruct_sub_file_in:
        in_sub_str = elstruct_sub_file_in.read()

    # Write the electronic structure sumbission script
    elstruct_sub_name = os.path.join(job_dir_path, sub_name)
    with open(elstruct_sub_name, 'w') as elstruct_sub_file:
        elstruct_sub_file.write(in_sub_str)

    # Make the auto1dmin.x and elstruct.x file executable
    subprocess.check_call(['chmod', '+x', elstruct_sub_name])


def submit_job(drive_path, run_path, njobs):
    """ submit the 1DMin job
    """

    submit_str = py1dmin.interface.writer.submission_script(
        drive_path, run_path, njobs)

    moldr.util.run_script(submit_str, run_path)


def obtain_overall_lj(run_path):
    """ get the lj params from each run and average them together
    """

    avg_sigma = None
    avg_epsilon = None

    job_dirs = [os.path.join(run_path, directory)
                for directory in os.listdir(run_path)
                if 'build' not in directory and 'yaml' not in directory]
    sigmas, epsilons = [], []
    for job_dir in job_dirs:
        lj_file_name = os.path.join(job_dir, 'lj.out')
        if os.path.exists(lj_file_name):
            with open(lj_file_name, 'r') as lj_file:
                output_string = lj_file.read()
            sigs, epss = py1dmin.interface.reader.lennard_jones(
                output_string)
            if sigs is not None and epss is not None:
                for sig, eps in zip(sigs, epss):
                    sigmas.append(sig)
                    epsilons.append(eps)

    assert len(sigmas) == len(epsilons)
    if sigmas and epsilons:
        avg_sigma = statistics.mean(sigmas)
        avg_epsilon = statistics.mean(epsilons)
        print('{0:<14s}{1:<16s}'.format('\nSigma (Ang)', 'Epsilon (cm-1)'))
        for sig, eps in zip(sigmas, epsilons):
            print('{0:<14.4f}{1:<16.4f}'.format(sig, eps))
        print('\nAverage Sigma =', avg_sigma)
        print('Average Epsilon =', avg_epsilon)
        print('Number of values = ', len(sigmas))
    else:
        print('No Sigma and Epsilon Values obtained =(')

    return avg_sigma, avg_epsilon


def obtain_geometries(run_path):
    """ get the miniumum geometries for each run
    """

    geo_str = ''

    job_dirs = [os.path.join(run_path, directory)
                for directory in os.listdir(run_path)
                if 'build' not in directory and 'yaml' not in directory]
    for job_dir in job_dirs:
        geo_file_name = os.path.join(job_dir, 'lj.out')
        if os.path.exists(geo_file_name):
            with open(geo_file_name, 'r') as geo_file:
                geo_str += geo_file.read()

    return geo_str


def obtain_zero_energy(run_path):
    """ get the zero-energy
    """

    ene = None

    job_dirs = [os.path.join(run_path, directory)
                for directory in os.listdir(run_path)
                if 'build' not in directory and 'yaml' not in directory]
    for job_dir in job_dirs:
        ene_file_name = os.path.join(job_dir, 'zero.ene')
        if os.path.exists(ene_file_name):
            with open(ene_file_name, 'r') as ene_file:
                ene_str = ene_file.read()
            ene = float(ene_str)
            # use autofile above to read the string
            break

    return ene


def read_lj_from_save(target_save_prefix, target_info, theory_level):
    """ lj
    """

    sigma = None
    epsilon = None

    # Set new theory level for the filesystem
    theory_level = [theory_level[1],
                    theory_level[2],
                    moldr.util.orbital_restriction(
                        target_info, theory_level)]

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

    # Set new theory level for the filesystem
    theory_level = [theory_level[1],
                    theory_level[2],
                    moldr.util.orbital_restriction(
                        target_info, theory_level)]

    # Search save file system for LJ params
    tgt_save_fs = autofile.fs.species(target_save_prefix)
    tgt_save_fs.leaf.create(target_info)
    tgt_save_path = tgt_save_fs.leaf.path(target_info)
    etrans_save_fs = autofile.fs.energy_transfer(tgt_save_path)
    etrans_save_fs.leaf.create(theory_level)
    sigma = etrans_save_fs.leaf.file.lennard_jones_sigma.write(
        sigma, theory_level)
    epsilon = etrans_save_fs.leaf.file.lennard_jones_epsilon.write(
        epsilon, theory_level)

    etrans_save_path = etrans_save_fs.leaf.path(theory_level)
    print('\nWriting Lennard-Jones parameters to Save FileSystem at\n')
    print(etrans_save_path)


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
