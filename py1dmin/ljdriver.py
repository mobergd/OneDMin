"""
Executes the automation part of 1DMin
"""

import os
import sys
import moldr
import autofile
from py1dmin import _routines as odm_routines
from py1dmin import _in_parser as odm_inparser


# Set paths to various working directories
DRIVE_PATH = os.getcwd()
GEOM_PATH = os.path.join(DRIVE_PATH, 'geoms')

# Read the input file into a string
with open(os.path.join(DRIVE_PATH, 'input.dat'), 'r') as infile:
    INPUT_STRING = infile.read()

# Check to see if the parameters are defined
odm_inparser.check_defined_sections(INPUT_STRING)
odm_inparser.check_defined_lennard_jones_keywords(INPUT_STRING)

# Read the targets and baths from the input
TARGET_DCTS = odm_inparser.read_targets(INPUT_STRING)
BATH_LST = odm_inparser.read_baths(INPUT_STRING)

# Read the run parameters
THEORY_LEVEL = odm_inparser.read_theory_level(INPUT_STRING)
POTENTIAL = odm_inparser.read_potential(INPUT_STRING)
NJOBS = odm_inparser.read_njobs(INPUT_STRING)
NSAMPS = odm_inparser.read_nsamps(INPUT_STRING)
SMIN = odm_inparser.read_smin(INPUT_STRING)
SMAX = odm_inparser.read_smax(INPUT_STRING)
CONFS = odm_inparser.read_confs(INPUT_STRING)
TARGET_SAVE_PREFIX = odm_inparser.read_save_prefix(INPUT_STRING)
BATH_SAVE_PREFIX = odm_inparser.read_save_prefix(INPUT_STRING)
RUN_PREFIX = odm_inparser.read_run_prefix(INPUT_STRING)

# Read the theory.dat file into a string
with open(os.path.join(DRIVE_PATH, 'theory.dat'), 'r') as theoryfile:
    THEORY_STRING = theoryfile.read()

# Check theory level parameters
odm_inparser.check_defined_theory_level_keywords(THEORY_STRING, THEORY_LEVEL)

# Read the theory level parameters
RUN_PROG = odm_inparser.read_program(THEORY_STRING, THEORY_LEVEL)
RUN_METHOD = odm_inparser.read_method(THEORY_STRING, THEORY_LEVEL)
RUN_BASIS = odm_inparser.read_basis(THEORY_STRING, THEORY_LEVEL)
RUN_ORB_REST = odm_inparser.read_orb_restrict(THEORY_STRING, THEORY_LEVEL)
RUN_MEMORY = odm_inparser.read_memory(THEORY_STRING, THEORY_LEVEL)
RUN_NPROCS = odm_inparser.read_nprocs(THEORY_STRING, THEORY_LEVEL)

# Set additional pieces of info using input keywords
THEORY_INFO = [RUN_PROG, RUN_METHOD, RUN_BASIS, RUN_ORB_REST]


# Loop over the species and launch all of the desired jobs
for target, target_info in TARGET_DCTS.items():

    # Set the new information
    RUN_CHG = target_info[1] + BATH_LST[1]
    RUN_MLT = max([target_info[2], BATH_LST[2]])
    BATH_INFO = BATH_LST

    # Search save file system for LJ params
    SIGMA, EPSILON = odm_routines.read_lj_from_save(
        TARGET_SAVE_PREFIX, target_info, THEORY_INFO)

    # Run 1DMin to calculate the LJ parameters
    if SIGMA is None and EPSILON is None:  # or NSAMP <= nsamp_run:

        # Obtain the geometry for the target and bath
        TARGET_GEO = odm_routines.get_geometry(
            target_info,
            THEORY_INFO,
            TARGET_SAVE_PREFIX,
            geom_dct=moldr.util.geometry_dictionary(GEOM_PATH),
            conf=CONFS)
        BATH_GEO = odm_routines.get_geometry(
            BATH_INFO,
            THEORY_INFO,
            BATH_SAVE_PREFIX,
            geom_dct=moldr.util.geometry_dictionary(GEOM_PATH),
            conf=CONFS)

        # Write the params to the run file system
        FS_THEORY_INFO = [THEORY_INFO[1],
                          THEORY_INFO[2],
                          moldr.util.orbital_restriction(
                              target_info, THEORY_INFO)]
        tgt_run_fs = autofile.fs.species(RUN_PREFIX)
        tgt_run_fs.leaf.create(target_info)
        tgt_run_path = tgt_run_fs.leaf.path(target_info)
        etrans_run_fs = autofile.fs.energy_transfer(tgt_run_path)
        etrans_run_path = etrans_run_fs.leaf.path(FS_THEORY_INFO)
        etrans_run_fs.leaf.create(FS_THEORY_INFO)

        # Run an instancw of 1DMin for each processor
        for i in range(NJOBS):

            # Build run directory
            job_dir_path = os.path.join(
                etrans_run_path, 'run{0}'.format(str(i+1)))
            os.mkdir(job_dir_path)
            print('\n\nWriting files to'+job_dir_path)

            # Write the 1DMin input file
            print('  Writing input files...')
            odm_routines.write_input(
                job_dir_path, NSAMPS,
                target_name='target.xyz', bath_name='bath.xyz',
                smin=SMIN, smax=SMAX)

            # Write the geometry files
            print('  Writing xyz files for target and bath....')
            odm_routines.write_xyz(
                job_dir_path, TARGET_GEO, BATH_GEO)

            # Write the electronic structure template file
            print('  Writing electronic structure submission inp template...')
            odm_routines.write_elstruct_inp(
                job_dir_path,
                RUN_CHG, RUN_MLT, RUN_METHOD, RUN_BASIS, THEORY_INFO,
                RUN_PROG, RUN_MEMORY)

            # Write the electronic structure sumbission script
            print('  Writing electronic structure submission script...')
            odm_routines.write_elstruct_sub(
                job_dir_path, DRIVE_PATH, RUN_PROG)

        # Submit the job
        print('\n\nRunning each OneDMin job...')
        odm_routines.submit_job(DRIVE_PATH, etrans_run_path, NJOBS)

        # Read the lj from all the outputs
        print('\n\nAll OneDMin jobs finished.')
        print('\nReading the Lennard-Jones parameters...')
        SIGMA, EPSILONS = odm_routines.obtain_overall_lj(etrans_run_path)
        if SIGMA is None and EPSILON is None:
            print('\n\nExiting OneDMin...')
            sys.exit()

        # Write the params to the save file system
        odm_routines.write_lj_to_save(
            SIGMA, EPSILON,
            TARGET_SAVE_PREFIX, target_info, THEORY_INFO)
        print('\n\nExiting OneDMin...')
