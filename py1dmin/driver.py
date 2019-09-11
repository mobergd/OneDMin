"""
Executes the automation part of 1DMin
"""

import os
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

# Read the targets and baths from the input
TARGET_DCTS = odm_inparser.read_targets(INPUT_STRING)
BATH_LST = odm_inparser.read_baths(INPUT_STRING)

# Read the run parameters
THEORY = odm_inparser.read_level(INPUT_STRING)
POTENTIAL = odm_inparser.read_potential(INPUT_STRING)
NJOBS = odm_inparser.read_njobs(INPUT_STRING)
NSAMPS = odm_inparser.read_nsamps(INPUT_STRING)
SMIN = odm_inparser.read_smin(INPUT_STRING)
SMAX = odm_inparser.read_smax(INPUT_STRING)
CONFS = odm_inparser.read_confs(INPUT_STRING)
TARGET_SAVE_PREFIX = odm_inparser.read_save_prefix(INPUT_STRING)
BATH_SAVE_PREFIX = odm_inparser.read_save_prefix(INPUT_STRING)
RUN_PREFIX = odm_inparser.read_run_prefix(INPUT_STRING)

# Read the the theory level parameters

# Set additional pieces of info using input keywords
# THEORY_LEVEL = [RUN_PROG, RUN_METHOD, RUN_BASIS, 'RR']
#
#
# FS_THEORY_LEVEL = [THEORY_LEVEL[1],
#                THEORY_LEVEL[2],
#                moldr.util.orbital_restriction(
#                    TARGET_INFO, THEORY_LEVEL)]
THEORY_LEVEL = ''
FS_THEORY_LEVEL = ''
BATH_INFO = BATH_LST

# Loop over the species and launch all of the desired jobs
for target in TARGET_DCTS:

    # Set the combined multiplicity and charge for the run
    RUN_CHG = target[0] + BATH_LST[1]
    RUN_MLT = max([target[1], BATH_LST[1]])
    TARGET_INFO = target

    # Search save file system for LJ params
    SIGMA, EPSILON = odm_routines.read_lj_from_save(
        TARGET_SAVE_PREFIX, TARGET_INFO, THEORY_LEVEL)

    # Run 1DMin to calculate the LJ parameters
    if SIGMA is None and EPSILON is None:  # or NSAMP <= nsamp_run:

        # Obtain the geometry for the target and bath
        TARGET_GEO = odm_routines.get_geometry(
            TARGET_INFO,
            THEORY_LEVEL,
            TARGET_SAVE_PREFIX,
            geom_dct=moldr.util.geometry_dictionary(GEOM_PATH),
            conf=CONFS)
        BATH_GEO = odm_routines.get_geometry(
            BATH_INFO,
            THEORY_LEVEL,
            BATH_SAVE_PREFIX,
            geom_dct=moldr.util.geometry_dictionary(GEOM_PATH),
            conf=CONFS)

        # Write the params to the run file system
        tgt_run_fs = autofile.fs.species(RUN_PREFIX)
        tgt_run_fs.leaf.create(TARGET_INFO)
        tgt_run_path = tgt_run_fs.leaf.path(TARGET_INFO)
        etrans_run_fs = autofile.fs.energy_transfer(tgt_run_path)
        etrans_run_path = etrans_run_fs.leaf.path(FS_THEORY_LEVEL)
        etrans_run_fs.leaf.create(FS_THEORY_LEVEL)

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
                job_dir_path, NSAMP,
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
                RUN_CHG, RUN_MLT, RUN_METHOD, RUN_BASIS, THEORY_LEVEL,
                RUN_PROG, MEMORY)

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
            TARGET_SAVE_PREFIX, TARGET_INFO, THEORY_LEVEL)
        print('\n\nExiting OneDMin...')
