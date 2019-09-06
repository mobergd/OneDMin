"""
Executes the automation part of 1DMin
"""

import os
import moldr
import autofile
from py1dmin import _routines as odm_routines


# Set paths to various working directories
RUN_DIR = os.getcwd()
GEOM_DIR = os.path.join(RUN_DIR, 'geoms')
SRC_PATH = '../../auto1dmin.x'

# Set species data
TGT_ICH = 'InChI=1S/CH4/h1H4'
TGT_CHG = 0
TGT_MLT = 1
TGT_INFO = [TGT_ICH, TGT_CHG, TGT_MLT]
ROUND_GEOM = False

# Set Bath data
BATH_ICH = 'InChI=1S/N2/c1-2'
BATH_CHG = 0
BATH_MLT = 1
BATH_INFO = [BATH_ICH, BATH_CHG, BATH_MLT]
RND_BATH = False

# Set Run Params
NPROCS = 10
MEMORY = 4
MIN_NSAMP = 0
NSAMP = 2
POT = 'lj126'
RUN_PROG = 'molpro'
RUN_METHOD = 'mp2'
RUN_BASIS = 'aug-cc-pvdz'
THRY_LVL = [RUN_PROG, RUN_METHOD, RUN_BASIS, 'RR']
RUN_MLT = max([TGT_MLT, BATH_MLT])
RUN_CHG = TGT_CHG + BATH_CHG
RND_RUN = False
SUB_NAME = 'm.x'
SMIN = 2
SMAX = 6
ELSTRUCT_SUB_STR = ("molpro -n 1 --nouse-logfile --no-xml-output"
                    "-L /soft/molpro/2015.1_170920/bebop/"
                    "molprop_2015_1_linux_x86_64_i8/lib/"
                    "-d /scratch/$USER -o qc.out -s qc.in")
TGT_XYZ_NAME = 'target.xyz'
BATH_XYZ_NAME = 'bath.xyz'
RUN_PATH = os.getcwd()
CONFS = 'low'

RUN_PREFIX = 'run'
TGT_SAVE_PREFIX = 'save'
BATH_SAVE_PREFIX = 'save'


# Search save file system for LJ params
SIGMA, EPSILON = odm_routines.read_lj_from_save(
        TGT_SAVE_PREFIX, TGT_INFO, THRY_LVL)

# Run 1DMin to calculate the LJ parameters
if SIGMA is None and EPSILON is None:  # or NSAMP <= nsamp_run:

    # Obtain the geometry for the target and bath
    TARGET_GEO = odm_routines.get_geometry(
        TGT_INFO,
        THRY_LVL,
        TGT_SAVE_PREFIX,
        geom_dct=moldr.util.geometry_dictionary(GEOM_DIR),
        conf=CONFS)
    BATH_GEO = odm_routines.get_geometry(
        BATH_INFO,
        THRY_LVL,
        BATH_SAVE_PREFIX,
        geom_dct=moldr.util.geometry_dictionary(GEOM_DIR),
        conf=CONFS)

    # Write the params to the run file system
    tgt_run_fs = autofile.fs.species(RUN_PREFIX)
    tgt_run_fs.leaf.create(TGT_INFO)
    tgt_run_path = tgt_run_fs.leaf.path(TGT_INFO)
    etrans_run_fs = autofile.fs.energy_transfer(tgt_run_path)
    etrans_run_path = etrans_run_fs.leaf.path(THRY_LVL)
    etrans_run_fs.leaf.create(THRY_LVL)

    # Run an instancw of 1DMin for each processor
    for i in range(NPROCS):

        # Build run directory
        job_dir_path = os.path.join(RUN_PATH, 'run{0}'.format(str(i+1)))
        os.mkdir(job_dir_path)

        # Write the 1DMin input file
        odm_routines.write_input(
            job_dir_path, NSAMP,
            target_name=TGT_XYZ_NAME, bath_name=BATH_XYZ_NAME,
            smin=SMIN, smax=SMAX)

        # Write the geometry files
        odm_routines.write_xyz(
            job_dir_path, TARGET_GEO, BATH_GEO)

        # Write the electronic structure template file
        odm_routines.write_elstruct_inp(
            job_dir_path,
            RUN_CHG, RUN_MLT, RUN_METHOD, RUN_METHOD, THRY_LVL,
            RUN_PROG, MEMORY)

        # Write the electronic structure sumbission script
        odm_routines.write_elstruct_sub(
            job_dir_path, ELSTRUCT_SUB_STR, SUB_NAME)

    # Submit the job
    odm_routines.submit_job(RUN_PATH, NPROCS)

    # Read the lj from all the outputs
    SIGMA, EPSILON = odm_routines.obtain_overall_lj(RUN_PATH)

    # Write the params to the save file system
    SIGMA, EPSILON = odm_routines.read_lj_from_save(
        TGT_SAVE_PREFIX, TGT_INFO, THRY_LVL)
