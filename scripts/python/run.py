"""
Executes the automation part of 1DMin
"""

import os
import subprocess
import random
import elstruct
import moldr
import autofile
import onedmin


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

TGT_SAVE_PREFIX = 'save'
BATH_SAVE_PREFIX = 'save'

# Search save file system for LJ params
tgt_save_fs = autofile.fs.species(TGT_SAVE_PREFIX)
if tgt_save_fs.leaf.exists(TGT_INFO):
    tgt_save_path = tgt_save_fs.leaf.path(TGT_INFO)

etrans_save_fs = autofile.fs.energy_transfer(tgt_save_path)
etrans_save_path = etrans_save_fs.leaf.path(THRY_LVL)
if etrans_save_fs.leaf.exists(THRY_LVL):
    epsilon = etrans_save_fs.leaf.file.lennard_jones_epsilon.read(THRY_LVL)
    sigma = etrans_save_fs.leaf.file.lennard_jones_sigma.read(THRY_LVL)

SIGMA = None
EPSILON = None

# Build a geometry dictionary to potentially get starting geometries
GEOM_DCT = moldr.util.geometry_dictionary(GEOM_DIR)
print('geom dict')
for key, val in GEOM_DCT.items():
    print(val)

# Run 1DMin to calculate the LJ parameters
if SIGMA is None and EPSILON is None or NSAMP >= MIN_NSAMP:

    # Obtain the geometry for the target and bath
    TGT_GEO = onedmin.get_geometry(
        TGT_INFO,
        THRY_LVL,
        TGT_SAVE_PREFIX,
        geom_dct=GEOM_DCT,
        conf='low',
        minmax=False)
    BATH_GEO = onedmin.get_geometry(
        BATH_INFO,
        THRY_LVL,
        BATH_SAVE_PREFIX,
        geom_dct=GEOM_DCT,
        conf='low',
        minmax=False)

    # Write the params to the run file system
    # tgt_save_fs = autofile.fs.species(TGT_RUN_PREFIX)
    # tgt_save_fs.leaf.exists(TGT_INFO)
    # tgt_save_path = tgt_save_fs.leaf.path(TGT_INFO)
    # etrans_run_fs = autofile.fs.energy_transfer(tgt_run_path)
    # etrans_run_path = etrans_run_fs.leaf.path(THRY_LVL)
    # etrans_run_fs.leaf.create(THRY_LVL)
    # etrans_fs.leaf.file.energy.read(locs)
    # etrans_fs.leaf.file.lennard_jones_epsilon.read(locs)
    # etrans_fs.leaf.file.lennard_jones_sigma.read(locs)

    # DO: Check the number of samples already run

    # For each processor, run an instance of 1Dmin
    for i in range(NPROCS):

        # Build run directory
        job_dir_path = os.path.join('run', 'run'+str(i+1))
        os.mkdir(job_dir_path)

        # Build the 1DMin input file string with a random 9-digit integer
        ranseed = random.randrange(1E8, 1E9)
        inp_str = onedmin.write_onedmin_inp(
            ranseed, NSAMP, BATH_GEO, SMIN, SMAX)

        # Write the 1DMin input file
        job_file_path = os.path.join(job_dir_path, 'input.dat')
        with open(job_file_path, 'w') as input_file:
            input_file.write(inp_str)

        # Build the electronic structure template file string
        elstruct_inp_str = elstruct.writer.energy(
            geom='GEOMETRY',
            charge=RUN_CHG,
            mult=RUN_MLT,
            method=RUN_METHOD,
            basis=RUN_BASIS,
            prog=RUN_PROG,
            mol_options=('nosym', 'noorient', 'angstrom'),
            memory=4,
            comment='SAMPLE GEOM',
            machine_options=(),
            orb_restricted=moldr.util.orbital_restriction(
                ['', RUN_CHG, RUN_MLT], THRY_LVL),
            scf_options=(),
            casscf_options=(),
            corr_options=(),
            gen_lines=()
        )

        # Write the electronic structure template file
        elstruct_inp_name = os.path.join(job_dir_path, 'qc.mol')
        with open(elstruct_inp_name, 'w') as elstruct_inp_file:
            elstruct_inp_file.write(elstruct_inp_str)

        # Write the electronic structure sumbission script
        elstruct_sub_name = os.path.join(job_dir_path, SUB_NAME)
        with open(elstruct_sub_name, 'w') as elstruct_sub_file:
            elstruct_sub_file.write(ELSTRUCT_SUB_STR)

        # Copy the auto1dmin.x executable
        auto1dmin_exe_name = os.path.join(job_dir_path, 'auto1dmin.x')
        onedmin.copy_exe(job_dir_path)

        # Make the auto1dmin.x and elstruct.x file executable
        subprocess.check_call(['chmod', '+x', elstruct_sub_name])
        subprocess.check_call(['chmod', '+x', auto1dmin_exe_name])

    # Write the 1DMin batch submission script and run
    SUBMIT_STR = onedmin.write_submission_script(NPROCS)

    # Write the file
    SUBMIT_NAME = os.path.join('run', 'onedmin.batch')
    with open(SUBMIT_NAME, 'w') as submit_file:
        submit_file.write(SUBMIT_STR)
    # moldr.run_script(1dmin_scr)

    # Write the params to the save file system
    # etrans_save_fs = autofile.fs.energy_transfer(tgt_save_path)
    # etrans_save_path = etrans_save_fs.leaf.path(THRY_LVL)
    # etrans_save_fs.leaf.create(THRY_LVL)
    # etrans_fs.leaf.file.energy.read(locs)
    # etrans_fs.leaf.file.lennard_jones_epsilon.read(locs)
    # etrans_fs.leaf.file.lennard_jones_sigma.read(locs)

