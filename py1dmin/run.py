"""
Executes the automation part of 1DMin
"""

import tempfile
import autodir


# Set variables needed for run
name = 'H2O'
geoid = 'O'
mult = 1
nsamp = 40
bath = 'N2'
pot = 'lj126'
iread = False
isphere = False

# NJOBs = sys.argv[1]
# while jobs_submit < NJOBS: (do something below)

# Check if the species is in the database
DBANK_DIR = '/path/to/dbnk'
dbank_dir_path = autodir.lj.directory_path(DBANK_DIR, bath, pot)
if os.path.exists(dbank_dir_path):
    eps = autodir.lj.read_epsilon_file(DBANK_DIR, bath, pot)
    sig = autodir.lj.read_sigma_file(DBANK_DIR, bath, pot)
    # continue loop
    # might just give message saying they are found for run

# Set temp directory name
TMP_DIR = tempfile.mkdtemp()

# Set the temp directory path
dir_path = autodir.lj.directory_path(TMP_DIR, bath, pot)

# Create the temp directory
autodir.lj.create(TMP_DIR, bath, pot)

# Get the xyz coordinates for the target geometry
geom_target = get_geometry(geoid, iread, isphere)
geom_bath = get_bath(bath)

# Run the LJ-126 with direct dynamics
run_pot(geom_target, geom_bath, mult, nsamp, pot='lj126')
# assert pot == 'lj126'

 # Open the input file and read in the all of the strings
 with open('input.dat', 'r') as input_file:
     INPUT_LINES = input_file.readlines()

 # Open the qc file to get the quantum chem information
 with open('qc.dat', 'r') as qc_file:
     QC_LINES = qc_file.readlines()

 # Loop through lines
 for line in INPUT_LINES:

