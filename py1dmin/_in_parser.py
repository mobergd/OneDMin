"""
parses the input file for keywords
"""

from autoparse.find import first_capture
from autoparse.find import all_captures
from autoparse.pattern import capturing
from autoparse.pattern import zero_or_more
from autoparse.pattern import one_or_more
from autoparse.pattern import escape
from autoparse.pattern import NONSPACE
from autoparse.pattern import SPACE
from autoparse.pattern import WILDCARD
from autoparse.pattern import INTEGER
from autoparse.pattern import LINE_FILL
from autoparse.pattern import NONNEWLINE
from autoparse.pattern import NEWLINE


ODM_SUPPORTED_SECTIONS = [
    'lennard_jones',
    'properties',
    'baths',
    'targets'
]
ODM_REQUIRED_SECTIONS = [
    'lennard_jones',
    'baths',
    'targets'
]

ODM_SUPPORTED_KEYWORDS = [
    'theory_level',
    'potential',
    'nsamps',
    'njobs',
    'smin',
    'smax',
    'run_prefix',
    'save_prefix',
    'confs'
]
ODM_REQUIRED_KEYWORDS = [
    'theory_level',
    'potential',
    'nsamps',
    'njobs',
    'run_prefix',
    'save_prefix',
    'confs'
]

ELSTRUCT_REQUIRED_KEYWORDS = [
    'program',
    'method',
    'basis',
    'orb_restrict'
]
ELSTRUCT_SUPPORTED_KEYWORDS = [
    'program',
    'method',
    'basis',
    'orb_restrict',
    'memory',
    'nprocs',
    'econv',
    'gconv'
]


# Parser functions for the OneDMin Input files

# Read the targets and baths

def read_targets(input_string):
    """ builds a dictionary containing all needed info for the targets
    """
    targets_section = _get_targets_section(input_string)

    targets_dct = {}
    for line in targets_section.splitlines():
        tmp = line.strip().split()
        assert len(tmp) >= 4
        name, ich, chg, mult = tmp[0], tmp[1], tmp[2], tmp[3]
        targets_dct[name] = [ich, int(chg), int(mult)]

    assert targets_dct

    return targets_dct


def read_baths(input_string):
    """ builds a dictionary containing all needed info for the baths
    """
    baths_section = _get_baths_section(input_string)

    baths_lst = []
    for line in baths_section.splitlines():
        tmp = line.strip().split()
        assert len(tmp) >= 4
        ich, chg, mult = tmp[1], tmp[2], tmp[3]
        baths_lst = [ich, int(chg), int(mult)]

    assert baths_lst

    return baths_lst


def _get_targets_section(input_string):
    """ grabs the section of text containing all of the targets
    """
    pattern = (escape('$targets') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


def _get_baths_section(input_string):
    """ grabs the section of text containing all of the baths
    """
    pattern = (escape('$baths') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


# Read the keywords from the lennard jones section

def read_potential(input_string):
    """ obtain the potential to be used
    """

    pattern = ('potential' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword == 'lj126'

    return keyword


def read_nsamps(input_string):
    """ obtain the nsamps to be used
    """

    pattern = ('nsamps' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_njobs(input_string):
    """ obtain the njobs to be used
    """

    pattern = ('njobs' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_smin(input_string):
    """ obtain the smin to be used
    """

    pattern = ('smin' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    if keyword is None:
        keyword = 2
    else:
        keyword = int(keyword)

    return keyword


def read_smax(input_string):
    """ obtain the smax to be used
    """

    pattern = ('smax' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    if keyword is None:
        keyword = 6
    else:
        keyword = int(keyword)

    return keyword


def read_confs(input_string):
    """ obtain the confs to be used
    """

    pattern = ('confs' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_run_prefix(input_string):
    """ obtain the run_prefix to be used
    """

    pattern = ('run_prefix' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_save_prefix(input_string):
    """ obtain the save_prefix to be used
    """

    pattern = ('save_prefix' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_theory_level(input_string):
    """ obtain the theory level
    """
    pattern = ('theory_level' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def _get_lennard_jones_options_section(input_string):
    """ grabs the section of text containing all of the job keywords
        for lennard jones calculations
    """
    pattern = (escape('$lennard_jones') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


# Parser functions for the theory.dat file

def read_program(input_string, level):
    """ obtain the theory level
    """
    pattern = ('program' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_method(input_string, level):
    """ obtain the method level
    """
    pattern = ('method' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_basis(input_string, level):
    """ obtain the basis level
    """
    pattern = ('basis' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_orb_restrict(input_string, level):
    """ obtain the orb_restricted level
    """
    pattern = ('orb_restrict' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_memory(input_string, level):
    """ obtain the memory level
    """
    pattern = ('memory' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    if keyword is None:
        keyword = 0.5
    else:
        keyword = float(keyword)

    return keyword


def read_nprocs(input_string, level):
    """ obtain the memory level
    """
    pattern = ('memory' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_level_section(input_string, level)

    keyword = first_capture(pattern, block)

    if keyword is None:
        keyword = 1
    else:
        keyword = int(keyword)

    assert keyword is not None

    return keyword


def _get_level_section(input_string, level):
    """ species input
    """
    pattern = ('level' + one_or_more(SPACE) + level +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               'end')
    section = first_capture(pattern, input_string)

    return section


# Functions to check for errors in the input file

def check_defined_sections(input_string):
    """ verify all defined sections have been defined
    """
    pattern = (escape('$') + capturing(one_or_more(NONNEWLINE)))

    matches = all_captures(pattern, input_string)

    # See if each section has an paired end and is a supported keywords
    defined_sections = []
    for i, match in enumerate(matches):
        if (i+1) % 2 == 0:
            if match != 'end':
                raise ValueError
        else:
            defined_sections.append(match)

    # Check if sections are supported
    if not all(section in ODM_SUPPORTED_SECTIONS
               for section in defined_sections):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(section in defined_sections
               for section in ODM_REQUIRED_SECTIONS):
        raise NotImplementedError


def check_defined_lennard_jones_keywords(input_string):
    """ obtains the keywords defined in the input by the user
    """
    section_string = _get_lennard_jones_options_section(input_string)
    defined_keywords = _get_defined_keywords(section_string)

    # Check if keywords are supported
    if not all(keyword in ODM_SUPPORTED_KEYWORDS
               for keyword in defined_keywords):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(keyword in defined_keywords
               for keyword in ODM_REQUIRED_KEYWORDS):
        raise NotImplementedError


def check_defined_theory_level_keywords(input_string, level):
    """ obtains the keywords defined in the input by the user
    """
    section_string = _get_level_section(input_string, level)
    defined_keywords = _get_defined_keywords(section_string)

    # Check if keywords are supported
    if not all(keyword in ELSTRUCT_SUPPORTED_KEYWORDS
               for keyword in defined_keywords):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(keyword in defined_keywords
               for keyword in ELSTRUCT_REQUIRED_KEYWORDS):
        raise NotImplementedError


def _get_defined_keywords(section_string):
    """ gets a list of all the keywords defined in a section
    """

    defined_keys = []
    for line in section_string.splitlines():
        if '=' in line:
            tmp = line.strip().split('=')[0]
            defined_keys.append(tmp.strip())

    return defined_keys
