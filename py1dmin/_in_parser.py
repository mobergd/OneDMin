"""
parses the input file for keywords
"""

from autoparse.find import first_capture
# from autoparse.find import all_captures
from autoparse.pattern import capturing
from autoparse.pattern import zero_or_more
from autoparse.pattern import one_or_more
# from autoparse.pattern import NONNEWLINE
from autoparse.pattern import NONSPACE
from autoparse.pattern import SPACE
from autoparse.pattern import WILDCARD
from autoparse.pattern import INTEGER
from autoparse.pattern import LINE_FILL
from autoparse.pattern import NEWLINE
# from autoparse.pattern import FLOAT


def read_targets(input_string):
    """ builds a dictionary containing all needed info for the targets
    """
    targets_section = _get_targets_section(input_string)

    targets_dct = {}
    for line in targets_section.splitlines():
        tmp = line.strip().split()
        assert len(tmp) >= 4
        name, ich, chg, mult = tmp[0], tmp[1], tmp[2], tmp[3]
        targets_dct[name] = [ich, chg, mult]

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
        baths_lst = [ich, chg, mult]

    assert baths_lst

    return baths_lst


def _get_targets_section(input_string):
    """ grabs the section of text containing all of the targets
    """
    pattern = ('targets' + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               'end')
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


def _get_baths_section(input_string):
    """ grabs the section of text containing all of the baths
    """
    pattern = ('baths' + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               'end')
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


def read_level(input_string):
    """ obtain the theory level
    """
    pattern = ('theory_level'
               + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_potential(input_string):
    """ obtain the potential to be used
    """

    pattern = ('potential'
               + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword == 'lj126'

    return keyword


def read_nsamps(input_string):
    """ obtain the nsamps to be used
    """

    pattern = ('nsamps'
               + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_njobs(input_string):
    """ obtain the njobs to be used
    """

    pattern = ('njobs'
               + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_smin(input_string):
    """ obtain the smin to be used
    """

    pattern = ('smin'
               + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
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

    pattern = ('smax'
               + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    if keyword is None:
        keyword = 2
    else:
        keyword = int(keyword)

    return keyword


def read_confs(input_string):
    """ obtain the confs to be used
    """

    pattern = ('confs'
               + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_run_prefix(input_string):
    """ obtain the run_prefix to be used
    """

    pattern = ('run_prefix'
               + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_save_prefix(input_string):
    """ obtain the save_prefix to be used
    """

    pattern = ('save_prefix'
               + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_lennard_jones_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def _get_lennard_jones_options_section(input_string):
    """ grabs the section of text containing all of the job keywords
        for lennard jones calculations
    """
    pattern = ('lennard_jones' + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               'end')
    section = first_capture(pattern, input_string)

    assert section is not None

    return section
