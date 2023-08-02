import os
import sys
import argparse
import logging
import json
from glob import glob


def read_kpp_config(kpp_dir):
    """
    Read all KPP config files in a directory

    Parameters
        (str) kpp_dir: KPP directory

    Returns
        (list of str): all lines from all config files
    """

    suffixes = ['.kpp', '.spc', '.eqn', '.def']

    lines = list()

    for suffix in suffixes:
        files = glob(os.path.join(kpp_dir, '*' + suffix))
        logging.debug(files)
        for filename in files:
            f = open(filename, 'r')
            lines.extend(f.readlines())

    # remove empty lines and tabs
    lines = [line.replace('\t', '') for line in lines if line.strip()] 

    for line in lines:
        logging.debug(line.strip())

    return lines


def split_by_section(lines):
    """
    Split KPP config lines by section

    Parameters
        (list of str) lines: all lines config files

    Returns
        (dict of list of str): lines in each section
    """

    sections = {'#ATOMS': [],
                '#DEFVAR': [],
                '#DEFFIX': [],
                '#EQUATIONS': []}

    joined_lines = ''.join(lines)
    section_blocks = joined_lines.split('#')

    for section in sections:
        for section_block in section_blocks:
            if section.replace('#', '') in section_block:
                sections[section].extend(section_block.split('\n')[1:-1])

    return sections


def read_species_template(template_file):

    with open(template_file, 'r') as f:
        species_template = json.load(f)
    logging.debug(species_template)


if __name__ == '__main__':

    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--logfile', type=str,
        default=sys.stdout,
        help='log file (default stdout)')
    parser.add_argument('--kpp_dir', type=str,
        default=os.path.join('..', 'configs', 'kpp'),
        help='KPP config directory')
    parser.add_argument('--species_template', type=str,
        default='species_template.json',
        help='MICM JSON species template file')
    parser.add_argument('--debug', action='store_true',
        help='set logging level to debug')
    args = parser.parse_args()

    """
    Setup logging
    """
    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(stream=args.logfile, level=logging_level)

    """
    Read KPP config files
    """
    lines = read_kpp_config(args.kpp_dir)


    """
    Split KPP config by section
    """
    sections = split_by_section(lines)
    for section in sections:
        logging.info('____ KPP section %s ____' % section)
        for line in sections[section]:
            logging.info(line)
        print('\n')

    """
    Read species template file
    """
    read_species_template(args.species_template) 

