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
            with open(filename, 'r') as f:
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


def micm_species_json(lines, fixed=False, tolerance=1.0e-12):
    """
    Generate MICM species JSON

    Parameters
        (list of str) lines: lines of species section

    Returns
        (list of dict): list of MICM species entries
    """

    species_json = list() # list of dict

    for line in lines:
        lhs, rhs = tuple(line.split('='))
        logging.debug((lhs, rhs))
        species_dict = {'name': lhs.strip().lstrip(), 'type': 'CHEM_SPEC'}
        if fixed:
            species_dict['tracer type'] = 'CONSTANT'
        else:
            species_dict['absolute tolerance'] = tolerance
        species_json.append(species_dict)

    return species_json


def micm_equation_json(lines):
    """
    Generate MICM equation JSON

    Parameters
        (list of str) lines: lines of equation section

    Returns
        (list of dict): list of MICM equation entries
    """

    equations = list() # list of dict

    for line in lines:
        lhs, rhs = tuple(line.split('='))
        reactants = lhs.split('+')
        products = rhs.split('+')

        # extract equation label delimited by < >
        label, reactants[0] = tuple(reactants[0].split('>'))

        # extract equation coefficients delimited by :
        products[-1], coeffs = tuple(products[-1].split(':'))
        coeffs = coeffs.replace('(', '').replace(')', '').replace(';', '')

        # remove trailing and leading whitespace
        reactants = [reactant.strip().lstrip() for reactant in reactants]
        products = [product.strip().lstrip() for product in products]

        equation_dict = dict()

        if 'SUN' in coeffs:
            equation_dict['type'] = 'PHOTOLYSIS' 
        else:
            equation_dict['type'] = 'ARRHENIUS' 
            equation_dict['A'] = float(coeffs)
            equation_dict['B'] = 0.0

        equation_dict['reactants'] = dict()
        equation_dict['products'] = dict()

        for reactant in reactants:
            if reactant[0].isdigit():
                equation_dict['reactants'][reactant[1:]] \
                    = {'yield': float(reactant[0])}
            elif 'hv' in reactant:
                pass
            else:
                equation_dict['reactants'][reactant] = dict()

        for product in products:
            if product[0].isdigit():
                equation_dict['products'][product[1:]] \
                    = {'yield': float(product[0])}
            else:
                equation_dict['products'][product] = dict()

        equations.append(equation_dict)

    return equations


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
    Generate MICM species JSON from KPP #DEFFIX section
    """
    deffix_json = micm_species_json(sections['#DEFFIX'], fixed=True)

    """
    Generate MICM species JSON from KPP #DEFVAR section
    """
    defvar_json = micm_species_json(sections['#DEFVAR'])

    """
    Generate MICM equations JSON from KPP #EQUATIONS section
    """
    equations_json = micm_equation_json(sections['#EQUATIONS'])

    """
    Assemble MICM JSON
    """
    micm_species_json = {'camp-data': deffix_json + defvar_json}
    micm_species_json_str = json.dumps(micm_species_json, indent=4)
    logging.info('____ MICM species ____')
    logging.info(micm_species_json_str)
    print('\n')

    micm_mechanism_json = {'camp-data': {'reactions': equations_json}}
    micm_mechanism_json_str = json.dumps(micm_mechanism_json, indent=4)
    logging.info('____ MICM reactions ____')
    logging.info(micm_mechanism_json_str)
    print('\n')

