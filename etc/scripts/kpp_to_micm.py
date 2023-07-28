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
        kpp_dir (str): KPP directory

    Returns
    """

    kpp_suffixes = ['.kpp', '.spc', '.eqn', '.def']

    for suffix in kpp_suffixes:
        files = glob(os.path.join(kpp_dir, '*' + suffix))
        logging.debug(files)


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

    read_kpp_config(args.kpp_dir)

