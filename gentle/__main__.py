"""GENTLE: Galactic ENvironment and TidaL Estimation

Use NED and HyperLeda databases to estimate the galactic
environment of EGIS galaxies."""


import argparse
import sys

from gentle.utility import get_nearby_galaxies


def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--egis-path',
        type = str,
        help = 'Path to EGIS .fit table.',
        default = None
    )
    
    parser.add_argument(
        '--leda-path',
        type = str,
        help = 'Path to HyperLeda .tsv table.',
        default = None
    )
    
    parser.add_argument(
        '--verbose',
        type = str,
        help = 'Whether to log while running.',
        action = argparse.BooleanOptionalAction,
        default = False
    )
    
    parser.add_argument(
        '--log-path',
        type = str,
        help = 'Where to save log.',
        default = None
    )

    subparsers = parser.add_subparsers()
    
    radius_parser = subparsers.add_parser(
        'radius',
        help = 'Extract all galaxies within a certain radius.'
    )
    radius_parser.set_defaults(func=get_nearby_galaxies)

    radius_parser.add_argument(
        'galaxy',
        type = str,
        help = 'The target galaxy.'
    )

    radius_parser.add_argument(
        'search_radius',
        type = float,
        help = 'The search radius around the galaxy (Mpc).'
    )

    radius_parser.add_argument(
        'search_distance',
        type = float,
        help = 'The spatial distance around the galaxy (Mpc).'
    )

    radius_parser.add_argument(
        '--output-path',
        type = str,
        help = 'Where to write the data. Will print instead if blank.',
        default = ''
    )

    radius_parser.add_argument(
        '--ned',
        help = 'Whether to search in the NED database.',
        action = argparse.BooleanOptionalAction,
        default = True
    )

    radius_parser.add_argument(
        '--leda',
        help = 'Whether to search in the HyperLeda database.',
        action = argparse.BooleanOptionalAction,
        default = False
    )

    radius_parser.add_argument(
        '--names',
        help = 'Whether to only display galaxy names.',
        action = argparse.BooleanOptionalAction,
        default = True
    )

    radius_parser.add_argument(
        '--angular-search',
        help = 'Whether the search radius should be in arcsec instead of Mpc.',
        action = argparse.BooleanOptionalAction,
        default = False
    )
    
    args = parser.parse_args()

    if 'func' in args:
        args.func(args)
    else:
        parser.print_help(sys.stderr)


if __name__ == '__main__':
    main()
