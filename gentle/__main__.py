"""GENTLE: Galactic ENvironment and TidaL Estimation

Use NED and HyperLeda databases to estimate the galactic
environment of EGIS galaxies."""


import argparse
import sys

from gentle.utility import get_nearby_galaxies, search_galaxy, convert_redshift, convert_size_distance


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

    search_parser = subparsers.add_parser(
        'search-name',
        help = 'Find data for a particular galaxy.'
    )
    search_parser.set_defaults(func=search_galaxy)

    search_parser.add_argument(
        'galaxy',
        type = str,
        help = 'The target galaxy.'
    )

    search_parser.add_argument(
        '--field',
        type = str,
        help = 'The target field.',
        default = None
    )

    search_parser.add_argument(
        '--ned',
        help = 'Whether to search in the NED database.',
        action = argparse.BooleanOptionalAction,
        default = True
    )

    search_parser.add_argument(
        '--leda',
        help = 'Whether to search in the HyperLeda database.',
        action = argparse.BooleanOptionalAction,
        default = True
    )

    redshift_parser = subparsers.add_parser(
        'redshift',
        help = 'Convert redshift to distance in Mpc.'
    )
    redshift_parser.set_defaults(func=convert_redshift)

    redshift_parser.add_argument(
        'z',
        type = float,
        help = 'The redshift.'
    )

    size_distance_parser = subparsers.add_parser(
        'DA',
        help = 'Convert `v` to angular size distance.'
    )
    size_distance_parser.set_defaults(func=convert_size_distance)

    size_distance_parser.add_argument(
        'v',
        type = float,
        help = 'The radial heliocentric velocity.'
    )
    
    args = parser.parse_args()

    if 'func' in args:
        args.func(args)
    else:
        parser.print_help(sys.stderr)


if __name__ == '__main__':
    main()
