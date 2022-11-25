from astropy.table import Table, hstack

from gentle.data import Data
from gentle.redshift import redshift_to_distance


def get_nearby_galaxies(args):
    """
    Interactive entry point to get nearby galaxies.
    """

    with Data(args.egis_path, args.leda_path) as data:
        data.set_log(args.verbose, args.log_path)

        data.log(f"""Started running `radius` command with parameters:
  - Galaxy:          {args.galaxy}
  - Radius:          {args.search_radius}
  - Distance:        {args.search_distance}
  - NED:             {args.ned}
  - Leda:            {args.leda}
  - Angular search:  {args.angular_search}""")

        galaxies = data.nearby_galaxies(
            args.galaxy,
            args.search_radius,
            args.search_distance,
            args.ned,
            args.leda,
            args.angular_search
        )

        results = None

        if args.names:
            ned_galaxies = galaxies['ned']
            leda_galaxies = galaxies['leda']

            if ned_galaxies is not None:
                ned_names = ned_galaxies['Object Name']
                ned_galaxies.rename_column('Object Name', 'NED')
            else:
                ned_names = Table()

            if leda_galaxies is not None:
                leda_names = leda_galaxies['objname']
                leda_galaxies.rename_column('objname', 'HyperLeda')
            else:
                leda_names = Table()

            results = hstack([ned_names, leda_names])

        if isinstance(results, Table):
            if args.output_path:
                results.write(args.output_path, overwrite=True)
                data.log(f'Wrote results to {args.output_path}.')
            else:
                data.log(f'Results:\n {results}')


def search_galaxy(args):
    """
    Interactive entry point to search for specific galaxies.
    """

    with Data(args.egis_path, args.leda_path) as data:
        data.set_log(args.verbose, args.log_path)

        data.log(f"""Started running `search` command with parameters:
  - Galaxy:  {args.galaxy}
  - NED:     {args.ned}
  - Leda:    {args.leda}""")

        results = data.search_galaxy(args.galaxy, args.ned, args.leda, args.field)

        data.log(f"""Results:
NED
---
{results['ned']}
HyperLeda
---------
{results['leda']}""")


def convert_redshift(args):
    """
    Interactive entry point to convert redshift to distance.
    """

    with Data(None, None) as data:
        data.set_log(args.verbose, args.log_path)

        data.log(f"""Started running `redshift` command with parameters:
  - z:  {args.z}""")

        distance = redshift_to_distance(args.z)

        data.log(f'Results:\n {distance}')


def convert_size_distance(args):
    """
    Interactive entry point to convert `v` to distance.
    """

    with Data(None, None) as data:
        data.set_log(args.verbose, args.log_path)

        data.log(f"""Started running `DA` command with parameters:
  - v:     {args.v} """)

        distance = data.angular_size_distance(args.v)

        data.log(f'Results:\n {distance}')
