from astropy.table import Table, hstack
from gentle.data import Data


def get_nearby_galaxies(args):
    """
    Interactive entry point to get nearby galaxies.
    """

    with Data(args.egis_path, args.leda_path) as data:
        data.set_log(args.verbose, args.log_path)

        galaxies = data.nearby_galaxies(
            args.galaxy,
            args.search_radius,
            args.search_distance,
            args.ned,
            args.leda,
            args.angular_search
        )

        data.log(f"""Started running `radius` command with parameters:
  - Galaxy:          {args.galaxy}\n
  - Radius:          {args.search_radius}\n
  - Distance:        {args.search_distance}\n
  - NED:             {args.ned}\n
  - Leda:            {args.leda}\n
  - Angular search:  {args.angular_search}
        """)

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
            else:
                print(results)
