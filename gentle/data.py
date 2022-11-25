"""
Abstraction over data processing.
"""

import datetime

from astroquery.ipac.ned import Ned
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord

import astropy.units as u

import numpy as np

from gentle.redshift import redshift_to_distance
from gentle.constants import *


class Data:
    """
    Abstraction over data processing.
    """

    def __init__(self, egis_path, leda_path):
        self.egis_path = egis_path
        self.leda_path = leda_path

        self.egis_data = None
        self.egis_hdul = None

        self.leda_data = None

        self.masks = set()

        self.log_enabled = False
        self.log_path = None
        self.log_file = None


    def __enter__(self):
        if self.egis_path is not None:
            self.egis_hdul = fits.open(self.egis_path)
            self.egis_data = self.egis_hdul[1].data

            self.log('Got data from EGIS.')
        else:
            self.egis_hdul = None

        if self.leda_path is not None:
            self.leda_data = ascii.read(
                self.leda_path,
                guess=False,
                format='tab'
            )

            self.log('Got data from HyperLeda.')

        return self


    def __exit__(self, exc_type, exc_value, exc_tb):
        self.log('Closing Data instance.')

        if self.egis_hdul is not None:
            self.egis_hdul.close()

        if self.log_file is not None:
            self.log_file.close()


    def set_log(self, enabled, path):
        """
        Set the state of the log.

        Paramaters
        ----------
        - enabled: bool -> whether to enable logging
        - path: str | None -> None to print, path to write to file
        """

        self.log_enabled = enabled
        self.log_path = path

        if path is not None:
            try:
                self.log_file = open(path, 'w')
            except FileNotFoundError:
                self.log_path = None

                self.log(f'Path invalid: {path}; will switch to printing mode.')


    def log(self, message):
        """
        Log messages as execution runs.

        Parameters
        ----------
        - message: str -> the message to be logged
        """

        log_message = f'[{datetime.datetime.now()}]: {message}\n'

        if self.log_enabled:
            if self.log_path is None:
                print(log_message, end='')
            else:
                self.log_file.write(log_message)


    def add_mask(self, mask):
        """
        For selecting rows that match mask.

        Parameters
        ----------
        - mask: (data) -> [bool]
        """

        self.masks.add(mask)


    def apply_masks(self, data):
        """
        Apply all masks on data.

        Parameters
        ----------
        - data: astropy.table.Table
        """

        for mask in self.masks:
            data = data[mask(data)]

        return data


    def angular_distance(self, dec_a, ra_a, dec_b, ra_b, unit_a, unit_b):
        """
        Get the angular distance (arcsec) between two objects.

        Parameters
        ----------
        - dec_a: float -> declension of object A
        - ra_a: float -> RA of object A
        - dec_b: float -> declension of object B
        - ra_b: float -> RA of object B
        """

        if unit_a == 'hms':
            coord_a = SkyCoord(
                dec=dec_a, ra=ra_a,
                frame='icrs',
                unit=(u.hourangle, u.degree)
            )
        else:
            coord_a = SkyCoord(
                dec=dec_a, ra=ra_a,
                frame='icrs',
                unit=u.degree
            )

        if unit_b == 'hms':
            coord_b = SkyCoord(
                dec=dec_b, ra=ra_b,
                frame='icrs',
                unit=(u.hourangle, u.degree)
            )
        else:
            coord_b = SkyCoord(
                dec=dec_b, ra=ra_b,
                frame='icrs',
                unit=u.degree
            )

        dec_a = coord_a.dec.radian
        ra_a = coord_a.ra.radian
        dec_b = coord_b.dec.radian
        ra_b = coord_b.ra.radian

        theta = np.arccos(
            np.sin(dec_a)*np.sin(dec_b) +\
            np.cos(dec_a)*np.cos(dec_b)*np.cos(ra_a - ra_b)
        )
        theta = np.degrees(theta) * 3600

        self.log(f'Got angular separation of {theta}.')

        return theta


    def angular_size_distance(self, cz):
        """
        Convert from cz to angular size distance (Mpc).

        Parameters
        ----------
        - cz: float -> heliocentric radial velocity (km/s)
        """

        D = cz/H0

        self.log(f'Got size distance of {D}.')

        return D


    def redshift_to_distance(self, z):
        """
        Wrapper for `redshift_to_distance`.

        Parameters
        ----------
        - z: float -> redshift
        """

        distance = redshift_to_distance(z)

        self.log(f'Got size distance of {distance}.')

        return distance


    def search_galaxy(self, galaxy, ned=True, leda=True, field=None):
        """
        Find a target galaxy in the NED and HyperLeda databases.

        Paramters
        ---------
        - galaxy: str -> the galaxy name
        - field: str -> the field to extract
        """

        tables = {
            'ned': None,
            'leda': None
        }

        if ned:
            obj = Ned.query_object(galaxy)

            if field is not None and field in obj.colnames:
                obj = obj[field]

            tables['ned'] = obj

            self.log('Got object from NED database.')

        if leda and self.leda_data is not None:
            obj = self.leda_data[self.leda_data['objname'] == galaxy]

            if field is not None and field in obj.colnames:
                obj = obj[field]

            tables['leda'] = obj

            self.log('Got object from HyperLeda database.')

        return tables


    def nearby_galaxies(self, galaxy, radius, distance, ned=True, leda=True, angular_search=False):
        """
        Find all galaxies within a cylindrical angular radius and distance.

        Parameters
        ----------
        - galaxy: str -> the galaxy name
        - radius: float -> the angular radius within which to find nearby
          galaxies (arcsec/Mpc; see `angular_search`)
        - distance: float -> the spatial distance within which to find
          nearby galaxies (Mpc)
        - ned: bool -> whether to include the NED database in the search
        - leda: bool -> whether to include the HyperLeda database in the
          search
        - angular_search: bool -> whether to use arcsec instead of Mpc
          for radius
        """

        tables = {
            'ned': None,
            'leda': None
        }

        ned_object = Ned.query_object(galaxy)
        object_position = (ned_object['DEC'], ned_object['RA'])

        self.log(f'Got position {object_position[0][0]} (DEC), {object_position[1][0]} (RA) from NED database.')

        object_velocity = ned_object['Velocity'][0]

        self.log(f'Got object velocity {object_velocity} from NED database.')

        object_distance = self.angular_size_distance(object_velocity)

        if not angular_search:
            radius = 206265 * radius/object_distance

            self.log(f'Switched to search angle {radius} arcsec.')

        if ned:
            self.log(f'Started getting nearby galaxies from NED database.')

            ned_table = Ned.query_region(galaxy, radius=radius * u.arcsec)

            angular_size_distance = np.vectorize(self.redshift_to_distance)

            ned_table = ned_table[np.abs(
                angular_size_distance(
                    ned_table['Redshift']
                ) - object_distance
            ) < distance]

            tables['ned'] = ned_table

            self.log(f'Finished getting nearby galaxies from NED database.')

        if leda and self.leda_path is not None:
            self.log(f'Started getting nearby galaxies from HyperLeda database.')

            leda_data = self.leda_data.copy()
            leda_table = leda_data.copy()
            leda_table.remove_rows(slice(0, len(leda_table)))

            leda_data['v'] = leda_data['v'].filled(np.nan)

            for galaxy in leda_data:
                cz = galaxy['v']

                if np.isnan(cz):
                    continue

                angular_distance = np.abs(self.angular_distance(
                    object_position[0][0],
                    object_position[1][0],
                    galaxy['de2000'],
                    galaxy['al2000'],
                    'degree', 'hms'
                ))

                if angular_distance > radius:
                    continue

                size_distance = np.abs(self.angular_size_distance(cz) - object_distance)

                if size_distance > distance:
                    continue

                leda_table.add_row(galaxy)

            if len(leda_table) > 0:
                tables['leda'] = leda_table

            self.log(f'Finished getting nearby galaxies from HyperLeda database.')

        return tables

