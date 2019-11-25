from astropy.io import ascii
from astropy.table import Table

from penquins import Kowalski


def get_lightcurve_alerts(username, password, list_names):
    """Query the light curve for a list of candidates"""

    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "find",
         "query": {
                   "catalog": "ZTF_alerts",
                   "filter": {
                              'objectId': {'$in': list(list_names)}
                              },
                   "projection": {
                                  "objectId": 1,
                                  "candidate.jd": 1,
                                  "candidate.ra": 1,
                                  "candidate.dec": 1,
                                  "candidate.magpsf": 1,
                                  "candidate.fid": 1,
                                  "candidate.sigmapsf": 1,
                                  "candidate.programid": 1,
                                  "candidate.magzpsci": 1,
                                  "candidate.magzpsciunc": 1
                                  }
                       },
         "kwargs": {"hint": "objectId_1"}
         }

    r = k.query(query=q)
    if r['result_data']['query_result'] == []:
        print("No candidates to be checked?")
        return None

    return r['result_data']['query_result']


def create_tbl_lc(light_curves, outfile):
    """Create a table with the light curves
    and write a CSV output file"""

    # fid -> filter
    filters = {'1': 'g', '2': 'r', '3': 'i'}

    tbl = Table([[], [], [], [], [], [], [], [], [], []],
                names=('name', 'ra', 'dec', 'jd', 'magpsf', 'sigmapsf',
                       'filter', 'magzpsci', 'magzpsciunc',
                       'programid'),
                dtype=('S12', 'double', 'double', 'double',
                       'f', 'f', 'S', 'f', 'f', 'i'))

    for l in light_curves:
        magzpsci = l["candidate"].get("magzpsci")
        magzpsciunc = l["candidate"].get("magzpsciunc")
        row = [l["objectId"], l["candidate"]["ra"], l["candidate"]["dec"],
               l["candidate"]["jd"], l["candidate"]["magpsf"],
               l["candidate"]["sigmapsf"], filters[str(l["candidate"]["fid"])],
               magzpsci, magzpsciunc,
               l["candidate"]["programid"]]
        tbl.add_row(row)
    tbl.sort("jd")
    tbl.write(outfile, format='csv', overwrite=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query kowalski to fetch \
transient light curves.')
    parser.add_argument('--n', dest='names', nargs='+', required=False,
                        help='Names of the ZTF candidates; if given, \
the coordinates will be queried from kowalski', default=None)
    parser.add_argument('--out', dest='out', type=str, required=False,
                        help='Output filename', default='lightcurves.csv')

    args = parser.parse_args()

    if args.names is None:
        print("No input candidates. Please use --n and provide ZTF names")
        exit()

    # Read the secrets
    secrets = ascii.read('secrets.csv', format='csv')
    username_kowalski = secrets['kowalski_user'][0]
    password_kowalski = secrets['kowalski_pwd'][0]

    # Get the light curves
    light_curves = get_lightcurve_alerts(username_kowalski, password_kowalski,
                                         args.names)
    # Create a table and output CSV file
    create_tbl_lc(light_curves, args.out)
