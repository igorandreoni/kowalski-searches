'''
import pdb
Query Kowalski searching for transients
given a set of constraints.
'''

import requests

import numpy as np
from astropy.time import Time
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
import healpy as hp
from ligo.skymap import postprocess
import matplotlib

from penquins import Kowalski


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def print_query_params(args, ra_center, dec_center, jd_trigger):
    '''Print a summary of the query parameters'''

    print("#-----")
    print("Cone search parameters:")
    print(f"A list of {len(ra_center)} coordinate pairs will be explored")
    print(f"Search radius {args.radius} arcmin")
    if args.after_trigger or jd_trigger > 0:
        print(f"Only sources detected for the first time after \
            {Time(jd_trigger, format='jd').iso} will be considered")
    print(f"Minimum time between the first and last alert {args.min_days} days"
          )
    print(f"Maximum time between the first and last alert {args.max_days} days"
          )
    print(f"Query divided in {args.slices} slices")
    print("#-----")
    print(" ")

    return


def get_programidx(program_name, username, password):
    '''Given a marshal science program name, it returns its programidx'''

    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_programs.cgi',
                      auth=(username, password))
    programs = r.json()
    program_dict = {p['name']: p['programidx'] for p in programs}

    try:
        return program_dict[program_name]
    except KeyError:
        print(f'The user {username} does not have access to the program \
              {program_name}')
        return None


def get_candidates_growth_marshal(program_name, username, password):
    '''Query the GROWTH db for the science programs'''

    programidx = get_programidx(program_name, username, password)
    if programidx is None:
        return None
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi',
                      auth=(username, password),
                      data={'programidx': str(programidx)})
    sources = r.json()
    sources_out = []
    for s in sources:
        coords = SkyCoord(ra=s['ra']*u.deg, dec=s['dec']*u.deg, frame='icrs')
        sources_out.append({"name": s['name'], "ra": coords.ra,
                            "dec": coords.dec,
                            "classification": s['classification'],
                            "redshift": s['redshift'],
                            "creation_date": s['creationdate']})

    return sources_out


def check_clu_transients(sources_kowalski, clu_sources):
    '''Check if the selected sources are present in the
    CLU science program.  If so, print out the relevant information.'''

    sources_in_clu = []
    sources_not_in_clu = []
    list_clu_sources = list(s['name'] for s in clu_sources)

    for source in sources_kowalski:
        if source['name'] in list_clu_sources:
            clu_source = clu_sources[np.where(
                                              np.array(list_clu_sources)
                                              == source['name'])[0][0]]
            try:
                for k in clu_source.keys():
                    print(f"{k}: {clu_source[k]}")
                sources_in_clu.append(source['name'])
            except (KeyError, ValueError, AttributeError):
                print("No CLU sources?")
        else:
            sources_not_in_clu.append(source['name'])
    print("Summary:")
    print(f"Sources saved in CLU: {sources_in_clu}")
    print(f"Sources not saved in CLU: {sources_not_in_clu}")


def read_skymap(skymap_filename):
    '''Read the healpix skymap'''

    hpx, header = hp.read_map(skymap_filename, h=True, verbose=False)
    header_dict = {}
    for x in header:
        header_dict[x[0]] = x[1]

    return hpx, header_dict


def tesselation_spiral(FOV, scale=0.80):
    FOV = np.pi*FOV*FOV*scale

    area_of_sphere = 4*np.pi*(180/np.pi)**2
    n = int(np.ceil(area_of_sphere/FOV))
    print("Using %d points to tile the sphere..."%n)

    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(n)
    z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)
    radius = np.sqrt(1 - z * z)

    points = np.zeros((n, 3))
    points[:, 0] = radius * np.cos(theta)
    points[:, 1] = radius * np.sin(theta)
    points[:, 2] = z

    ra, dec = hp.pixelfunc.vec2ang(points, lonlat=True)

    return ra, dec


def angular_distance(ra1, dec1, ra2, dec2):

    delt_lon = (ra1 - ra2)*np.pi/180.
    delt_lat = (dec1 - dec2)*np.pi/180.
    dist = 2.0*np.arcsin(np.sqrt(np.sin(delt_lat/2.0)**2 +
                                 np.cos(dec1 * np.pi/180.) *
                                 np.cos(dec2 * np.pi/180.) *
                                 np.sin(delt_lon/2.0)**2))

    return dist/np.pi*180.


def do_getfields(healpix, FOV=60/3600.0, ra=None, dec=None, radius=None,
                 level=None):

    ras, decs = tesselation_spiral(FOV, scale=0.80)
    if (ra is not None) and (dec is not None):
        dist = angular_distance(ras, decs, ra, dec)
        idx = np.where(dist <= radius)[0]
        ras, decs = ras[idx], decs[idx]
    elif (level is not None):
        cls = 100 * postprocess.find_greedy_credible_levels(healpix)
        paths = postprocess.contour(cls, [level], degrees=True, simplify=True)
        paths = paths[0]

        pts = np.vstack((ras, decs)).T
        idx = np.zeros((len(ras)))
        for path in paths:
            polygon = matplotlib.path.Path(path)
            check = polygon.contains_points(pts)
            check = list(map(int, check))
            idx = np.maximum(idx, check)
        idx = np.where(idx == 1)[0]
        ras, decs = ras[idx], decs[idx]

    return ras, decs


def select_sources_in_contour(sources, skymap, level=90):
    """Check that the selected sources lie within a given integrated
    probability of the skymap usinng the pixels directly"""

    skymap_prob = skymap.flat_2d
    sort_idx = np.argsort(skymap_prob)[::-1]
    csm = np.empty(len(skymap_prob))
    csm[sort_idx] = np.cumsum(skymap_prob[sort_idx])
    ipix_keep = sort_idx[np.where(csm <= level/100.)[0]]
    nside = hp.pixelfunc.get_nside(skymap_prob)
    sources_contour = list(s for s in sources
                           if ("ra" in s)
                           and (hp.ang2pix(
                                           nside,
                                           0.5 * np.pi -
                                           np.deg2rad(s["dec"].value),
                                           np.deg2rad(s["ra"].value)
                                           ) in ipix_keep
                                )
                           )

    return sources_contour


def query_kowalski(username, password, ra_center, dec_center, radius,
                   jd_trigger, min_days, max_days, slices, ndethist_min,
                   within_days, after_trigger=True):
    '''Query kowalski and apply the selection criteria'''

    k = Kowalski(username=username, password=password, verbose=False)
    # Initialize a set for the results
    set_objectId_all = set([])
    slices = slices + 1
    for slice_lim, i in zip(
                           np.linspace(0, len(ra_center), slices)[:-1],
                           np.arange(len(np.linspace(0, len(ra_center),
                                                     slices)[:-1]))):
        try:
            ra_center_slice = ra_center[int(slice_lim):
                                        int(np.linspace(0, len(ra_center),
                                                        slices)[:-1][i+1])]
            dec_center_slice = dec_center[int(slice_lim):
                                          int(np.linspace(0, len(dec_center),
                                                          slices)[:-1][i+1])]
        except IndexError:
            ra_center_slice = ra_center[int(slice_lim):]
            dec_center_slice = dec_center[int(slice_lim):]
        coords_arr = []
        for ra, dec in zip(ra_center_slice, dec_center_slice):
            try:
                # Remove points too far south for ZTF.
                # Say, keep only Dec>-40 deg to be conservative
                if dec < -40.:
                    continue
                coords = SkyCoord(ra=float(ra)*u.deg, dec=float(dec)*u.deg)
                coords_arr.append((coords.ra.deg, coords.dec.deg))
            except ValueError:
                print("Problems with the galaxy coordinates?")
                continue

        # Correct the minimum number of detections
        ndethist_min_corrected = int(ndethist_min - 1)

        # Correct the jd_trigger if the user specifies to query
        # also before the trigger
        if after_trigger is False:
            jd_trigger = 0
        try:
            print(f"slice: {int(slice_lim)}:{int(np.linspace(0,len(ra_center),slices)[:-1][i+1])}")
        except:
            print(f"slice: {int(slice_lim)}:{int(len(ra_center))}")
        q = {"query_type": "cone_search",
             "object_coordinates": {
                                    "radec": f"{coords_arr}",
                                    "cone_search_radius": f"{radius}",
                                    "cone_search_unit": "arcmin"
                                    },
             "catalogs": {
                          "ZTF_alerts": {
                                         "filter": {
                                                    "candidate.jd":
                                                    {'$gt': jd_trigger},
                                                    "candidate.drb":
                                                    {'$gt': 0.5},
                                                    "candidate.ndethist":
                                                    {'$gt':
                                                     ndethist_min_corrected},
                                                    "candidate.jdstarthist":
                                                    {'$gt': jd_trigger,
                                                     '$lt': jd_trigger +
                                                     within_days}
                                                     },
                                         "projection": {
                                                         "objectId": 1,
                                                         "candidate.rcid": 1,
                                                         "candidate.ra": 1,
                                                         "candidate.dec": 1,
                                                         "candidate.jd": 1,
                                                         "candidate.ndethist":
                                                         1,
                                                         "candidate.jdstarthist": 1,
                                                         "candidate.jdendhist":
                                                         1,
                                                         "candidate.jdendhist":
                                                         1,
                                                         "candidate.magpsf": 1,
                                                         "candidate.sigmapsf":
                                                         1,
                                                         "candidate.fid": 1,
                                                         "candidate.programid":
                                                         1,
                                                         "candidate.isdiffpos":
                                                         1,
                                                         "candidate.ndethist":
                                                         1,
                                                         "candidate.ssdistnr":
                                                         1,
                                                         "candidate.rb": 1,
                                                         "candidate.drb": 1,
                                                         "candidate.distpsnr1": 1,
                                                         "candidate.sgscore1": 1,
                                                         "candidate.srmag1": 1,
                                                         "candidate.distpsnr2": 1,
                                                         "candidate.sgscore2": 1,
                                                         "candidate.srmag2": 1,
                                                         "candidate.distpsnr3": 1,
                                                         "candidate.sgscore3": 1,
                                                         "candidate.srmag3": 1
                                                         }
                                           }
                           },
             "kwargs": {"hint": "gw01"}
             }

        # Perform the query
        r = k.query(query=q)
        print('Search completed for this slice.')

        objectId_list = []
        with_neg_sub = []
        old = []
        out_of_time_window = []
        stellar_list = []

        try:
            if r['result_data']['ZTF_alerts'] == []:
                continue
            keys_list = list(r['result_data']['ZTF_alerts'].keys())
        except AttributeError:
            print("Error in the keys list?? Check 'r' ")
            # Try one more time
            print("Trying to query again the same slice")
            try:
                r = k.query(query=q)
                keys_list = list(r['result_data']['ZTF_alerts'].keys())
            except AttributeError:
                print("The query failed again, skipping slice..")
                continue

        for key in keys_list:
            all_info = r['result_data']['ZTF_alerts'][key]

            for info in all_info:
                if info['objectId'] in old:
                    continue
                if info['objectId'] in stellar_list:
                    continue
                if np.abs(info['candidate']['ssdistnr']) < 10:
                    continue
                if info['candidate']['isdiffpos'] in ['f', 0]:
                    with_neg_sub.append(info['objectId'])
                if (info['candidate']['jdendhist'] -
                info['candidate']['jdstarthist']) < min_days:
                    continue
                if (info['candidate']['jdendhist'] -
                info['candidate']['jdstarthist']) > max_days:
                    old.append(info['objectId'])
                if (info['candidate']['jdstarthist'] -
                jd_trigger) > within_days:
                    old.append(info['objectId'])
                if after_trigger is True:
                    if (info['candidate']['jdendhist'] -
                    jd_trigger) > max_days:
                        out_of_time_window.append(info['objectId'])
                else:
                    if (info['candidate']['jdendhist'] -
                    info['candidate']['jdstarthist']) > max_days:
                        out_of_time_window.append(info['objectId'])
                try:
                    if (np.abs(info['candidate']['distpsnr1']) < 2. and
                    info['candidate']['sgscore1'] >= 0.50):
                        stellar_list.append(info['objectId'])
                except (KeyError, ValueError):
                    pass
                try:
                    if (np.abs(info['candidate']['distpsnr1']) < 15. and
                    info['candidate']['srmag1'] < 15. and
                    info['candidate']['sgscore1'] >= 0.5):
                        continue
                except (KeyError, ValueError):
                    pass
                try:
                    if (np.abs(info['candidate']['distpsnr2']) < 15. and
                    info['candidate']['srmag2'] < 15. and
                    info['candidate']['sgscore2'] >= 0.5):
                        continue
                except (KeyError, ValueError):
                    pass
                try:
                    if (np.abs(info['candidate']['distpsnr3']) < 15. and
                    info['candidate']['srmag3'] < 15. and
                    info['candidate']['sgscore3'] >= 0.5):
                        continue
                except (KeyError, ValueError):
                    pass

                objectId_list.append(info['objectId'])

        set_objectId = set(objectId_list)

        # Remove those objects with negative subtraction
        for n in set(with_neg_sub):
            try:
                set_objectId.remove(n)
            except (ValueError, KeyError):
                pass

        # Remove stellar objects
        for n in set(stellar_list):
            try:
                set_objectId.remove(n)
            except (ValueError, KeyError):
                pass

        # Remove those objects considered old
        for n in set(old):
            try:
                set_objectId.remove(n)
            except (ValueError, KeyError):
                pass

        # Remove those objects whole alerts go bejond jd_trigger+max_days
        for n in set(out_of_time_window):
            try:
                set_objectId.remove(n)
            except (ValueError, KeyError):
                pass
        print(set_objectId)

        set_objectId_all = set_objectId_all | set_objectId
        print("Cumulative:", set_objectId_all)

    return set_objectId_all


def select_sources_in_level(sources, skymap_filename, level=90):
    hpx, header = hp.read_map(skymap_filename, h=True, verbose=False)
    i = np.flipud(np.argsort(hpx))
    sorted_credible_levels = np.cumsum(hpx[i])
    credible_levels = np.empty_like(sorted_credible_levels)
    credible_levels[i] = sorted_credible_levels
    npix = len(hpx)
    nside = hp.npix2nside(npix)
    sources_within = list(s for s in sources
                          if (credible_levels[hp.ang2pix(
                                                        nside,
                                                        0.5 * np.pi -
                                                        np.deg2rad(s["dec"]),
                                                        np.deg2rad(s["ra"]))
                                              ] <= level/100.
                              )
                          )

    return sources_within


def query_kowalski_coords(username, password, names):
    '''Query kowalski to get the coordinates of given ZTF sources.'''

    names = list(names)
    k = Kowalski(username=username, password=password, verbose=False)

    q = {"query_type": "find",
         "query": {
                   "catalog": "ZTF_alerts",
                   "filter": {"objectId": {"$in": names}},
                   "projection": {"_id": 0,
                                  "candid": 1,
                                  "objectId": 1,
                                  "candidate.ra": 1,
                                  "candidate.dec": 1
                                  },
                   }
         }
    results_all = k.query(query=q)
    results = results_all['result_data']['query_result']
    sources = []
    for n in names:
        source = {}
        source["name"] = n
        source["ra"] = list(r["candidate"]["ra"] for r in results if
                            r["objectId"] == n)[0]
        source["dec"] = list(r["candidate"]["dec"] for r in results if
                             r["objectId"] == n)[0]
        source["candid"] = list(r["candid"] for r in results if
                                r["objectId"] == n)[0]
        sources.append(source)

    return sources


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query kowalski to find \
counterpart candidates.')
    parser.add_argument('--skymap', dest='skymap_filename', type=str,
                        required=False,
                        help='Skymap filename', default=None)
    parser.add_argument('--level', dest='level', type=float, required=False,
                        help='Enclosed probability', default=90)
    parser.add_argument('--fov', dest='fov', type=float, required=False,
                        help='Field of view of each cone (radius, in arcmin)',
                        default=60)
    parser.add_argument('--ra-center', dest='ra_center', nargs='+',
                        required=False,
                        help='Right ascension of the center (list, in deg)')
    parser.add_argument('--dec-center', dest='dec_center', nargs='+',
                        required=False,
                        help='Declination of the center (list, in deg)')
    parser.add_argument('--radius', dest='radius', type=float, required=False,
                        help='Search radius (min), by default radius = fov',
                        default=None)
    parser.add_argument('--after-trigger', dest='after_trigger', type=str2bool,
                        required=False,
                        help='Query only alerts whose first detection \
                        occurred after a certain date. If this boolean \
                        value is True, then --jd-trigger must be given.',
                        default=True)
    parser.add_argument('--jd-trigger', dest='jd_trigger', type=float,
                        required=False,
                        help='Julian Day of the trigger', default=-1.)
    parser.add_argument('--min-days', dest='min_days', type=float,
                        required=False,
                        help='Minimum time (days) between the first \
                        and last alert', default=0.)
    parser.add_argument('--max-days', dest='max_days', type=float,
                        required=False,
                        help='Maximum time (days) between the first \
                        and last alert', default=10000.)
    parser.add_argument('--within-days', dest='within_days', type=float,
                        required=False,
                        help='Maximum time (days) between the jd-trigger \
                        and the first alert', default=1000.)
    parser.add_argument('--ndethist', dest='ndethist_min', type=int,
                        required=False,
                        help='Minimum number of detections', default=2)
    parser.add_argument('--slices', dest='slices', type=int, required=False,
                        help='Number (integer) of slices in which the query \
                        will be devided', default=10)
    parser.add_argument('--out', dest='out', type=str, required=False,
                        help='Output filename', default='results.txt')
    parser.add_argument('--ingest', dest='ingest', type=str2bool,
                        required=False,
                        help='Use with caution! Ingest the candidates to the \
                        scanning page on the GROWTH marshal',
                        default=False)
    parser.add_argument('--ingest-program',
                        dest='ingest_program', type=str, required=False,
                        help='Ingest the candidates to the scanning page of \
                        the given GROWTH marshal program',
                        default='Electromagnetic Counterparts to Gravitational Waves')

    args = parser.parse_args()

    # Read the skymap and create the tessellation
    if args.skymap_filename is not None:
        healpix, header = read_skymap(args.skymap_filename)
        ra_center, dec_center = do_getfields(healpix, FOV=args.fov/60.,
                                             level=args.level)
    else:
        ra_center, dec_center = args.ra_center, args.dec_center

    if args.radius is None:
        args.radius = args.fov
    # Pre-select coordinates above Dec = -30
    ra_center = ra_center[np.where(dec_center > -30)]
    dec_center = dec_center[np.where(dec_center > -30)]

    if args.jd_trigger == -1:
        jd_trigger = Time(header['MJD-OBS'], format='mjd').jd
    else:
        jd_trigger = args.jd_trigger

    # Print a summary of the query input
    print_query_params(args, ra_center, dec_center, jd_trigger)
    if args.after_trigger is True:
        print("JD trigger:", jd_trigger)
        print("DATE trigger:", Time(jd_trigger, format='jd').iso)
    # Read the secrets
    secrets = ascii.read('secrets.csv', format='csv')

    username = secrets['kowalski_user'][0]
    password = secrets['kowalski_pwd'][0]

    # Query kowalski
    sources_kowalski = query_kowalski(username, password, ra_center,
                                      dec_center, args.radius, jd_trigger,
                                      args.min_days, args.max_days,
                                      args.slices, args.ndethist_min,
                                      args.within_days,
                                      after_trigger=args.after_trigger)

    # Check which sources are actually inside the contour
    sources = query_kowalski_coords(username, password, sources_kowalski)
    sources_within = select_sources_in_level(sources,
                                             args.skymap_filename,
                                             args.level)

    # Print results to an output text file
    with open(args.out, 'a') as f:
        f.write(f"{args} \n")
        for s in sources_within:
            f.write(f"{s['name']}, {s['ra']}, {s['dec']} \n")

    # GROWTH marshal credentials
    username_marshal = secrets['marshal_user'][0]
    password_marshal = secrets['marshal_pwd'][0]

    # Ingest into the marshal scanning page
    if args.ingest:
        tbl_ingested = ascii.read('ingested.csv')
        sources_ingestion = list({"name": s["name"], "candid": s["candid"]}
                                 for s in sources_within
                                 if not (s["name"] in
                                 list(tbl_ingested["name"])))
        # Check if the sources were already ingested under another program name
        for s in sources_within:
            if s["name"] in list(x["name"] for x in sources_ingestion):
                continue
            else:
                if not (args.ingest_program in
                        list(tbl_ingested["program"][tbl_ingested["name"] ==
                             s["name"]])):
                    print(f"Re-ingesting {s['name']} under \
                          {args.ingest_program}")
                    sources_ingestion.append({"name": s["name"],
                                              "candid": s["candid"]})
        if sources_ingestion == []:
            print("No sources to ingest")
        else:
            programidx = get_programidx(args.ingest_program, username_marshal,
                                        password_marshal)
            with open("ingested.csv", 'a') as ingested_file:
                for s in sources_ingestion:
                    r = requests.post("http://skipper.caltech.edu:8080/cgi-bin/growth/ingest_avro_id.cgi",
                                      auth=(username_marshal,
                                            password_marshal),
                                      data={'programidx': programidx,
                                            'avroid': s['candid']})
                    try:
                        r.raise_for_status()
                        print(f"Successfully ingested {s['name']} in \
                              {args.ingest_program}")
                        out_string = f"{s['name']}, {s['candid']}, {args.ingest_program}\n"
                        ingested_file.write(out_string)
                    except requests.exceptions.HTTPError as e:
                        print("Error ingesting {s['name']}\
                               in {args.ingest_program}")
                        print(e)

    # Check the CLU science program on the Marshal
    program_name = 'Census of the Local Universe'
    clu_sources = get_candidates_growth_marshal(program_name,
                                                username_marshal,
                                                password_marshal)

    # For each transient check if it is present in the CLU science program
    check_clu_transients(sources_within, clu_sources)

    print("Done.")
