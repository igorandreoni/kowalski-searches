'''
Query Kowalski and apply the CLU filter.
'''

import requests
from datetime import timedelta

import numpy as np
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii

from penquins import Kowalski


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


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


def check_clu_transients(sources_kowalski, clu_sources, program):
    '''Check if the selected sources are present in a
     science program.  If so, print out the relevant information.'''

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
                print(f"No {program} sources?")
        else:
            sources_not_in_clu.append(source['name'])
    print("Summary:")
    print(f"Sources saved in {program}: {sources_in_clu}")
    print(f"Sources not saved in {program}: {sources_not_in_clu}")


def angular_distance(ra1, dec1, ra2, dec2):

    delt_lon = (ra1 - ra2)*np.pi/180.
    delt_lat = (dec1 - dec2)*np.pi/180.
    dist = 2.0*np.arcsin(np.sqrt(np.sin(delt_lat/2.0)**2 +
                                 np.cos(dec1 * np.pi/180.) *
                                 np.cos(dec2 * np.pi/180.) *
                                 np.sin(delt_lon/2.0)**2))

    return dist/np.pi*180.


def get_kowalski_index_info(catalog):
    '''Get information about indexes to favor queries
    of the specified catalog.'''

    q = {"query_type": "info",
         "query": {
             "command": "index_info",
             "catalog": catalog
         }
         }
    r = k.query(query=q)
    print(r)


def check_history(list_sources, radius=1.):
    '''Query ZTF and ZUDS alerts and select
    only those sources without a negative detection.
    '''

    # Get the coordinates of the candidates
    sources = query_kowalski_coords(username, password, list_sources, catalog='ZUDS_alerts')
    coords_arr = list((c['ra'], c['dec']) for c in sources)

    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "cone_search",
         "object_coordinates": {
                                "radec": f"{coords_arr}",
                                "cone_search_radius": f"{radius}",
                                "cone_search_unit": "arcsec"
                                },
         "catalogs": {
                      "ZTF_alerts": {
                                     "filter": {'candidate.isdiffpos': {'$in': ['f', 0]}
                                                },
                                     "projection": {"objectId": 1,
                                                    "candidate.ra": 1,
                                                    "candidate.dec": 1
                                                    }
                                     },
                      "ZUDS_alerts": {
                                     "filter": {'candidate.isdiffpos': {'$in': ['f', 0]}
                                                },
                                     "projection": {"objectId": 1,
                                                    "candidate.ra": 1,
                                                    "candidate.dec": 1
                                                    }
                                     }
                     }
        }

    r = k.query(query=q)

    with_neg_sub = []
    for s, i in zip(sources, r['result_data']['ZTF_alerts'].keys()):
        if len(r['result_data']['ZTF_alerts'][i]) != 0 or len(r['result_data']['ZUDS_alerts'][i]) != 0:
            with_neg_sub.append(s['name'])

    selected_list = list(s['name'] for s in sources if not(s['name'] in with_neg_sub))

    return selected_list


def match_kowalski_clu(username, password, list_in,
                       catalog='ZUDS_alerts_aux'):
    '''Query kowalski and apply the CLU filter'''

    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "find",
         "query": {
                   "catalog": catalog,
                   "filter": {"_id": {'$in': list(list_in)}},
                   "projection": {"objectId": 1,
                                  "cross_matches.CLU_20190625": 1}
                   }
         }

    r = k.query(query=q)

    query_results = r['result_data']['query_result']
    list_out = list(t['_id'] for t in query_results if len(t['cross_matches']['CLU_20190625']) != 0)

    return list_out


def query_kowalski_alerts(username, password, date_start, date_end,
                       catalog='ZUDS_alerts', min_days=None, starthist=None):
    '''Query alerts with kowalski and apply the selection criteria'''

    k = Kowalski(username=username, password=password, verbose=False)

    # Initialize a set for the results
    set_objectId_all = set([])
    print(date_start.jd, date_end.jd, starthist.jd)
    q = {"query_type": "find",
         "query": {"catalog": catalog,
                   "filter": {"candidate.jd":
                              {'$gt': date_start.jd,
                              '$lt': date_end.jd},
                              "candidate.drb": {'$gt': 0.6},
                              "classifications.braai": {'$gt': 0.6},
                              "candidate.jdstarthist_single": {'$gt': starthist.jd},
                              "candidate.fwhm": {'$gt': 0.5, '$lt': 8},
                              },
                   "projection": {"objectId": 1,
                                  "candid": 1,
                                  "candidate.rcid": 1,
                                  "candidate.ra": 1,
                                  "candidate.dec": 1,
                                  "candidate.jd": 1,
                                  "candidate.ndethist": 1,
                                  "candidate.jdstarthist_single": 1,
                                  "candidate.jdstarthist_stack": 1,
                                  "candidate.jdendhist_single": 1,
                                  "candidate.jdendhist_stack": 1,
                                  "candidate.magpsf": 1,
                                  "candidate.sigmapsf": 1,
                                  "candidate.fid": 1,
                                  "candidate.programid": 1,
                                  "candidate.isdiffpos": 1,
                                  "candidate.ndethist": 1,
                                  "candidate.ssdistnr": 1,
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
                                  "candidate.srmag3": 1,
                                  "candidate.fwhm": 1,
                                  "candidate.lstype1": 1,
                                  "candidate.lszspec1": 1,
                                  "candidate.lsz1": 1,
                                  "candidate.lszphotl681": 1,
                                  "candidate.alert_type": 1
                                  }
                   },
             "kwargs": {} # {"limit": 3}
             }

        # Perform the query
    r = k.query(query=q)
    result = r['result_data']['query_result']

    set_names = set(list(c['objectId'] for c in r['result_data']['query_result']))
    print(f"There are {len(set_names)} sources found")

    # Match with CLU
    list_clu = match_kowalski_clu(username, password, set_names)
    print(f"{len(list_clu)} sources were found matched with CLU galaxies")

    # Stricted selection
    reject = []
    done = []
    done_names = []

    for info in r['result_data']['query_result']:
        if info['objectId'] == 'ZUDS20esmwf':
            import pdb
            pdb.set_trace()
        if not (info['objectId'] in list_clu) or info['objectId'] in done_names or info['objectId'] in reject:
            continue
        # If single alerts (not stack), check that there is enough
        # time separation between first and last observation
        try:
            if info['candidate']['alert_type'] == 'single':
                if info['candidate']['jdendhist_single'] - info['candidate']['jdstarthist_single'] > min_days:
                    pass
                else:
                    reject.append(info['objectId'])
                    continue
            else:
                pass
        except (KeyError, ValueError, TypeError):
            pass
             
        try:
            if (np.abs(info['candidate']['distpsnr1']) < 1. and
            info['candidate']['sgscore1'] >= 0.60):
                reject.append(info['objectId'])
                continue
        except (KeyError, ValueError, TypeError):
            pass
        #try:
        #    if (np.abs(info['candidate']['distpsnr1']) < 1. and
        #    info['candidate']['lstype1'] == 'PSF'):
        #        reject.append(info['objectId'])
        #        continue
        #except (KeyError, ValueError, TypeError):
        #    pass
        try:
            if (np.abs(info['candidate']['distpsnr1']) < 1. and
            info['candidate']['lsz1'] < 999. and
            info['candidate']['lszspec1'] > 0.1):
                reject.append(info['objectId'])
                continue
        except (KeyError, ValueError, TypeError):
            pass
        try:
            if (np.abs(info['candidate']['distpsnr1']) < 1. and
            info['candidate']['lsz1'] < 21. and
            info['candidate']['lszphotl681'] > 0.1):
                reject.append(info['objectId'])
                continue
        except (KeyError, ValueError, TypeError):
            pass
        try:
            if (np.abs(info['candidate']['distpsnr1']) < 15. and
            info['candidate']['srmag1'] < 15. and
            info['candidate']['srmag1'] > 0. and
            info['candidate']['sgscore1'] >= 0.49):
                reject.append(info['objectId'])
                continue
        except (KeyError, ValueError, TypeError):
            pass
        try:
            if (np.abs(info['candidate']['distpsnr2']) < 15. and
            info['candidate']['srmag2'] < 15. and
            info['candidate']['srmag2'] > 0. and
            info['candidate']['sgscore2'] >= 0.49):
                reject.append(info['objectId'])
                continue
        except (KeyError, ValueError, TypeError):
            pass
        try:
            if (np.abs(info['candidate']['distpsnr3']) < 15. and
            info['candidate']['srmag3'] < 15. and
            info['candidate']['srmag3'] > 0. and
            info['candidate']['sgscore3'] >= 0.49):
                reject.append(info['objectId'])
                continue
        except (KeyError, ValueError, TypeError):
            pass
        done_names.append(info['objectId'])
        done.append((info['objectId'], info['candid']))

    # Check that no source was kept among the rejected ones
    checked = list(d for d in done if not (d[0] in reject))
    checked_names = list(c[0] for c in checked)

    print(f"{len(done)} sources survived stellarity and bright sources cuts")

    # Check history for negative subtractions in ZUDS and ZTF alerts
    list_selected = check_history(checked_names)
    print(f"{len(list_selected)} sources have no historical negative detections")

    sources = list({'name': c[0], 'candid': c[1]} for c in checked if c[0] in list_selected)
    return sources


def query_kowalski_coords(username, password, names, catalog='ZTF_alerts'):
    '''Query kowalski to get the coordinates of given ZTF sources.'''

    names = list(names)
    k = Kowalski(username=username, password=password, verbose=False)

    q = {"query_type": "find",
         "query": {
                   "catalog": catalog,
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

    parser = argparse.ArgumentParser(description='Find CLU transients.')
    parser.add_argument('--date-start', dest='date_start', type=str,
                        required=False,
                        help="Start date of the query, in ISO format. \
                        Example: '2017-08-17 12:41:04.4'", default=None)
    parser.add_argument('--date-end', dest='date_end', type=str,
                        required=False,
                        help="End date of the query, in ISO format. \
                        Example: '2017-08-18 12:00:00.0'", default=None)
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
                        default='ZUDS and CLU')
    parser.add_argument('--starthist',
                        dest='starthist', type=str, required=False,
                        help="Earliest date of the first detection; \
                        Example: '2017-08-18 12:00:00.0'",
                        default=None)
    parser.add_argument('--min-days',
                        dest='min_days', type=str, required=False,
                        help="Minimum time (days) between first and \
                        last detection",
                        default=0.01)

    args = parser.parse_args()

    # Define start and end date as Time objects
    if args.date_start is None:
        date_start = Time.now() - timedelta(days=1)
    else:
        try:
            date_start = Time(args.date_start, format='iso')
        except  ValueError:
            print("Invalid start date. It must be a string in ISO format.")
            print("Example: '2017-08-17 12:41:04.4'")
            exit()
    if args.date_end is None:
        date_end = Time.now()
    else:
        try:
            date_end = Time(args.date_end, format='iso')
        except ValueError:
            print("Invalid end date. It must be a string in ISO format.")
            print("Example: '2017-08-18 12:00:00.0'")
            exit()
    print(f"Query start date: {date_start}")
    print(f"Query End date: {date_end}")

    try:
        if args.starthist is None:
            starthist = Time.now() - timedelta(days=30) 
        else:
            starthist = Time(args.starthist, format='iso')
    except ValueError:
        print("Invalid starthist date. It must be a string in ISO format.")
        print("Example: '2017-08-15 12:00:00.0'")


    # Print a summary of the query input

    # Read the secrets
    secrets = ascii.read('../kowalski/secrets.csv', format='csv')

    username = secrets['kowalski_user'][0]
    password = secrets['kowalski_pwd'][0]

    sources = query_kowalski_alerts(username, password,
                                    date_start, date_end,
                                    catalog='ZUDS_alerts',
                                    min_days=args.min_days,
                                    starthist=starthist)

    # Print results
    print(sources)
    print(list(s['name'] for s in sources))

    '''
    # Print results to an output text file
    with open(args.out, 'a') as f:
        f.write(f"{args} \n")
        for s in sources_within:
            f.write(f"{s['name']}, {s['ra']}, {s['dec']} \n")
    '''

    # GROWTH marshal credentials
    username_marshal = secrets['marshal_user'][0]
    password_marshal = secrets['marshal_pwd'][0]

    # Ingest into the marshal scanning page
    if args.ingest:
        """
        tbl_ingested = ascii.read('ingested.csv')
        sources_ingestion = list({"name": s["name"], "candid": s["candid"]}
                                 for s in sources
                                 if not (s["name"] in
                                 list(tbl_ingested["name"])))
        """
        # Check if the sources were already ingested under another program name
        for s in sources:
            """
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
            sources_ingestion.append({"name": s["name"],
                                      "candid": s["candid"]})
            """
        sources_ingestion = sources
        if sources_ingestion == []:
            print("No sources to ingest")
        else:
            with open("ingested.csv", 'a') as ingested_file:
                for s in sources_ingestion:
                    if args.ingest_program == 'ZUDS and CLU':
                        cgi_script = "http://skipper.caltech.edu:8080/cgi-bin/growth/ingest_zuds_id.cgi"
                        r = requests.post(cgi_script,
                                          auth=(username_marshal,
                                                password_marshal),
                                          data={'programid': 87,
                                                'candid': s['candid'],
                                                'save': 'false'
                                                })
                    else:
                        programidx = get_programidx(args.ingest_program,
                                                    username_marshal,
                                                    password_marshal)
                        cgi_script = "http://skipper.caltech.edu:8080/cgi-bin/growth/ingest_avro_id.cgi"
                        r = requests.post(cgi_script,
                                          auth=(username_marshal,
                                                password_marshal),
                                          data={'programidx': programidx,
                                                'avroid': s['candid'],
                                                'save': 'false'})
                    try:
                        r.raise_for_status()
                        print(f"Successfully ingested {s['name']} in {args.ingest_program}")
                        out_string = f"{s['name']}, {s['candid']}, {args.ingest_program}\n"
                        ingested_file.write(out_string)
                    except requests.exceptions.HTTPError as e:
                        print("Error ingesting {s['name']}\
                               in {args.ingest_program}")
                        print(e)

    exit()

    # Check the CLU science program on the Marshal
    program_name = args.ingest_program
    clu_sources = get_candidates_growth_marshal(program_name,
                                                username_marshal,
                                                password_marshal)

    # For each transient check if it is present in the CLU science program
    #for s in sources:
    #    print(s['name'])
    check_clu_transients(sources, clu_sources, program_name)

    print("Done.")
