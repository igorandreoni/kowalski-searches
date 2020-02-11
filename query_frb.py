'''
Query Kowalski searching for counterparts to FRBs
Author: Igor Andreoni
'''

import numpy as np
import json
from collections import OrderedDict
import pdb

from astropy.time import Time
from astropy.table import Table, unique
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord

from ztfquery import query
from penquins import Kowalski


def query_metadata(ra, dec, username, password,
                   start_jd=None, end_jd=None,
                   out_csv=None):
    """Use ZTFquery to get more reliable upper limits"""

    zquery = query.ZTFQuery()
    if start_jd == None and end_jd==None:
        zquery.load_metadata(kind='sci',
                             radec=[str(ra.deg), str(dec.deg)],
                             size=0.003,
                             auth=[username, password])
    else:
        if start_jd!=None and end_jd==None:
            sql_query='obsjd>'+repr(start_jd)
        elif start_jd==None and end_jd!=None:
            sql_query='obsjd<'+repr(end_jd)
        elif start_jd!=None and end_jd!=None:
            sql_query='obsjd<'+repr(end_jd)+'+AND+'+'obsjd>'+repr(start_jd)
        zquery.load_metadata(kind='sci',
                             radec=[str(ra.deg), str(dec.deg)],
                             size=0.003,
                             sql_query=sql_query,
                             auth=[username, password])
    out = zquery.metatable
    final_out = out.sort_values(by=['obsjd'])
    if out_csv is not None:
        final_out.to_csv(out_csv)

    return final_out


def create_tbl_lc(light_curves, outfile):
    """Create a table with the light curves
    and write a CSV output file"""

    # fid -> filter
    filters = {'1': 'g', '2': 'r', '3': 'i'}

    tbl = Table([[], [], [], [], [], [], [], [], [], [], [], [], [], [], [],
                 [], [], [], []],
                names=('name', 'ra', 'dec', 'jd', 'magpsf', 'sigmapsf',
                       'filter', 'magzpsci', 'magzpsciunc',
                       'programid', 'field', 'rcid', 'pid',
                       'sgscore1', 'sgscore2', 'sgscore3',
                       'distpsnr1', 'distpsnr2', 'distpsnr3'),
                dtype=('S12', 'double', 'double', 'double',
                       'f', 'f', 'S', 'f', 'f', 'i', 'i', 'i', 'int_',
                       'f', 'f', 'f', 'f', 'f', 'f'))

    for l in light_curves:
        magzpsci = l["candidate"].get("magzpsci")
        magzpsciunc = l["candidate"].get("magzpsciunc")
        try:
            row = [l["objectId"], l["candidate"]["ra"], l["candidate"]["dec"],
               l["candidate"]["jd"], l["candidate"]["magpsf"],
               l["candidate"]["sigmapsf"], filters[str(l["candidate"]["fid"])],
               magzpsci, magzpsciunc,
               l["candidate"]["programid"], l["candidate"]["field"],
               l["candidate"]["rcid"], l["candidate"]["pid"],
               l["candidate"]["sgscore1"], l["candidate"]["sgscore2"],
               l["candidate"]["sgscore3"], l["candidate"]["distpsnr1"],
               l["candidate"]["distpsnr2"], l["candidate"]["distpsnr3"]]
        except KeyError:
            row = [l["objectId"], l["candidate"]["ra"], l["candidate"]["dec"],
               l["candidate"]["jd"], l["candidate"]["magpsf"],
               l["candidate"]["sigmapsf"], filters[str(l["candidate"]["fid"])],
               magzpsci, magzpsciunc,
               l["candidate"]["programid"], l["candidate"]["field"],
               l["candidate"]["rcid"], l["candidate"]["pid"], np.nan,
               np.nan, np.nan, np.nan, np.nan, np.nan]
        tbl.add_row(row)
    # Remove exact duplicates
    tbl = unique(tbl)
    tbl.sort("jd")
    if args.out_lc is not None:
        tbl.write(args.out_lc, format='csv', overwrite=True)

    return tbl


def get_lightcurve_alerts_aux(username, password, list_names):
    """Query the light curve for a list of candidates"""

    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "find",
         "query": {
                   "catalog": "ZTF_alerts_aux",
                   "filter": {
                              '_id': {'$in': list(list_names)}
                              },
                   "projection": {}
                   },
         "kwargs": {"hint": "_id_"}
         }
    r = k.query(query=q)
    if r['result_data']['query_result'] == []:
        print("No candidates to be checked?")
        return None
    out = []
    for l in r['result_data']['query_result']:
        with_det = list({'objectId': l['_id'], 'candidate': s}
                        for s in l['prv_candidates']
                        if 'magpsf' in s.keys())
        out = out + with_det

    return out


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
                                  "candidate.magzpsciunc": 1,
                                  "candidate.sgscore1": 1,
                                  "candidate.sgscore2": 1,
                                  "candidate.sgscore3": 1,
                                  "candidate.distpsnr1": 1,
                                  "candidate.distpsnr2": 1,
                                  "candidate.distpsnr3": 1,
                                  "candidate.field": 1,
                                  "candidate.rcid": 1,
                                  "candidate.pid": 1
                                  }
                       },
         "kwargs": {"hint": "objectId_1"}
         }

    r = k.query(query=q)
    if r['result_data']['query_result'] == []:
        print("No candidates to be checked?")
        return None

    return r['result_data']['query_result']


def query_kowalski_frb(args, t):
    """Query kowalski with cone searches centered at given locations"""

    # Prepare a dictionary for each source
    dict_sources = {}
    for s in t:
        if args.frb_names is not None and not(s['frb_name'] in args.frb_names):
            continue
        try:
            coords=SkyCoord(ra=s["rop_raj"], dec=s["rop_decj"],
                            unit=(u.hourangle, u.deg), frame='icrs')
        except ValueError:
            pdb.set_trace()
        id_ra = f"{str(coords.ra.deg).replace('.','_')}"
        id_dec = f"{str(coords.dec.deg).replace('.','_')}"
        id_coords = f"({id_ra}, {id_dec})"
        date = Time(s['utc'].replace('/','-'), format='iso')
        dict_sources[s['frb_name']] = {
                                       'ra': coords.ra.deg,
                                       'dec': coords.dec.deg,
                                       'id_coords': id_coords,
                                       'jd': date.jd,
                                       'candidates': []
                                       }
    # Check that there is at least one source
    if len(dict_sources.keys()) == 0:
        print("No FRBs correspond to the given input.")
        if args.frb_names is not None:
            print(f"No FRB among {args.frb_names} are present in {args.cat_name}")
        return None

    # coords_arr.append((coords.ra.deg,coords.dec.deg))
    coords_arr = list((dict_sources[k]['ra'],
                        dict_sources[k]['dec'])
                       for k in dict_sources.keys())

    k = Kowalski(username=username, password=password, verbose=False)

    q = {"query_type": "cone_search",
         "object_coordinates": {
                                "radec": f"{coords_arr}",
                                "cone_search_radius": args.search_radius,
                                "cone_search_unit": "arcmin"
                                },
         "catalogs": {
                      "ZTF_alerts": {
                                     "filter": {
                                                "candidate.drb": {'$gt': 0.5},
                                                "candidate.ndethist": {'$gte': args.ndethist},
                                                "classifications.braai": {'$gt': 0.5},
                                                "candidate.ssdistnr": {'$gt': 10},
                                                "candidate.magpsf": {'$gt': 10}
                                                },
                                     "projection": {
                                                    "objectId": 1,
                                                    "candidate.rcid": 1,
                                                    "candidate.drb": 1,
                                                    "candidate.ra": 1,
                                                    "candidate.dec": 1,
                                                    "candidate.jd": 1,
                                                    "candidate.magpsf": 1,
                                                    "candidate.sigmapsf": 1,
                                                    "candidate.fid": 1,
                                                    "candidate.sgscore1": 1,
                                                    "candidate.distpsnr1": 1,
                                                    "candidate.sgscore2": 1,
                                                    "candidate.distpsnr2": 1,
                                                    "candidate.sgscore3": 1,
                                                    "candidate.distpsnr3": 1,
                                                    "candidate.ssdistnr": 1,
                                                    "candidate.isdiffpos": 1
                                                    }
                                     }
                      },
          "kwargs": {"hint": "gw01"}
         }

    r = k.query(query=q)

    for idcoords in r['result_data']['ZTF_alerts'].keys():
        #Identify 'candid' for all relevant candidates
        objectId_list = []
        with_neg_sub = []
        stellar_list = []

        # No sources        
        if len(r['result_data']['ZTF_alerts'][idcoords]) == 0:
            key = list(k for k in dict_sources.keys()
                       if dict_sources[k]['id_coords']==idcoords)[0]
            dict_sources[key]['candidates'] = []
            print(f"No candidates for {key}")
            continue

        for i in np.arange(len(r['result_data']['ZTF_alerts'][idcoords])):
            info = r['result_data']['ZTF_alerts'][idcoords][i]
            if info['objectId'] in stellar_list or (info['objectId'] in with_neg_sub):
                continue
            if info['candidate']['isdiffpos'] in ['f',0]:
                with_neg_sub.append(info['objectId'])
            try:
                if (np.abs(info['candidate']['distpsnr1']) < 2.
                and info['candidate']['sgscore1'] >= 0.5):
                    stellar_list.append(info['objectId'])
            except:
                pass            
            try:
                if (np.abs(info['candidate']['distpsnr1']) < 15. and
                           info['candidate']['srmag1'] < 15. and
                           info['candidate']['srmag1'] > 0. and
                           info['candidate']['sgscore1'] >= 0.5):
                    continue
            except:
                pass
            try:
                if (np.abs(info['candidate']['distpsnr2']) < 15. and
                           info['candidate']['srmag2'] < 15. and
                           info['candidate']['srmag2'] > 0. and
                           info['candidate']['sgscore2'] >= 0.5):
                    continue
            except:
                pass
            try:
                if (np.abs(info['candidate']['distpsnr3']) < 15. and
                           info['candidate']['srmag3'] < 15. and
                           info['candidate']['srmag3'] > 0. and
                           info['candidate']['sgscore3'] >= 0.5):
                    continue
            except:
                pass
            objectId_list.append(info['objectId'])
        set_objectId = set(objectId_list)

        # Remove objects with negative subtraction
        if args.reject_neg:
            for n in set(with_neg_sub):
                try:
                    set_objectId.remove(n)
                except:
                    pass

        # Remove stellar objects
        for n in set(stellar_list):
            try:
                set_objectId.remove(n)
            except:
                pass

        # Add the list of ZTF candidates to the FRB list
        key = list(k for k in dict_sources.keys() 
                   if dict_sources[k]['id_coords']==idcoords)[0]
        dict_sources[key]['candidates'] = list(set(set_objectId))
        tot_sources = len(r['result_data']['ZTF_alerts'][idcoords])
        print(f"{len(set_objectId)}/{tot_sources} candidates selected for {key}")

    return dict_sources


def plot_results_frb(sources, t, t_lc, args):
    """Fetch the light curves of the candidates and
       plot the results"""

    filters = ['g', 'r', 'i']
    filters_id = {'1': 'g', '2': 'r', '3': 'i'}
    colors = {'g': 'g', 'r': 'r', 'i': 'y'}

    for frb in sources.keys():
        if len(sources[frb]['candidates']) == 0:
            continue
        # A different plot for each candidate
        jd0 = sources[frb]['jd']
        for cand in set(sources[frb]['candidates']):      
            plt.clf()
            plt.subplot(1, 1, 1)
            t_cand = t_lc[t_lc['name']==cand]
            for f in filters:
                tf = t_cand[t_cand['filter'] == f]
                tf["jd"] = tf["jd"] - jd0
                mag = np.array(tf["magpsf"])
                magerr = np.array(tf["sigmapsf"])
                plt.errorbar(np.array(tf["jd"]),
                             np.array(tf["magpsf"]),
		             fmt=colors[f]+'o',
                             yerr=np.array(tf["sigmapsf"]),
		             markeredgecolor='k', markersize=8,
                             label=f)
                # Upper limits
                if args.use_metadata is True:
                    username_ztfquery = secrets['ztfquery_user'][0]
                    password_ztfquery = secrets['ztfquery_pwd'][0]
                    coords = SkyCoord(ra=np.mean(t_cand[t_cand['name']==cand]['ra']*u.deg),
                                      dec=np.mean(t_cand[t_cand['name']==cand]['dec']*u.deg))
                    metadata = query_metadata(coords.ra, coords.dec,
                                              username_ztfquery,
                                              password_ztfquery,
                                              start_jd=None,
                                              end_jd=None,
                                              out_csv=None)
                    t_ul = Table([[],[],[],[],[],[],[],[]],
                                 names=('jd', 'magpsf', 'sigmapsf', 'filter',
                                        'snr', 'ul', 'seeing', 'programid'),
                                 dtype=('double','f','f','S','f','f','f','int')
                                 )
                    for j, ml, fid, s, pid in zip(metadata['obsjd'],
                                                  metadata['maglimit'],
                                                  metadata['fid'],
                                                  metadata['seeing'],
                                                  metadata['pid']):
                        #if not (pid in t['pid']):
                        if not (j in t_cand['jd']):
                            new_row = [j, 99.9, 99.9, filters_id[str(fid)],
                                       np.nan, ml, s, 0]
                            t_ul.add_row(new_row)
                    for f in filters:
                        tf_ul = t_ul[t_ul['filter'] == f]
                        if len(tf_ul) > 0:
                            tf_ul["jd"] = tf_ul["jd"] - jd0
                            plt.plot(np.array(tf_ul["jd"]),
                                     np.array(tf_ul["ul"]),
                                     colors[f]+'v',
                                     markeredgecolor=colors[f],
                                     markerfacecolor='w')
                            plt.plot([],[], 'kv', label='UL')
            # Plot the FRB detection time
            plt.plot([0, 0], [22, 15], 'b--', label=frb)
            # Legend
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys())
            # Labels
            plt.ylabel("Apparent magnitude [AB]")
            plt.xlabel(f"Days since {Time(jd0, format='jd').iso}")
            plt.title(f"{frb}; {cand}")
            # Invert the y axis
            plt.gca().invert_yaxis()
            # Save the plot automatically?
            if args.saveplot:
                plt.savefig(f"lc_{frb}_{cand}.png")
            plt.show()


def get_index_info(catalog):
    """List which indexes are available on Kowalski to query a catalog
       more quickly"""
    q = {"query_type": "info",
         "query": {
             "command": "index_info",
             "catalog": catalog
         }
         }
    k = Kowalski(username=username, password=password, verbose=False)
    r = k.query(query=q)
    indexes = r['result_data']['query_result']
    for ii, (kk, vv) in enumerate(indexes.items()):
        print(f'index #{ii+1}: "{kk}"\n{vv["key"]}\n')


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Query kowalski.')
    parser.add_argument('--cat', dest='cat_name', type=str, required=True,
                        help='CSV file downloaded from the FRBcat page',
                        default=None)
    parser.add_argument('--frb', dest='frb_names', nargs='+', required=False,
                        help='Names of specific FRBs to check', default=None)
    parser.add_argument('--out', dest='out', type=str,
                        required=False,
                        help='output JSON file name; example: results.json',
                        default=None)
    parser.add_argument('--out-lc', dest='out_lc', type=str,
                        required=False,
                        help='output CSV file name; example: lc.csv',
                        default=None)
    parser.add_argument('--r', dest='search_radius', type=float,
                        required=False,
                        help='Cone search radius in arcmin (default=15)',
                        default=15)
    parser.add_argument('--ndethist', dest='ndethist', type=int,
                        required=False,
                        help='Minimum number of detections (default=2)',
                        default=2)
    parser.add_argument('--p', dest='plot', type=str2bool,
                        required=False,
                        help='Plot the results? (boolean)',
                        default=True)
    parser.add_argument('--um', dest='use_metadata', type=str2bool,
                        required=False,
                        help='Plot upper limits using ztfquery metadata (boolean)',
                        default=True)
    parser.add_argument('--sp', dest='saveplot', type=str2bool,
                        required=False,
                        help='Save the plot of the results? (boolean)',
                        default=False)
    parser.add_argument('--reject-neg', dest='reject_neg', type=str2bool,
                        required=False,
                        help='Reject candidates with negative detections? (boolean)',
                        default=True)
    args = parser.parse_args()

    # Read the table of FRBs
    t=ascii.read(args.cat_name, format='csv')

    # Correct the first column name if necessary
    if '\ufeff"frb_name"' in t.colnames:
        t.rename_column('\ufeff"frb_name"', "frb_name")

    # Read the secrets
    secrets = ascii.read('secrets.csv', format = 'csv')

    username = secrets['kowalski_user'][0]
    password = secrets['kowalski_pwd'][0]

    if args.use_metadata is True:
        username_ztfquery = secrets['ztfquery_user'][0]
        password_ztfquery = secrets['ztfquery_pwd'][0]

    sources = query_kowalski_frb(args, t)

    # Write the results into an output JSON file
    if args.out is not None:
        with open(args.out, 'w') as j:
            json.dump(sources, j)

    # Get the light curves
    list_names = []
    for k in sources.keys():
        list_names += sources[k]['candidates']

    if len(list_names) == 0:
        print("No ZTF sources selected. Exiting..")
        exit()
    light_curves_alerts = get_lightcurve_alerts(username, password,
                                         list_names)

    # Add prv_candidates photometry to the light curve
    light_curves_aux = get_lightcurve_alerts_aux(username, password,
                                         list_names)

    light_curves = light_curves_alerts + light_curves_aux

    # Create a table and output CSV file
    t_lc = create_tbl_lc(light_curves, args)

    # Plot the results
    if args.plot is True:
        plot_results_frb(sources, t, t_lc, args)
