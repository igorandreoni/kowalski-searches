'''
Query Kowalski searching for transients
given a set of constraints.

'''

import json
import requests
import pdb

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import ascii
from astropy.table import Table

from penquins import Kowalski

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def print_query_params(args, ra_center, dec_center):
    '''Print a summary of the query parameters'''

    print("#-----")
    print("Cone search parameters:")
    print(f"A list of {len(ra_center)} coordinate pairs will be explored")
    print(f"Search radius {args.radius} arcmin")
    if args.after_trigger or args.jd_trigger > 0:
        print(f"Only sources detected for the first time \
after {Time(args.jd_trigger, format='jd').iso} will be considered")
    print(f"Minimum time between the first and last alert {args.min_days} days")
    print(f"Maximum time between the first and last alert {args.max_days} days")    
    print(f"Query divided in {args.slices} slices")
    print("#-----")
    print(" ")

    return


def get_programidx(program_name, username, password):
    ''' Given a marshal science program name, it returns its programidx'''

    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_programs.cgi',
                      auth=(username, password))
    programs=json.loads(r.text)
    program_dict={p['name']:p['programidx'] for i,p in enumerate(programs)}

    try:
        return program_dict[program_name]
    except KeyError:
        print(f'The user {username} does not have access to \
the program {program_name}')
        return None


def get_candidates_growth_marshal(program_name, username, password):
    ''' Query the GROWTH db for the science programs '''

    programidx=get_programidx(program_name, username, password)
    if programidx==None:
        return None
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi', \
        auth=(username, password), data={'programidx':str(programidx)})
    sources=json.loads(r.text)
    sources_out=[]
    for s in sources:
            coords=SkyCoord(ra=s['ra']*u.deg, dec=s['dec']*u.deg, frame='icrs')
            sources_out.append({"name":s['name'],
                                "ra":coords.ra, "dec":coords.dec,
	                        "classification":s['classification'],
                                "redshift":s['redshift'],
                                "creation_date":s['creationdate']})

    return sources_out


def check_clu_transients(sources_kowalski, clu_sources):
    '''Check if the selected sources are present in the 
    CLU science program.  If so, print out the relevant information.'''

    sources_in_clu = []
    sources_not_in_clu = []
    list_clu_sources = list(s['name'] for s in clu_sources)

    for source in sources_kowalski:
        print("-------")
        if source in list_clu_sources:
            clu_source = clu_sources[np.where(np.array(list_clu_sources) == source)[0][0]]
            try:
                for k in clu_source.keys():
                    print(f"{k}: {clu_source[k]}")
                sources_in_clu.append(source)
            except:
                pdb.set_trace()
        else:
            print(f"{source} was not saved in CLU")
            sources_not_in_clu.append(source)
        print("-------")
    print("Summary:")
    print(f"Sources saved in CLU: {sources_in_clu}")
    print(f"Sources not saved in CLU: {sources_not_in_clu}")

    return


def match_kowalski_clu(username, password, list_in,
                       catalog='ZTF_alerts_aux'):
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


def check_lightcurve_alerts(username, password, list_names, min_days, max_days):
    """Re-query light curve info for a list of candidates\
    and check that their full/updated duration is consistent\
    with the time limits provided"""

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
                                  "candidate.ndethist": 1,
                                  "candidate.jdstarthist": 1,
                                  "candidate.jdendhist": 1,
                                  "candidate.jdendhist": 1,
                                  "candidate.magpsf": 1,
                                  "candidate.sigmapsf": 1,
                                  "candidate.programid": 1,
                                  }
                       },
            "kwargs": {"hint": "objectId_1"}
             }

    r = k.query(query=q)
    if r['result_data']['query_result'] == []:
        print("No candidates to be checked?")
        return None

    old = []
    objectid_list = []
    for info in r['result_data']['query_result']:
        if info['objectId'] in old:
            continue
        if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) < min_days:
            continue
        if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) > max_days:
            old.append(info['objectId'])
        objectid_list.append(info['objectId'])
    clean_set = set(objectid_list)
    #Remove those objects considered old
    for n in set(old):
        try:
            clean_set.remove(n)
        except:
            do = 'do nothing'

    return clean_set



def match_file_clu(tbl, clu):
    '''Crossmatch the candidates with CLU'''

    names = set(list(n for n in list(tbl['name'])))

    ra = list(np.mean(tbl['ra'][tbl['name'] == name]) for name in list(names))
    dec = list(np.mean(tbl['dec'][tbl['name'] == name]) for name in list(names))
    coords = SkyCoord(ra=np.array(ra)*u.deg, dec=np.array(dec)*u.deg)
    # Read CLU
    clu_coords = SkyCoord(ra=clu['ra']*u.deg, dec=clu['dec']*u.deg)
    list_match = []
    for name, coord in zip(names, coords):
        sep = clu_coords.separation(coord)
        dist_kpc = clu['distmpc']*(10**3)*np.sin(sep)/np.cos(sep)
        condition0 = dist_kpc >= 0
        clu_match = clu[condition0]
        #sep_match = sep[condition0]
        dist_kpc = dist_kpc[condition0]
        condition = dist_kpc < 100
        clu_match = clu_match[condition]
        #sep_match = sep_match[condition]
        #dist_kpc = dist_kpc[condition]

        if len(clu_match) > 0:
            list_match.append(name)

    return list_match

def query_kowalski(username, password, list_fields, min_days, max_days,
                   ndethist_min, jd, jd_gap=50., clu=None):
    '''Query kowalski and apply the selection criteria'''

    k = Kowalski(username=username, password=password, verbose=False)

    # Correct the minimum number of detections
    ndethist_min_corrected = int(ndethist_min - 1)

    jd_start = jd
    jd_end = jd + jd_gap

    #Initialize a set for the results
    set_objectId_all = set([])
    for field in list_fields:
        set_objectId_field = set([])
        q = {"query_type": "find",
             "query": {
                       "catalog": "ZTF_alerts",      
                       "filter": {
                                  'candidate.jd': {'$gt': jd_start, '$lt': jd_end},
                                  'candidate.field': int(field),
                                  'candidate.rb': {'$gt': 0.2},
                                  #'candidate.drb': {'$gt': 0.5},
                                  'classifications.braai': {'$gt': 0.5},
                                  'candidate.ndethist': {'$gt': ndethist_min_corrected},
                                  #'candidate.ndethist': ndethist_min_corrected,
                                  'candidate.magpsf': {'$gt': 15}
                                  #'candidate.isdiffpos': 't'
                                   },
                       "projection": {
                                      "objectId": 1,
                                      "candidate.rcid": 1,
                                      "candidate.ra": 1,
                                      "candidate.dec": 1,
                                      "candidate.jd": 1,
                                      "candidate.ndethist": 1,
                                      "candidate.jdstarthist": 1,
                                      "candidate.jdendhist": 1,
                                      "candidate.jdendhist": 1,
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
                                      "candidate.srmag3": 1
                                       }
                       },
            "kwargs": {"hint": "jd_field_rb_drb_braai_ndethhist_magpsf_isdiffpos"}
             }

        #Perform the query
        r = k.query(query=q)
        print(f"Search completed for field {field}, \
{Time(jd, format='jd').iso} + {jd_gap} days.")

        try:
            if r['result_data']['query_result'] == []:
                print("No candidates")
                continue
        except KeyError:
            print(f"ERROR! jd={jd}, field={field}" ) 
            #pdb.set_trace()
            continue

        #Identify 'candid' for all relevant candidates
        objectId_list = []
        with_neg_sub = []
        old = []
        out_of_time_window = []
        stellar_list = []

        for info in r['result_data']['query_result']:    
            #if info['objectId'] == 'ZTF19abyfbii':
 	    #    pdb.set_trace()
            if info['objectId'] in old:
                continue
            if info['objectId'] in stellar_list:
                continue
            if np.abs(info['candidate']['ssdistnr']) < 10:
                continue
            try:
                if info['candidate']['drb'] < 0.5:
                    continue
            except KeyError:
                pass
            if info['candidate']['isdiffpos'] in ['f',0]:
                with_neg_sub.append(info['objectId'])
            if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) < min_days:
                continue
            if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) > max_days:
                old.append(info['objectId'])
            try:
                if (np.abs(info['candidate']['distpsnr1']) < 2. and info['candidate']['sgscore1'] >= 0.5):
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

        #Remove those objects with negative subtraction
        for n in set(with_neg_sub):
            try:
                set_objectId.remove(n)
            except:
                do = 'do nothing'

        #Remove stellar objects
        for n in set(stellar_list):
            try:
                set_objectId.remove(n)
            except:
                do = 'do nothing'

        #Remove those objects considered old
        for n in set(old):
            try:
                set_objectId.remove(n)
            except:
                do = 'do nothing'

        #Remove those objects whole alerts go bejond jd_trigger+max_days
        for n in set(out_of_time_window):
            try:
                set_objectId.remove(n)
            except:
                do = 'do nothing'

        if clu is not None and len(set_objectId) != 0:
            print(f"Field {field}, {len(set_objectId)} candidates before CLU cut")
            name_list = list(info['objectId'] for info 
                            in r['result_data']['query_result'] if info['objectId'] in set_objectId)
            ra_list = list(info['candidate']['ra'] for info 
                            in r['result_data']['query_result'] if info['objectId'] in set_objectId)
            dec_list = list(info['candidate']['dec'] for info 
                            in r['result_data']['query_result'] if info['objectId'] in set_objectId)
            
            # Create a table
            name_empty = ['']*len(name_list)
            num_empty = np.zeros(len(name_list))
            tt = Table([name_empty, num_empty, num_empty],
                       names=('name', 'ra', 'dec'))
            tt['name'] = name_list
            tt['ra'] = ra_list
            tt['dec'] = dec_list
            # Match with CLU
            list_clu = match_file_clu(tt, clu)
            set_objectId = set(list_clu)
            print(f"Field {field}, {len(set_objectId)} candidates after CLU cut")
        #print(set_objectId)
        set_objectId_all = set_objectId_all | set_objectId
        #print("Cumulative:", set_objectId_all)

        print("Field", field, len(set_objectId_all))

        '''
        print('----stats-----')
        print('Number of sources with negative sub: ', len(set(with_neg_sub)))
        print('Number of sources with only pos subtraction: ', len(set_objectId))
        print(f"Number of sources older than {max_days} days: {len(set(old))}, specifically {set(old)}")
        '''
    return set_objectId_all


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query kowalski.')
    parser.add_argument('--min-days', dest='min_days', type=float, required=False, \
    help='Minimum time (days) between the first and last alert', default = 0.)
    parser.add_argument('--max-days', dest='max_days', type=float, required=False, \
    help='Maximum time (days) between the first and last alert', default = 30.)
    parser.add_argument('--within-days', dest='within_days', type=float, required=False, \
    help='Maximum time (days) between the jd-trigger and the first alert', default = 1000.)
    parser.add_argument('--ndethist', dest='ndethist_min', type=int, required=False, \
    help='Minimum number of detections', default=2)
    parser.add_argument('--out', dest='out', type=str, required=False, \
    help='Output filename', default = 'results.txt')
    parser.add_argument('--clu', dest='clu_filename', type=str, required=False, \
    help='If assigned, crossmatched with the given CLU file', default=None)
    
    args = parser.parse_args()

    t = ascii.read('/data/ia/paper_kn_ZTF/selected_fields_ebv03.csv')
#    t01 = ascii.read('../myfolder/selected_fields.csv')

    # REMOVE!!!
    # Select only fields that do not have E(B-V) <0.1 but have E(B-V) <0.3 
    #list_fields = list(set(f for f in t['field'] if not (f in t01['field'])))
    list_fields = list(set(f for f in t['field'] if f > 156 and f < 1000))
    
    #Read the secrets
    secrets = ascii.read('/data/ia/kowalski_GW/secrets.csv', format = 'csv')

    username = secrets['kowalski_user'][0]
    password = secrets['kowalski_pwd'][0]

    # Iterate over JDs from July 1 2019 to November 15 2019
    date_start = Time('2018-03-01 00:00:00', format='iso')
    date_finish = Time('2020-03-30 00:00:00', format='iso')
    sources_kowalski_all = []
    jd_gap = 30.
    list_jd = np.linspace(date_start.jd, date_finish.jd,
                          int((date_finish.jd - date_start.jd)/jd_gap)+1)
    jd_gap = list_jd[1] - list_jd[0] + 1

    # Read the CLU catalog
    if args.clu_filename is not None:
        # Read the FITS file
        clu=Table.read(args.clu_filename)
        clu = clu[clu['distmpc'] > 5]
    else:
        clu = None 

    for jd in list_jd:    
        #Query kowalski
        sources_kowalski = query_kowalski(username, password, list_fields,
                                          args.min_days, args.max_days,
                                          args.ndethist_min,
                                          jd, jd_gap=jd_gap, clu=clu)

        sources_kowalski_all += list(sources_kowalski)
    sources_kowalski_all = set(sources_kowalski_all)

    # Check full light curve duration (alerts)
    print("Checking durations.....")
    clean_set = check_lightcurve_alerts(username, password,
                                        sources_kowalski_all,
                                        args.min_days, args.max_days)
    print("...Done.")
    print("Final set:")
    print(clean_set)
    #for s in clean_set:
    #    print(s)

    #Print results to an output text file
    with open(args.out, 'a') as f:
        f.write(f"#{args} \n")
        f.write("name \n")
        for n in clean_set:
            f.write(f"{n} \n")

    '''
    #Check the CLU science program on the Marshal
    username_marshal = secrets['marshal_user'][0]
    password_marshal= secrets['marshal_pwd'][0]
    
    program_name='Census of the Local Universe'
    clu_sources = get_candidates_growth_marshal(program_name,
                                                username_marshal,
                                                password_marshal)    

    #For each transient check if it is present in the CLU science program
    check_clu_transients(clean_set, clu_sources)
    '''
    print("Done.")
