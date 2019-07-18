'''
Query Kowalski searching for transients
given a set of constraints.

'''


def print_query_params(args):
    '''Print a summary of the query parameters'''

    print("#-----")
    print("Cone search parameters:")
    print(f"Search radius {args.radius} arcmin")
    print(f"Minimum time between the first and last alert {args.min_days} days")
    print(f"Maximum time between the first and last alert {args.max_days} days")    
    print(f"CLU galaxies selected with distance between {args.min_dist}Mpc and {args.max_dist}Mpc, with Dec > {args.min_dec}") 
    print(f"Query divided in {args.slices} slices")
    print("#-----")
    print(" ")

    return


def get_programidx(program_name, username, password):
    ''' Given a marshal science program name, it returns its programidx'''

    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_programs.cgi', auth=(username, password))
    programs=json.loads(r.text)
    program_dict={p['name']:p['programidx'] for i,p in enumerate(programs)}

    try:
        return program_dict[program_name]
    except KeyError:
        print(f'The user {username} does not have access to the program {program_name}')
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
            sources_out.append({"name":s['name'], "ra":coords.ra, "dec":coords.dec, \
	        "classification":s['classification'], "redshift":s['redshift'], "creation_date":s['creationdate']})

    return sources_out


def query_kowalski_clu(username, password, clu):
    '''Query kowalski to get a table of CLU galaxies. '''

    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "general_search", 
        "query": "db['CLU_20180513'].find({},{'distmpc': 1})" 
        }
    r = k.query(query=q)
    pdb.set_trace()
    return


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


def query_kowalski(username, password, clu, args):
    '''Query kowalski and apply the selection criteria'''

    k = Kowalski(username=username, password=password, verbose=False)
    #Initialize a set for the results
    set_objectId_all = set([])

    for slice_lim,i in zip(np.linspace(0,len(clu),args.slices)[:-1], np.arange(len(np.linspace(0,len(clu),args.slices)[:-1]))):
        try:
            t=clu[int(slice_lim):int(np.linspace(0,len(clu),args.slices)[:-1][i+1])]
        except IndexError:
            t=clu[int(slice_lim):]
        coords_arr = []
        galaxy_names_arr = []
        for galaxy,ra, dec in zip(t["name"],t["ra"], t["dec"]):
            try:
                coords=SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
                coords_arr.append((coords.ra.deg,coords.dec.deg))
            except ValueError:
                print("Problems with the galaxy coordinates?")
                pdb.set_trace()
                continue
            galaxy_names_arr.append(galaxy)
        try: 
            print(f"slice: {int(slice_lim)}:{int(np.linspace(0,len(clu),args.slices)[:-1][i+1])}" )
        except:
            print(f"slice: {int(slice_lim)}:{int(len(clu))}" )
        q = {"query_type": "cone_search",
        "object_coordinates": {
             "radec": f"{coords_arr}", 
             "cone_search_radius": f"{args.radius}",
             "cone_search_unit": "arcmin"
         },
         "catalogs": {
             "ZTF_alerts": {
                 "filter": {"candidate.ndethist": {'$gt': 1}},
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
             }
         }
         }

        #Perform the query
        r = k.query(query=q)

        print('Search completed for this slice.')
#        #Dump the results in a json file
#        with open(f'results_clu25Mpc_1week_{i+1}.json', 'w') as j:
#            json.dump(r, j)

        #Identify 'candid' for all relevant candidates
        objectId_list = []
        with_neg_sub = []
        old = []
        stellar_list = []
        try:
            keys_list=list(r['result_data']['ZTF_alerts'].keys())
        except:
            print("Error in the keys list?? Check 'r' ")
            pdb.set_trace()
        for key in keys_list:
            all_info=r['result_data']['ZTF_alerts'][key]
            for info in all_info:
#                #Stop at a certain candidId for debugging
#                if info['objectId'] == 'ZTF19aanfkyc':
#                    pdb.set_trace()
                if info['objectId'] in old:
                    continue
                if info['objectId'] in stellar_list:
                    continue
                if info['candidate']['rb'] < 0.2:
                    continue
                try:
                    if info['candidate']['drb'] < 0.8:
                        continue
                except:
                    do = 'do nothing.'
                if np.abs(info['candidate']['ssdistnr']) < 10:
                    continue
                if info['candidate']['isdiffpos'] in ['f',0]:
                    with_neg_sub.append(info['objectId'])
                if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) < args.min_days:
                    continue
                if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) > args.max_days:
                    old.append(info['objectId'])
                try:
                    if (np.abs(info['candidate']['distpsnr1']) < 3. and info['candidate']['sgscore1'] > 0.0):
                        stellar_list.append(info['objectId'])
                except:
                    do = 'do nothing.'
                try:
                    if (np.abs(info['candidate']['distpsnr1']) < 15. and info['candidate']['srmag1'] < 15. and info['candidate']['sgscore1'] >= 0.5):
                        continue
                except:
                    do = 'do nothing.'
                try:
                    if (np.abs(info['candidate']['distpsnr2']) < 15. and info['candidate']['srmag2'] < 15. and info['candidate']['sgscore2'] >= 0.5):
                        continue
                except:
                    do = 'do nothing.'
                try:
                    if (np.abs(info['candidate']['distpsnr3']) < 15. and info['candidate']['srmag3'] < 15. and info['candidate']['sgscore3'] >= 0.5):
                        continue
                except:
                    do = 'do nothing.'

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

        print(set_objectId)

        set_objectId_all = set_objectId_all | set_objectId
        print("Cumulative:", set_objectId_all)

        '''
        print('----stats-----')
        print('Number of sources with negative sub: ', len(set(with_neg_sub)))
        print('Number of sources with only pos subtraction: ', len(set_objectId))
        print(f"Number of sources older than {args.max_days} days: {len(set(old))}, specifically {set(old)}")
        '''

    return set_objectId_all


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query kowalski.')
    parser.add_argument('--radius', dest='radius', type=float, required=False, \
    help='Search radius (arcmin)', default = 1.)    
    parser.add_argument('--min-days', dest='min_days', type=float, required=False, \
    help='Minimum time (days) between the first and last alert', default = 3.)
    parser.add_argument('--max-days', dest='max_days', type=float, required=False, \
    help='Maximum time (days) between the first and last alert', default = 14.)
    parser.add_argument('--min-dist', dest='min_dist', type=float, required=False, \
    help='Minimum distance(Mpc) of the CLU galaxies to explore', default = 0.)
    parser.add_argument('--max-dist', dest='max_dist', type=float, required=False, \
    help='Maximum distance(Mpc) of the CLU galaxies to explore', default = 200.)
    parser.add_argument('--min-dec', dest='min_dec', type=float, required=False, \
    help='Minimum declination (celestial, deg) of the CLU galaxies to explore', default = -30.)
    parser.add_argument('--slices', dest='slices', type=int, required=False, \
    help='Number (integer) of slices in which the query will be devided', default = 40)

    args = parser.parse_args()

    from penquins import Kowalski
    import requests

    import numpy as np
    import json
    import pdb

    from astropy.time import Time
    from astropy.io import ascii
    from astropy.io import fits
    from astropy.table import Table

    from astropy import units as u
    from astropy.coordinates import Angle
    from astropy.coordinates import SkyCoord

    #Print a summary of the query input
    print_query_params(args)

    #Read the CLU catalog
    clu=Table.read('CLU_20181213V2.fits')
    clu=clu[clu['dec'] > args.min_dec]
    clu=clu[clu['distmpc'] >= args.min_dist]     
    clu=clu[clu['distmpc'] <= args.max_dist]
    print(f"There are {len(clu)} CLU galaxies in this sample.")

    #Read the secrets
    secrets = ascii.read('secrets.csv', format = 'csv')

    username = secrets['kowalski_user'][0]
    password = secrets['kowalski_pwd'][0]

    #Query kowalski
    sources_kowalski = query_kowalski(username, password, clu, args)

    #Check the CLU science program on the Marshal
    username_marshal = secrets['marshal_user'][0]
    password_marshal= secrets['marshal_pwd'][0]
    
    program_name='Census of the Local Universe'
    clu_sources = get_candidates_growth_marshal(program_name, username_marshal, password_marshal)    

    #For each transient check if it is present in the CLU science program
    check_clu_transients(sources_kowalski, clu_sources)

    print("Done.")



'''
#Plot the data
for galaxy_name, idcoords in zip(galaxy_names_arr, r['result_data']['ZTF_alerts'].keys()):
    all_info=r['result_data']['ZTF_alerts'][idcoords]
    pdb.set_trace()

    jd_arr=[]
    mag_arr=[]
    magerr_arr=[]
    filter_arr=[]
    for info in all_info:
        if info['candidate']['isdiffpos'] != 't':
            continue
        magpsf=info['candidate']['magpsf']
        sigmapsf=info['candidate']['sigmapsf']
        jd=info['candidate']['jd']
        fil=info['candidate']['fid']
        filter_arr.append(fil)
        jd_arr.append(jd)
        mag_arr.append(magpsf)
        magerr_arr.append(sigmapsf)
    if mag_arr!=[]:
        print(info['candidate']['programid'])
        jd0=min(jd_arr)
        jd0_time=Time(jd0, format='jd')
        for i in np.arange(len(jd_arr)):
            jd_arr[i]=jd_arr[i]-jd0
        plt.figure()
        print(galaxy_name, info['objectId'], info['candidate']['ra'], info['candidate']['dec'])
        plt.title(galaxy_name + ' '+info['objectId'])
        #plt.errorbar(jd_arr, mag_arr, yerr=magerr_arr, color='blue', linestyle=' ', marker='o')
        for jj,mm,me,ff in zip(jd_arr, mag_arr, magerr_arr, filter_arr):
            if ff==1:
                fcolor='b'
            if ff==2:
                fcolor='r'
            if ff==3:
                fcolor='y'
            plt.errorbar(jj, mm, yerr=me, color=fcolor, linestyle=' ', marker='o')
        plt.xlabel(f'Days since {jd0_time.iso}')
        plt.gca().invert_yaxis()
        plt.show()

'''
