# kowalski-searches
Search for transients using the kowalski database.


## Required software
Python 3.6 or later version

### Required packages
numpy <br>
astropy <br>
matplotlib <br>
healpy <br>
ligo.skymap (0.1.11, or any version that properly supports multi-resolution skymaps)

## Query skymap 

usage: query_skymap.py [-h] [--skymap SKYMAP_FILENAME] [--level LEVEL] <br>
                       [--fov FOV] [--ra-center RA_CENTER [RA_CENTER ...]] <br>
                       [--dec-center DEC_CENTER [DEC_CENTER ...]] <br>
                       [--radius RADIUS] [--after-trigger AFTER_TRIGGER] <br>
                       [--jd-trigger JD_TRIGGER] [--min-days MIN_DAYS] <br>
                       [--max-days MAX_DAYS] [--ndethist NDETHIST_MIN] <br>
                       [--slices SLICES] [--out OUT] <br>
 <br>
Query kowalski with cone searches within a LIGO/Virgo skymap contour. <br>
 <br>
optional arguments: <br>
  -h, --help            show this help message and exit <br>
  --skymap SKYMAP_FILENAME 
                        Skymap filename <br>
  --level LEVEL         Enclosed probability <br>
  --fov FOV             Field of view of each cone (radius, in arcmin) <br>
  --ra-center RA_CENTER [RA_CENTER ...]
                        Right ascension of the center (array, in degrees) <br>
  --dec-center DEC_CENTER [DEC_CENTER ...]
                        Declination of the center (array, in degrees) <br>
  --radius RADIUS       Search radius (min), by default radius = fov <br>
  --after-trigger AFTER_TRIGGER 
                        Query only alerts whose first detection occurred after 
                        a certain date. If this boolean value is True, then
                        --jd-trigger must be given (default=True).  <br>
  --jd-trigger JD_TRIGGER 
                        Julian Day of the trigger. If not given, it will be <br>
                        read from the skymap header. <br>
  --min-days MIN_DAYS   Minimum time (days) between the first and last alert <br>
  --max-days MAX_DAYS   Maximum time (days) between the first and last alert <br>
  --within-days WITHIN_DAYS
                        Maximum time (days) between the jd-trigger and the  <br>
                        first alert  <br>
  --ndethist NDETHIST_MIN
                        Minimum number of detections <br>
  --slices SLICES       Number (integer) of slices in which the query will be
                        devided <br>
  --out OUT             Output filename <br>
  --ingest INGEST       Use with caution! Ingest the candidates to the <br>
                        scanning page on the GROWTH marshal <br>
  --ingest-program INGEST_PROGRAM
                        Ingest the candidates to the scanning page of the <br>
                        given GROWTH marshal program <br>
  --phi PHI             Phi angle rotation of the skymap <br>
  --theta THETA         Theta angle rotation of the skymap <br>

example: python query_skymap.py --skymap LALInference.fits.gz --level 90 --slice 20 --min-days 0.02 --max-days 30. --within-days 2.0 <br>


## Check skymap
WIP <br>

## Query CLU
usage: query_clu.py [-h] [--radius RADIUS] [--min-days MIN_DAYS] <br>
                    [--max-days MAX_DAYS] [--min-dist MIN_DIST] <br>
                    [--max-dist MAX_DIST] [--min-dec MIN_DEC] <br>
                    [--slices SLICES]
 <br>
Query kowalski with cone searces centred on CLU galaxies. <br>
 <br>
optional arguments: <br>
  -h, --help           show this help message and exit <br>
  --radius RADIUS      Search radius (arcmin) <br>
  --min-days MIN_DAYS  Minimum time (days) between the first and last alert <br>
  --max-days MAX_DAYS  Maximum time (days) between the first and last alert <br>
  --min-dist MIN_DIST  Minimum distance(Mpc) of the CLU galaxies to explore <br>
  --max-dist MAX_DIST  Maximum distance(Mpc) of the CLU galaxies to explore <br>
  --min-dec MIN_DEC    Minimum declination (celestial, deg) of the CLU <br>
                       galaxies to explore <br>
  --slices SLICES      Number (integer) of slices in which the query will be <br>
                       devided <br>


