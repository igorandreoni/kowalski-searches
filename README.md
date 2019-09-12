# kowalski-searches
Search for transients using the kowalski database.


## Required software
Python 3.6 or later version

###Requited packages
numpy <br>
astropy <br>
matplotlib <br>
healpy <br>
ligo.skymap (0.1.11, or any version that properly supports multi-resolution skymaps)

## Query skymap 

usage: query_skymap.py [-h] [--skymap SKYMAP_FILENAME] [--level LEVEL]
                       [--fov FOV] [--ra-center RA_CENTER [RA_CENTER ...]]
                       [--dec-center DEC_CENTER [DEC_CENTER ...]]
                       [--radius RADIUS] [--after-trigger AFTER_TRIGGER]
                       [--jd-trigger JD_TRIGGER] [--min-days MIN_DAYS]
                       [--max-days MAX_DAYS] [--ndethist NDETHIST_MIN]
                       [--slices SLICES] [--out OUT]

Query kowalski with cone searches within a LIGO/Virgo skymap contour.

optional arguments:
  -h, --help            show this help message and exit
  --skymap SKYMAP_FILENAME
                        Skymap filename
  --level LEVEL         Enclosed probability
  --fov FOV             Field of view of each cone (radius, in arcmin)
  --ra-center RA_CENTER [RA_CENTER ...]
                        Right ascension of the center (array, in degrees)
  --dec-center DEC_CENTER [DEC_CENTER ...]
                        Declination of the center (array, in degrees)
  --radius RADIUS       Search radius (min), by default radius = fov
  --after-trigger AFTER_TRIGGER
                        Query only alerts whose first detection occurred after
                        a certain date. If this boolean value is True, then
                        --jd-trigger must be given.
  --jd-trigger JD_TRIGGER
                        Julian Day of the trigger
  --min-days MIN_DAYS   Minimum time (days) between the first and last alert
  --max-days MAX_DAYS   Maximum time (days) between the first and last alert
  --ndethist NDETHIST_MIN
                        Minimum number of detections
  --slices SLICES       Number (integer) of slices in which the query will be
                        devided
  --out OUT             Output filename

example: python query_skymap.py --skymap LALInference.fits.gz --level 90 --slice 100 --jd-trigger 2458736.5599421295 --min-days 0.02 --max-days 30.


## Check skymap
WIP

## Query CLU
usage: query_clu.py [-h] [--radius RADIUS] [--min-days MIN_DAYS]
                    [--max-days MAX_DAYS] [--min-dist MIN_DIST]
                    [--max-dist MAX_DIST] [--min-dec MIN_DEC]
                    [--slices SLICES]

Query kowalski with cone searces centred on CLU galaxies.

optional arguments:
  -h, --help           show this help message and exit
  --radius RADIUS      Search radius (arcmin)
  --min-days MIN_DAYS  Minimum time (days) between the first and last alert
  --max-days MAX_DAYS  Maximum time (days) between the first and last alert
  --min-dist MIN_DIST  Minimum distance(Mpc) of the CLU galaxies to explore
  --max-dist MAX_DIST  Maximum distance(Mpc) of the CLU galaxies to explore
  --min-dec MIN_DEC    Minimum declination (celestial, deg) of the CLU
                       galaxies to explore
  --slices SLICES      Number (integer) of slices in which the query will be
                       devided


