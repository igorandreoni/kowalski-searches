__author__ = "Igor Andreoni"
__license__ = "MIT"
__email__ = "andreoni@caltech.edu"

import numpy as np
import matplotlib.pyplot as plt


def get_dist(m, M):
    """
    Function that calculates the distance
    of an object with absolute magnitude M
    and apparent magnitude m

    ---
    Parameters

    m float
        apparent mag
    M float
        absolute mag

    ---
    Returns

    distance float
        distance in pc
    """
    D = 10. # 10pc
    distance = D * 10**((m-M)/5.)

    return distance


def get_exptime(n, t_tot, overhead=10):
    """
    Function that calculates the exposure time
    per image

    ---
    Parameters

    n int or array
        number of fields
    t_tot float
        total ovservation time (s)
    overhead float
        overhead time per exposure (s)

    ---
    Returns

    exptime float
        exposure time per image
    """

    exptime = (t_tot - n*overhead)/n

    return exptime


def get_maglim(maglim_known, exptime_known, exptime):
    """
    Function that calculates the mag limit for
    a given exposure time, starting from a known limit

    ---
    Parameters

    maglim_known float or array
        known magnitude limit for exptime_known
    exptime_known float or array
        known exposure time
    exptime float or array
        new expoure time

    ---
    Returns

    maglim float
        mag limit rescaled for new exposure time
    """

    maglim = maglim_known + 2.5*np.log10(exptime/exptime_known)

    return maglim


def get_volume(n, dist, area_field=47.):
    """
    Function that calculates the volume based
    on the number of fields and the total
    observation time available

    ---
    Parameters

    n int or array
        number of fields
    dist float or array
        distance to the source
    area_field float
        area per field (deg2)

    ---
    Returns

    vol float or array
        volume

    """

    # Calculate the volume (sphere)
    vol_sphere = 4./3. * np.pi * dist**3 

    # ..just for the given area
    vol = vol_sphere * n * area_field / (4*np.pi*(180/np.pi)**2)

    return vol


if __name__ == '__main__':
    """
    Snippet to plot the explored volume varying the number
    of fields in a given time frame.
    """

    # Say that the total observation time is fixed (s)
    t_tot = 10.*60*60

    # Min exposure time
    exptime_min = 30

    # Min magnitude limit
    maglim_min = 20.5

    # Overhead between exposures
    overhead = 10.

    # Absolute mag of your favorite object
    M = -16

    # Area of your individual field (deg2)
    area_field = 47.

    # Let's vary the number of fields in a range
    n_fields = np.arange(1, int(t_tot/(exptime_min + overhead)))

    # For each given number of fields get the exposure
    # time per field
    exptimes = get_exptime(n_fields, t_tot, overhead=overhead)
    
    # Get the magnitude limits
    maglim = get_maglim(maglim_min, exptime_min, exptimes)

    # Get the distances
    distance = get_dist(maglim, M)

    # Get volumes
    volume = get_volume(n_fields, distance, area_field=area_field)

    print("n_fields, exptime, maglim, distance, volume (pc^3)")
    for n, e, mm, dd, vv in zip(n_fields, exptimes, maglim, distance, volume):
        print(n, e, mm, dd, vv) 

    # Plot the results
    fig, ax = plt.subplots()
    ax.plot(n_fields, volume/10**18, label=f"M={M}")
    ax.set_title(f"total time {t_tot/3600.} hr, {area_field:.0f} deg2")
    ax.set_xlabel("Number of fields")
    ax.set_ylabel("Volume (Mpc$^3$)")
    ax.set_yscale('log')

    plt.legend()
    plt.show()
