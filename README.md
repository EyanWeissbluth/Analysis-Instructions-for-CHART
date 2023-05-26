# Analysis-Instructions-for-CHART
To take the raw data we have collected and turn it into a rotation curve, we will have to run the data through some code. 

Step 1: Import package  
We need to import a few packages so that python can perform certain mathematical and visual functions. 



    %matplotlib inline
    import numpy as np
    import matplotlib.pyplot as plt
    import chart
    from scipy import signal
    from scipy.ndimage import median_filter
    %matplotlib inline
    import numpy as np
    import matplotlib.pyplot as plt
    import chart
    from astropy import units as u
    from astropy.coordinates import SpectralCoord, EarthLocation, SkyCoord
    from astropy.time import Time
    import pandas as pd
    import scipy.stats as stays
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import math

    f0 = 1.420405751768  # GHz
    speed_of_light = 299792458  # m/s
    def plot_f0(lims=[20, 50], xval=f0):

    plt.plot([xval, xval], lims, '--k', lw=0.5)
    def f2v(freq):
            return -(np.array(freq)-f0) * speed_of_light / np.array(freq)
            
