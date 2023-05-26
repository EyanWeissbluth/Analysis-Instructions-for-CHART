# Analysis-Instructions-for-CHART
Analysis Instructions for CHART
---bash
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
---bash
