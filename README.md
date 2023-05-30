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
            
Step 2: Import data  
We also need to import the data we collected into the notebook.  
Example:  
data_dir = '/home/lmberkhout/CHARTDATA/'  
paths = ['PHX_BABY1_11_11_2022','PHX_BABY2_11_11_2022'] 

    data_dir = '(directory_name)'
    paths = ['file names'] 

Step 3: Making a simple plot with our data  
To visulise the data we have collected, we will make a plot. Each plot will corrospond to a file that was imported from the previous step. 

    ntrials = len(paths)
    print(ntrials)
    data = []
    mdata = []
    spectra = [[] for _ in range(ntrials)]
    freqs = [[] for _ in range(ntrials)]

    for i in range(ntrials):
        d, m = chart.analysis.read_run(directory=data_dir + paths[i])
        d = np.array(d)
        data.append(d)
        mdata.append(m)
        nchans = m[0]['vector_length']
        nremove = nchans // 32
    

    for j in range(ntrials):
        for d, m in zip(data[j], mdata[j]):
            spectrum = np.mean(d, axis=0) 
            spectrum = spectrum[nremove:-nremove]
            spectrum = 10 * np.log10(spectrum)
            frequencies = ((np.arange(m['vector_length']) - m['vector_length'] / 2)
                               * m['samp_rate'] / m['vector_length'] + m['frequency'])
            frequencies = 1e-9 * frequencies[nremove:-nremove]
            spectra[j].append(spectrum)
            freqs[j].append(frequencies)

    for k in range(len(spectra[j]) - 1):
        spec1 = spectra[j][k]
        spec2 = spectra[j][k + 1]
        freq1 = freqs[j][k]
        freq2 = freqs[j][k + 1]
        ncommon = np.sum([1 if f in freq2 else 0 for f in freq1])
        spec2 += np.median(spec1[-ncommon:]) - np.median(spec2[:ncommon])
        spectra[j][k + 1] = spec2
    
    plt.figure(figsize=(10,10))

    for j in range(ntrials):
        if j == 0:
            ax = plt.subplot(ntrials, 1, j + 1)
            ax.set_title('Control Scan')
        else:
            plt.subplot(ntrials, 1, j + 1, sharex=ax, sharey=ax)
            plt.title('Data Scan')
        for f, s in zip(freqs[j], spectra[j]):
            plt.plot(f, s)
        plt.ylabel('[dB]')
        plot_f0()
    
    plt.ylim(30,45)
    plt.xlabel('Frequency [GHz]')

Step 4: Bandpass Calibration   
Now we need to apply the Bandpass filter to our scan so we can isolate the wavelengths we want to observe.   

    ntrials = len(paths)
    data = []
    mdata = []
    bps= []
    spectra = [[] for _ in range(ntrials)]
    freqs = [[] for _ in range(ntrials)]

    for i in range(ntrials):
        d, m = chart.analysis.read_run(directory=data_dir + paths[i])
        d = np.array(d)
        data.append(d)
        mdata.append(m)
        nchans = m[0]['vector_length']
        nremove = nchans // 32
    for i in range(ntrials):
        d, m = chart.analysis.read_run(directory=data_dir + paths[i])
        d = np.array(d)
        data.append(d)
        mdata.append(m)
        # Rough estimate for bandpass
        nchans = m[0]['vector_length']
        levels = np.median(d[:, :, nchans // 4:(-nchans // 4)], axis=(1, 2))
        rescaled = d / levels.reshape(-1, 1, 1)
        bp = np.median(rescaled, axis=(0, 1))
        bps.append(bp)
    
    for j in range(ntrials):
        for d, m in zip(data[j], mdata[j]):
            spectrum = np.mean(d, axis=0) / bps[j]
            spectrum = spectrum[nremove:-nremove]
            spectrum = 10 * np.log10(spectrum)
            frequencies = ((np.arange(m['vector_length']) - m['vector_length'] / 2)
                               * m['samp_rate'] / m['vector_length'] + m['frequency'])
            frequencies = 1e-9 * frequencies[nremove:-nremove]
            spectra[j].append(spectrum)
            freqs[j].append(frequencies)

        for k in range(len(spectra[j]) - 1):
            spec1 = spectra[j][k]
            spec2 = spectra[j][k + 1]
            freq1 = freqs[j][k]
            freq2 = freqs[j][k + 1]
            ncommon = np.sum([1 if f in freq2 else 0 for f in freq1])
            spec2 += np.median(spec1[-ncommon:]) - np.median(spec2[:ncommon])
            spectra[j][k + 1] = spec2  
            
Step 5: Flagging out Radio Frequencey Interference  
