#! /usr/bin/env python
##########################################################################
# Nsap - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import json
import nibabel
import numpy
import matplotlib

matplotlib.use("AGG")

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def snr_percent_fluctuation_and_drift(image_file, repetition_time, roi_size,
                                      output_directory):
    """ Compute the fluctuation and drift on a central roi.

    <process>
        <return name="snap_fluctuation_drift" type="File" desc="A functional
        volume."/>
        <return name="fluctuation_drift_file" type="File" desc="A functional
        volume."/>
        <input name="image_file" type="File" desc="A functional volume."/>
        <input name="repetition_time" type="Float" desc="The fMRI sequence
            repetition time (in seconds)."/>
        <input name="roi_size" type="Int" desc="The central ROI size used to
            compute the fluuctation and drift."/>
        <input name="output_directory" type="Directory" desc="The destination
            folder."/>
    </process>
    """
    # Load the functional volume array
    array = load_fmri_dataset(image_file)

    # Compute the drift and fluctuation
    (average_intensity, polynomial, residuals, fluctuation,
     drift) = get_snr_percent_fluctuation_and_drift(array, roi_size=roi_size)
    spectrum = get_residuals_spectrum(residuals, repetition_time * 10e-3)

    # Compute the snr
    signal_array = get_fmri_signal(array)
    ssn_array = get_static_spatial_noise(array)
    snr = get_spatial_noise_ratio(signal_array, ssn_array, array.shape[3],
                                  roi_size=roi_size)

    # Compute a weisskoff analysis on the fluctuation
    fluctuations, theoretical_fluctuations = get_weisskoff_analysis(
        array, max_roi_size=roi_size)

    # Save the result in a json
    fluctuation_drift_file = os.path.join(
        output_directory, "snr_fluctuation_drift.json")
    results = {
        "drift": drift * 100,
        "snr": 10 * numpy.log10(snr)
    }
    with open(fluctuation_drift_file, "w") as json_data:
        json.dump(results, json_data)

    # Display result in a pdf
    snap_fluctuation_drift = os.path.join(
        output_directory, "snr_fluctuation_drift.pdf")
    pdf = PdfPages(snap_fluctuation_drift)
    try:
        # Plot one figure per page
        fig = time_series_figure(average_intensity, polynomial, drift, snr)
        pdf.savefig(fig)
        plt.close()
        fig = spectrum_figure(spectrum)
        pdf.savefig(fig)
        plt.close()
        fig = weisskoff_figure(fluctuations, theoretical_fluctuations)
        pdf.savefig(fig)
        plt.close()

        # Close the pdf
        pdf.close()
    except:
        pdf.close()
        raise

    return snap_fluctuation_drift, fluctuation_drift_file


def load_fmri_dataset(fmri_file):
    """ Load a functional volume

    Parameters
    ----------
    fmri_file: str (mandatory)
        the path to the functional volume.

    Returns
    -------
    array: array [X,Y,Z,N]
        the functional volume.
    """
    img = nibabel.load(fmri_file)
    return img.get_data()


def get_fmri_signal(array):
    """ The fmri signal is the the average voxel intensity across time

    A signal summary value is obtained from an ROI placed in the center.

    Parameters
    ----------
    array: array [X,Y,Z,N]
        the functional volume.

    Returns
    -------
    signal: array [X,Y,Z]
        the fmri signal.

    Notes
    -----
    FMRI Data Quality, Pablo Velasco, 2014
    SINAPSE fMRI Quality Assurance, Katherine Lymer
    """
    return numpy.average(array, 3)


def get_fmri_temporal_fluctuation_noise(array):
    """ The fluctuation noise image is produced by substracting for each voxel
    a trend line estimated from the data (the BIRN function uses a second
    order polynomial - It's effectively a standard deviation image).

    By removing the trend from the data, the fluctuation noise image is an
    image of the standard deviation of the residuals, voxel by voxel.

    A linear trend may indicate a systematic increase or decrease in the data
    (caused by eg sensor drift)

    Parameters
    ----------
    array: array [X,Y,Z,N]
        the functional volume.

    Returns
    -------
    tfn: array [X,Y,Z]
        the temporal fluctuation noise array.
    """
    # Flatten input array
    flat_array = array.ravel()
    voxels_per_volume = reduce(lambda x, y: x * y, array.shape[: -1], 1)

    # Compute the temporal fluctuation noise
    tfn = numpy.ndarray((voxels_per_volume,), dtype=numpy.single)
    x = numpy.arange(array.shape[3])
    for i in range(voxels_per_volume):

        # Get the temporal signal at one voxel
        y = flat_array[i::voxels_per_volume]

        # Estimate a second order polynomial trend
        polynomial = numpy.polyfit(x, y, 2)
        model = numpy.polyval(polynomial, x)

        # Compute the residuals
        residuals = y - model

        # Compute the emporal fluctuation noise
        tfn[i] = numpy.std(residuals)

    return tfn.reshape(array.shape[:-1])


def get_signal_to_fluctuation_noise_ratio(signal_array, tfn_array,
                                          roi_size=10):
    """ The SFNR image is is obtained by dividing, voxel by voxel,
    the mean fMRI signal image by the temporal fluctuation image.

    A 21 x 21 voxel ROI, placed in the center of the image, is created.
    The average SFNR across these voxels is the SFNR summary value.

    Parameters
    ----------
    signal_array: array [X,Y,Z]
        the fmri signal.
    tfn_array: array [X,Y,Z]
        the temporal fluctuation noise array.
    roi_size: int (default 10)
        the size of the central roi used to get the summary indice.

    Returns
    -------
    sfnr_array: array [X,Y,Z]
        the signal to fluctuation noise ratio array.
    sfnr_summary: float
        the signal to fluctuation noise ratio average on a central roi.
    """
    # Compute the signal to fluctuation noise ratio
    sfnr_array = signal_array / (tfn_array + numpy.finfo(float).eps)

    # Create a central roi
    center = numpy.round(sfnr_array.shape / 2)
    roi = sfnr_array[center[0] - roi_size: center[0] + roi_size,
                     center[1] - roi_size: center[1] + roi_size,
                     center[2]]

    # Compute the signal to fluctuation noise ratio summary
    sfnr_summary = numpy.average(roi)

    return sfnr_array, sfnr_summary


def get_static_spatial_noise(array):
    """ In order to measure the spatial noise, the sum of the odd and even
    numbered volumes are calculated separately.
    The static patial noise is then approximated by the differences of
    these two sums.

    If there is no drift in either amplitude or geometry across time series,
    then the difference image will show no structure.

    Parameters
    ----------
    array: array [X,Y,Z,N]
        the functional volume.

    Returns
    -------
    ssn_array: array [X,Y,Z]
        the static spatial noise.   
    """
    shape_t = array.shape[3]
    odd_array = array[..., range(1, shape_t, 2)]
    odd_sum_array = numpy.sum(odd_array, 3)
    even_array = array[..., range(0, shape_t, 2)]
    even_sum_array = numpy.sum(even_array, 3)
    ssn_array = odd_sum_array - even_sum_array
    return ssn_array


def get_spatial_noise_ratio(signal_array, ssn_array, nb_time_points,
                            roi_size=10):
    """ A central ROI is placed in the center of the static spatial noise
    image. The SNR is the signal summary value divided by by the square root of the
    variance summary value divided by the numbre of time points

    SNR = (signal summary value)/sqrt((variance summary value)/#time points).

    Parameters
    ----------
    signal_array: array [X,Y,Z]
        the fmri signal.
    ssn_array: array [X,Y,Z]
        the static spatial noise.  
    nb_time_points: int
        the number of time points in the fMRI serie.
    roi_size: int (default 10)
        the size of the central roi used to get the summary indice.

    Returns
    -------
    snr: float
        the fMRI image signal to noise ratio.
    """
    # Create a central roi
    center = numpy.round(numpy.asarray(signal_array.shape) / 2)
    signal_roi_array = signal_array[
        center[0] - roi_size: center[0] + roi_size,
        center[1] - roi_size: center[1] + roi_size,
        center[2]]
    ssn_roi_array = ssn_array[
        center[0] - roi_size: center[0] + roi_size,
        center[1] - roi_size: center[1] + roi_size,
        center[2]]

    # Get the signal and variance summaries
    signal_summary = numpy.average(signal_roi_array)
    variance_summary = numpy.var(ssn_roi_array)

    # Compute the SNR
    snr = signal_summary / numpy.sqrt(variance_summary / nb_time_points)

    return snr


def get_snr_percent_fluctuation_and_drift(array, roi_size=10):
    """ A time-series of the average intensity within a 21 x 21 voxel ROI
    centered in the image is calculated.

    A second-order polynomial trend is fit to the volume number vs
    average intensity.

    The mean signal intensity of the time-series (prior to detrending) and SD
    of the residuals after subtracting the fit line from the data
    are calculated.

    fluctuation = 100*(SD of the residuals)/(mean signal intensity)
    drift = ((max fit value) - (min fit value)) /
             (mean signal intensity signal)

    Parameters
    ----------
    array: array [X,Y,Z,N]
        the functional volume.
    roi_size: int (default 10)
        the size of the central roi used to get the summary indice.

    Returns
    -------
    average_intensity: array [N]
        the average voxel intensity across time.
    polynomial: array [N]
        the second-order polynomial used to remove the slow drift.
    residuals: array [N]
        the signal/model residuals.
    fluctuation: float
        the fluctuation value computed on a central roi.
    drift: float
        the signal temporal drift on a central roi.
    """
    # Create a central roi
    center = numpy.round(numpy.asarray(array.shape) / 2)
    roi = array[center[0] - roi_size: center[0] + roi_size,
                center[1] - roi_size: center[1] + roi_size,
                center[2]]
    shape = roi.shape

    # Compute the mean signal intensity
    mean_signal_intensity = numpy.average(roi)

    # Compute the average voxel intensity across time
    average_intensity = numpy.sum(numpy.sum(roi, 0), 0)
    average_intensity /= shape[0] * shape[1]

    # Compute the temporal fluctuation noise
    # > a second-order polynomial detrending to remove the slow drift
    x = numpy.arange(array.shape[3])
    polynomial = numpy.polyfit(x, average_intensity, 2)
    average_intensity_model = numpy.polyval(polynomial, x)
    # > a fluctuation value is calculated as the temporal standard deviation
    # of the residual variance of each voxel after the detrending 
    residuals = average_intensity - average_intensity_model
    fluctuation = numpy.std(residuals) / mean_signal_intensity

    # Compute the drift
    fluctuation = numpy.std(residuals) / mean_signal_intensity
    drift = (average_intensity_model.max() -
            average_intensity_model.min()) / mean_signal_intensity

    return average_intensity, polynomial, residuals, fluctuation, drift


def get_residuals_spectrum(residuals, repetition_time):
    """ Residuals of the mean signal intensity fit are submitted to
    a fast Fourier transform (FFT).

    Parameters
    ----------
    residuals: array [N]
        the residuals between the time serie values and the model.
    repetition_time: float
        the repetition time in ms.
    """
    fft = numpy.fft.rfft(residuals)
    if residuals.size % 2 == 0:
        # Discard real term for frequency n/2
        fft = fft[:-1]
    fftfreq = numpy.fft.fftfreq(len(residuals), repetition_time)
    return numpy.vstack((fftfreq[:len(fft)], numpy.abs(fft)))


def get_weisskoff_analysis(array, max_roi_size=30):
    """ The Weisskoff analysis provides another measure of scanner
    stability included in the GSQAP.

    Parameters
    ----------
    array: array [X,Y,Z,N]
        the functional volume.
    max_roi_size: int
        the max size of the central size that will be used to compute the
        drift and fluctuation.

    Returns
    -------
    fluctuation: 2-uplet
        the roi sizes array and the fluctuation array.
    theoretical_fluctuation: 2-uplet
        the roi sizes array and the theoretical fluctuation array.
    """
    # Generate all the roi sizes we want to consider
    roi_sizes = numpy.arange(1, max_roi_size + 1)

    # Compute the fluctuation for each roi size
    fluctuation = [get_snr_percent_fluctuation_and_drift(array, s)[-2]
                   for s in roi_sizes]

    # The theorical fluctuation is given by the ine voxel roi fluctuation
    theoretical_fluctuation = fluctuation[0] / roi_sizes

    return (numpy.vstack((roi_sizes, fluctuation)),
            numpy.vstack((roi_sizes, theoretical_fluctuation)))


def weisskoff_figure(fluctuations, theoretical_fluctuations):
    """ Return a matplotlib figure containing the Weisskoff analysis.

    Parameters
    ----------
    fluctuation: 2-uplet
        the roi sizes array and the fluctuation array.
    theoretical_fluctuation: 2-uplet
        the roi sizes array and the theoretical fluctuation array.
    """
    figure = plt.figure()
    plot = figure.add_subplot(111)
    plot.grid(True, which="both", ls="-")
    plt.title("Weisskoff analysis")
    plot.plot(fluctuations[0], 100 * fluctuations[1], "ko-", fillstyle="full")
    plot.plot(theoretical_fluctuations[0], 100 * theoretical_fluctuations[1],
              "ko-", markerfacecolor="w")
    plot.axes.loglog()
    plot.axes.set_xlabel("ROI width (pixels)")
    plot.axes.set_ylabel("Fluctuation (%)")
    plot.xaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
    plot.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
    plot.legend(("Measured", "Theoretical"), "upper right")
    return figure


def spectrum_figure(spectrum):
    """ Return a matplotlib figure containing the Fourier spectrum, without its
    DC coefficient.
    """
    figure = plt.figure()
    plot = figure.add_subplot(111)
    plot.grid()
    plt.title("Residual Spectrum")
    plot.plot(spectrum[0, 1:], spectrum[1, 1:], "k-")
    plot.axes.set_xlabel("Frequency (Hz)")
    plot.axes.set_ylabel("Magnitude")
    return figure


def time_series_figure(time_series, polynomial, drift, snr):
    """ Return a matplotlib figure containing the time series and its
    polynomial model.
    """
    figure = plt.figure()
    plot = figure.add_subplot(111)
    plot.grid()
    plt.title("Drift: {0: .1f}% - SNR: {1: .1f}dB".format(
        drift * 100, 10 * numpy.log10(snr)))
    x = numpy.arange(2, 2 + len(time_series))
    model = numpy.polyval(polynomial, x)
    plot.plot(x, time_series, "k-")
    plot.plot(x, model, "k-")
    plot.axes.set_xlabel("Volume number")
    plot.axes.set_ylabel("Intensity")
    return figure

