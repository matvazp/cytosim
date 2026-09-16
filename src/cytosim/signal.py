"""From pulse shape to ADC sampling-rate requirement.

Two complementary criteria are offered; the final number is their maximum.

1. Time-domain: N samples across the pulse FWHM
       f_s = samples_per_fwhm / FWHM
2. Frequency-domain: oversample the Nyquist rate of the pulse bandwidth
       f_s = oversampling * 2 * f_3dB
   For a Gaussian pulse of width sigma_t, |S(f)|^2 = exp(-4 pi^2 sigma_t^2 f^2),
   so f_3dB = sqrt(ln 2) / (2 pi sigma_t).
"""

import math

import numpy as np


def gaussian_bandwidth_3db(sigma_t: float) -> float:
    """-3 dB (power) bandwidth [Hz] of a Gaussian pulse."""
    return math.sqrt(math.log(2.0)) / (2.0 * math.pi * sigma_t)


def pulse_sigma_t(t: np.ndarray, s: np.ndarray) -> float:
    """RMS width [s] of a sampled pulse (second moment)."""
    s = np.clip(s, 0.0, None)
    mean = np.trapezoid(t * s, t) / np.trapezoid(s, t)
    var = np.trapezoid((t - mean) ** 2 * s, t) / np.trapezoid(s, t)
    return math.sqrt(var)


def sample_rate_from_fwhm(fwhm: float, samples_per_fwhm: float) -> float:
    return samples_per_fwhm / fwhm


def sample_rate_from_bandwidth(f_3db: float, oversampling: float) -> float:
    return oversampling * 2.0 * f_3db
