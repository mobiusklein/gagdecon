import re

import numpy as np

from ms_deisotope import deconvolution, scoring
from ms_peak_picker import PeakIndex, PeakSet, pick_peaks


class PeakListReader(object):
    def __init__(self, path):
        self.path = path
        self.peaklist = []
        self.parse()

    def parse(self):
        path = self.path
        for line in open(path):
            if not line:
                break
            if re.match(r"\d+\.?[0-9]*", line):
                mz, intensity = map(float, re.findall(r"\d+\.?[0-9]*", line))

                self.peaklist.append(deconvolution.FittedPeak(mz, intensity, 0, 0, 0, 0, 0))
        self.peaklist = PeakIndex(np.array([]), np.array([]), PeakSet(self.peaklist))


class PointListReader(object):
    def __init__(self, path):
        self.path = path
        self.mzs = []
        self.intensities = []
        self.parse()

    def parse(self):
        path = self.path
        for line in open(path):
            if not line:
                break
            if re.match(r"\d+\.?[0-9]*", line):
                mz, intensity = map(float, re.findall(r"\d+\.?[0-9]*", line))
                self.mzs.append(mz)
                self.intensities.append(intensity)
        self.mzs = np.array(self.mzs)
        self.intensities = np.array(self.intensities)


def pick_peaks_from_file(path):
    reader = PointListReader(path)
    return pick_peaks(reader.mzs, reader.intensities, "lorenztian")


def read_peaklist(path):
    peakreader = PeakListReader(path)
    return peakreader.peaklist


def deconvolute(peaklist, compositions, charge_range=(-1, -10), mass_error_tolerance=2e-5):
    scorer = scoring.ScaledGTestFitter()
    scorer.select.minimum_score = 1.5
    dec = deconvolution.CompositionListDeconvoluter(peaklist, compositions, scorer=scorer)
    return DeconvolutedResult.wrap(dec.deconvolute(charge_range=charge_range, error_tolerance=mass_error_tolerance))


class DeconvolutedResult(object):
    @classmethod
    def wrap(cls, output):
        return [cls(*args) for args in output]

    def __init__(self, composition, peak, fit):
        self.composition = composition
        self.peak = peak
        self.fit = fit

    @property
    def ppm_error(self):
        acc = []
        for obs, exp in zip(self.fit.experimental, self.fit.theoretical):
            acc.append((obs.mz - exp.mz) / exp.mz)
        return sum(acc) / float(len(acc))

    def __getitem__(self, key):
        try:
            if key == "mz":
                return self.peak.mz
            elif key == "charge":
                return self.peak.charge
            elif key == "score":
                return self.peak.score
            elif key == "chain_length":
                return self.composition.chain_length
            elif key == "intensity":
                return self.peak.intensity
            elif key == "ppm_error":
                return self.ppm_error
            else:
                return self.composition.glycan_composition[key]
        except:
            return self.composition.glycan_composition[key]

    def __iter__(self):
        yield self.composition
        yield self.peak
        yield self.fit
