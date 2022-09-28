'''
Created on Oct. 14, 2021

@author: Ryan
'''
from __future__ import annotations
import pyms
import scipy.stats
import numpy
import logging
import sqlalchemy.orm

from . import sql_interface
from pandas._libs import index

DEBUG_STATS = False


class MassSpectrumEntry(sql_interface.Base):
    __tablename__ = "ion_peak"

    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    mass = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    intensity = sqlalchemy.Column(sqlalchemy.Float, nullable=False)
    area = sqlalchemy.Column(sqlalchemy.Float, nullable=False)
    spectrum_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey('mass_spectrum.id'))
    spectrum = sqlalchemy.orm.relationship("MassSpectrum_SQL", back_populates="values")


class MassSpectrum_SQL(sql_interface.Base):
    __tablename__ = "mass_spectrum"

    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    values = sqlalchemy.orm.relationship("MassSpectrumEntry", back_populates="spectrum")


class MassSpectrum(pyms.Spectrum.MassSpectrum):
    ''' '''
    def __init__(self, mass_spec):
        '''
        Convert a pyms.MassSpectrum to a database
        '''
        if mass_spec is None:
            super().__init__([], [])
        else:
            super().__init__(mass_spec.mass_list, mass_spec.intensity_list)
        self.sql_ref = MassSpectrum_SQL()

    @classmethod
    def merge_mass_spectrum_list(cls, mass_spectrum_list: list, required_correlation: float=0.50) -> ('MassSpectrum', numpy.array):
        ''' Checks the relationship of a list of mass spectra, creating a mapping of which are most similar to the average and provides a summary average ms object '''
        mass_lists = [x.mass_list for x in mass_spectrum_list]
        if not (numpy.array(mass_lists) == mass_lists[0]).all():
            raise ValueError(f"Mass lists are not equal")
        spectrum_lists = numpy.array([x.intensity_list for x in mass_spectrum_list])
        input_size = len(mass_spectrum_list)
        exclude_mapping = numpy.full(input_size, True)
        for _pass_count in range(input_size):
            assert(sum(exclude_mapping) >= 1)
            worst_index = -1
            worst_corr = numpy.Inf
            averaged_spectrum = numpy.mean(spectrum_lists[exclude_mapping], axis=0)
            for index, spectrum in enumerate(spectrum_lists[exclude_mapping]):
                corr = cls.single_compare(averaged_spectrum, spectrum)
                if corr < required_correlation and corr < worst_corr:
                    worst_index = index
                    worst_corr = corr

            if worst_index != -1:
                exclude_mapping[worst_index] = False
            else:
                break

        pyms_ms = pyms.Spectrum.MassSpectrum(mass_lists[0], averaged_spectrum)
        return MassSpectrum(pyms_ms), exclude_mapping

    @classmethod
    def compare_to_target_ms(cls, mass_spectrum_target: MassSpectrum, comparison_ms_list: list, verify_highest=False):
        '''
        Not symmetric, mass_spectrum_target is used to define the valid mass range for comparison
        '''
        # TODO: on verify highest, check that no ions in ms_list are substantially higher than lowest target ion
        target_range_low = mass_spectrum_target.mass_list[0]
        target_range_high = mass_spectrum_target.mass_list[-1]
        target_ms_dict = dict(zip(mass_spectrum_target.mass_list, mass_spectrum_target.intensity_list))
        result = []
        for mass_spectrum_compare in comparison_ms_list:
            try:
                compare_ms_dict = dict(zip(mass_spectrum_compare.mass_list, mass_spectrum_compare.intensity_list))
                compare_masses = [x for x in compare_ms_dict.keys() if target_range_low < x < target_range_high]
                combined_keys = set(target_ms_dict.keys()).union(compare_masses)
                target_intensities = cls._dict_align(combined_keys, target_ms_dict)
                compare_intensities = cls._dict_align(combined_keys, compare_ms_dict)
                corr = cls.single_compare(target_intensities, compare_intensities)
                result.append(corr)
            except Exception as _e:
                logging.warning("Failed to compare ms")
                result.append(0.0)
        return result

    @staticmethod
    def convert_from_dict(input_dict) -> 'MassSpectrum':
        result = MassSpectrum(None)
        masses = list(input_dict.keys())
        masses.sort()
        result.mass_list = masses
        result.intensity_list = [input_dict[x] for x in masses]
        return result

    @staticmethod
    def convert_from_lists(mass_list, intensity_list) -> 'MassSpectrum':
        result = MassSpectrum(None)
        result.mass_list = mass_list
        result.intensity_list = intensity_list
        return result

    @staticmethod
    def _dict_align(full_keys, input_dict):
        result = []
        for single_key in full_keys:
            result.append(input_dict.get(single_key, 0.0))
        return result

    @staticmethod
    def single_compare(intensity_left: list, intensity_right: list) -> float:  # numpy.typing.NDArray
        if DEBUG_STATS:
            print("Mass Spec relation")
            print(scipy.stats.pearsonr(intensity_left, intensity_right))    # Pearson's r
            print(scipy.stats.spearmanr(intensity_left, intensity_right))   # Spearman's rho
            print(scipy.stats.kendalltau(intensity_left, intensity_right))  # Kendall's tau
            print(scipy.stats.linregress(intensity_left, intensity_right))  # linear regression - see r2
        # normalize

        return scipy.stats.pearsonr(intensity_left, intensity_right)[0]
