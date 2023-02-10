'''
Created on Oct 4, 2020

@author: Ryan
'''
from __future__ import annotations
import numpy
import scipy.stats
import pathlib
import pickle
import logging
import csv
import gzip

from . import settings
from . import sample
from . import align
from . import sql_interface
from . import peak_detect
from . import mass_spectrum
from . import identification


class Analysis_Set(object):
    ''' Currently supports a single analysis for alignment and processing, may be extended for multiple '''
    def __init__(self, test_sample_set: sample.SampleSet, analysis_name: str, config: settings.Setting, sql: sql_interface.SQL_Interface):
        self._sample_set = test_sample_set
        self.config = config
        self.sql = sql
        self._results = []
        self.feature_plots = None
        self.method_analysis_name = analysis_name
        self._alignment_data = None
        self.alignment_data_groups = None

        cache_directory = pathlib.Path(self.config.temp_path)
        pickle_path = cache_directory / (str(analysis_name)+'_align.pickle')
        pickle_zip_path = cache_directory / (str(analysis_name)+'_align.pickle.cmp')
        if self.config.load_align_pickle and pickle_path.is_file():
            logging.info(f"Loading aligned values from {pickle_path}")
            with open(pickle_path, "rb") as pf:
                self._alignment_data = pickle.load(pf)
                self._alignment_data.pickle_load_finish(test_sample_set)
        elif self.config.load_align_pickle and pickle_zip_path.is_file():
            logging.info(f"Loading aligned values from zipped {pickle_path}")
            with gzip.open(pickle_zip_path, "rb") as pf_cmp:
                self._alignment_data = pickle.load(pf_cmp)
                self._alignment_data.pickle_load_finish(test_sample_set)
        else:
            logging.info(f"Beginning alignment of peaks for {analysis_name}")
            output_directory = pathlib.Path(self.config.analysis_directory)
            self._alignment_data = align.Align(analysis_name, [], test_sample_set, cache_directory, output_directory, self.config)

    def merge_aligned_peaks(self):
        '''
        Alignment pass outside of pymassspec library to check for suspect alignment failures
        Significant merging at this step should be taken as an error and suggest a rerun with
        different settings.
        '''
        # Note that peaks can be included in multiple aligned sets!
        aligned_data = self._alignment_data.peaks_aligned
        grouped_data = []
        for row_index in range(aligned_data.shape[0]):
            aligned_set = aligned_data.iloc[row_index, :]
            grouped_data.append(align.AlignedPeakSet(aligned_set, self.config, self._sample_set))

        candidate_rt_groups = []
        for outer_index, outer_entry in enumerate(grouped_data):
            start_average_rt = outer_entry.average_rt()
            current_group = []
            for aligned_item in grouped_data[outer_index:]:  # includes outer_entry as the first entry
                if aligned_item.average_rt() > start_average_rt + self.config.consider_rt_range:
                    break
                current_group.append(aligned_item)
            if len(current_group) > 1:
                candidate_rt_groups.append(current_group)

        for rt_group in candidate_rt_groups:
            target_aligned_peaks = rt_group[0]
            first_group = target_aligned_peaks.merged_ion_areas()
            related = []
            for other_group in rt_group[1:]:
                if other_group == target_aligned_peaks:
                    logging.debug(f"Duplicate entry in rt group {other_group}")
                other_peaks = [x for x in other_group if x is not None]
                other_area_list = [x.ion_areas for x in other_peaks]
                mapping_scores = self.compare_ion_areas(first_group, other_area_list)
                if numpy.mean(mapping_scores) > 0.7:
                    related.append(other_group)
            target_aligned_peaks.relations = related

        print(f"Merging alignment data {len(grouped_data)}")
        self.merged_data = []
        try:
            for grouped_item in grouped_data:
                if grouped_item.considered:
                    continue
                new_item = grouped_item.merge_relations()
                self.merged_data.append(new_item)
        except Exception as _e:
            # Triggered by memory reduction
            self.merged_data = grouped_data

        for grouped_item in grouped_data:
            aps = sql_interface.AlignedPeakSet_SQL(average_rt=grouped_item.average_rt(), total_area=grouped_item.total_area(), count=grouped_item.count)
            grouped_item.sql_ref = aps
            for pyms_peak in grouped_item:
                if pyms_peak is not None:
                    new_peak = peak_detect.Peak(pyms_peak)
                    new_peak.set_ions(pyms_peak.ion_areas)
                    aps.gc_peaks.append(new_peak)
        self.alignment_data_groups = grouped_data
        print(f"Merge complete, changed {len(self.merged_data)}  :  {len(self.alignment_data_groups)}")

    def compare_ion_areas(self, first: dict, other_list: list[dict]) -> list[float]:
        ''' Compares ion area dicts in a list against a target set and returns a list of scores matching the other_list '''
        result_scores = []
        for comparison_item in other_list:
            result_scores.append(self._single_compare_ion_dict(first, comparison_item))
        return result_scores

    def _single_compare_ion_dict(self, left: dict, right: dict) -> float:
        key_set = set(left).union(right)
        intensity_left = []
        intensity_right = []
        for key in key_set:
            intensity_left.append(left.get(key, 0.0))
            intensity_right.append(right.get(key, 0.0))
        if self.config.merge_compare_normalize:
            adjust_left = numpy.array(intensity_left)/max(intensity_left)
            adjust_right = numpy.array(intensity_right)/max(intensity_right)
        else:
            adjust_left = intensity_left
            adjust_right = intensity_right
        return scipy.stats.pearsonr(adjust_left, adjust_right)[0]  # normalization is not required

    def add_library_matches(self, source_library: identification.ID_Library, sql_session: sql_interface.SQL_Interface, use_ion_areas: bool =True):
        existing_db_entries = []
        max_rt_difference = self.config.rt_sensitivity_s * self.config.library_rt_scale
        try:
            existing_db_entries = sql_session.search_compounds()
        except Exception as _e:
            logging.warning(f"Could not open existing compound entries from {sql_session}")
        existing_compounds = {}
        for row in existing_db_entries:
            entry = row['Compound_SQL']
            existing_compounds[entry.key] = entry
        for aligned_set in self.alignment_data_groups:
            peak_data = aligned_set.highest_peak()
            if use_ion_areas:
                ms = mass_spectrum.MassSpectrum.convert_from_dict(peak_data.ion_areas)
            else:
                ms = peak_data.mass_spectrum
            target_rt = aligned_set.average_rt()
            if ms.mass_spec is None:
                logging.error(f"No mass spec for {ms}")
                continue
            found_reference_data = source_library.search(ms, target_rt, max_hits=5)
            if found_reference_data is None or len(found_reference_data) < 1:
                continue

            # self.mass_plot.show_peak_mass_spectrum(selected_peak, selected_sample, normalize=self.NORM_MAX, first_plot=True)
            library_ms_list = [x[1].mass_spec for x in found_reference_data]
            compare_scores_ms = mass_spectrum.MassSpectrum.compare_to_target_ms(peak_data.mass_spectrum, library_ms_list)
            try:
                library_rt_list = 1-(numpy.array([abs(x[1].rt-target_rt) if x[1].rt is not None else numpy.NaN for x in found_reference_data])/max_rt_difference)
            except AttributeError as _e:
                library_rt_list = numpy.zeros(len(library_ms_list))*numpy.nan
                # logging.exception("RT compare failed")

            for (_cmpd_name, ref_data), match_score, rt_difference in zip(found_reference_data, compare_scores_ms, library_rt_list):
                if not numpy.isnan(rt_difference):
                    current_match_type = "rt&ms"
                    match_score = max(rt_difference, 0) * match_score
                else:
                    current_match_type = "ms"
                if match_score < self.config.minimum_compound_score:
                    continue
                new_matcher = sql_interface.CompoundMatch_SQL(score=match_score, source=source_library.type.name, match_type=current_match_type)
                new_matcher.compound = ref_data
                aligned_set.add_match(new_matcher)  # sql_ref.compound_match.append(new_matcher)

    def write_aligned_csv(self, input_data, filepath: pathlib.Path, output_compounds: bool =False):
        '''
        Create an output beyond pymassspec library to summarize results, including compounds
        '''
        max_compound_output = 10
        name_skip_length = len(self.method_analysis_name)+1
        sample_names = [x[name_skip_length:] for x in self._alignment_data.peaks_aligned.columns]

        field_names = ["UID", "RT_s", "RT_m", "Mass1", "Portion1", "Mass2", "Portion2", "Mass3", "Portion3",
                       "Ion M1", "Ion I1", "Ion M2", "Ion I2", "Ion M3", "Ion I3", "Total area", "Count", "Related RT"]  # see AlignedPeakSet.get_output_dict
        field_names.extend(sample_names)
        if output_compounds:
            for index in range(max_compound_output):
                field_names.extend([f"C{index} name", f"C{index} formula", f"C{index} mw", f"C{index} score"])

        with open(filepath, 'w', newline='') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=field_names, dialect="excel")
            writer.writeheader()
            for aligned_peaks in input_data:
                data_dict = aligned_peaks.get_output_dict(sample_names, output_compounds, max_compound_output)
                writer.writerow(data_dict)
