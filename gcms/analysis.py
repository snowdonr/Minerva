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
        if self.config.load_align_pickle and pickle_path.is_file():
            logging.info(f"Loading aligned values from {pickle_path}")
            with open(pickle_path, "rb") as pf:
                self._alignment_data = pickle.load(pf)
                self._alignment_data.pickle_load_finish(test_sample_set)
        else:
            logging.info(f"Beginning alignment of peaks for {analysis_name}")
            output_directory = pathlib.Path(self.config.analysis_directory)
            self._alignment_data = align.Align(analysis_name, [], test_sample_set, cache_directory, output_directory, self.config)

    def merge_aligned_peaks(self):
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

        print(f"test point begin {len(grouped_data)}")
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
                    aps.gc_peaks.append(new_peak)
        self.alignment_data_groups = grouped_data
        print(f"test point end {len(self.merged_data)}  :  {len(self.alignment_data_groups)}")

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
                cmpd_key = sql_interface.Compound_SQL.dict_key(ref_data)
                try:
                    flavor_string = ref_data.odor
                except AttributeError as _e:
                    flavor_string = ""  # absent for non-user entries, i.e. NIST data
                if cmpd_key in existing_compounds:
                    new_compound = existing_compounds[cmpd_key]
                else:
                    new_compound = sql_interface.Compound_SQL(name=ref_data.name, formula=ref_data.formula, mw=ref_data.mw, cas=ref_data.cas, flavors=flavor_string)
                    existing_compounds[cmpd_key] = new_compound
                new_matcher = sql_interface.CompoundMatch_SQL(score=match_score, source=source_library.type.name, match_type=current_match_type)
                new_matcher.compound = new_compound
                aligned_set.sql_ref.compound_match.append(new_matcher)

#     def scan_for_ms_matches(self, aligned_data):
#         consider_mass_count = 3
#         consider_mass_fraction = 0.3
#         consider_rt_range = 5
#         candidate_mass_groups = collections.defaultdict(list)
#         for row_index in range(aligned_data.shape[0]):
#             aligned_set = aligned_data.iloc[row_index, :]
#             peaks = [x for x in aligned_set if x is not None]
#             if len(peaks) < 1:
#                 return None
#             exemplar_peak = peaks[0]
#             mass_spectrum_list = [x.mass_spectrum for x in peaks]
#             merged_values = mass_spectrum.MassSpectrum.merge_mass_spectrum_list(mass_spectrum_list)
#             if merged_values is None:
#                 continue
#             average_rt = aligned_set.average_rt
#             intensity_sort = numpy.argsort(merged_values.intensity_list)
#             for index in range(1, consider_mass_count+1):
#                 mass_index = intensity_sort[-index]
#                 mass_value = merged_values.mass_list[mass_index]
#                 candidate_mass_groups[mass_value].append((average_rt, merged_values, aligned_set))
#
#         candidate_mass_rt_groups = []
#         for mass_related_set in candidate_mass_groups.values():
#             if len(mass_related_set) < 2:
#                 continue
#             current_group = []
#             last_rt = -numpy.Inf
#             for average_rt, merged_ms, set_list in mass_related_set:
#                 if average_rt-last_rt < consider_rt_range:
#                     current_group.append((average_rt, merged_ms, set_list))
#                 elif len(current_group) > 0:
#                     candidate_mass_rt_groups.append(current_group)
#                     last_rt = average_rt  # may get updated if partial group is used forward
#                     for index in range(1, len(current_group)):
#                         previous_entry_rt = current_group[index][0]
#                         if average_rt-previous_entry_rt < consider_rt_range:
#                             current_group = current_group[index:]
#                             last_rt = previous_entry_rt
#                             break
#                 else:
#                     last_rt = average_rt

    def write_aligned_csv(self, input_data, filepath: pathlib.Path, output_compounds: bool =False):
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


#         dose_list_y = numpy.array(dose_list_y, dtype=float)
#
#         self.feature_sources = []
#         if self.config.analyze_mass_sums:
#             feature_data_X = self.raw_data.mass_grid
#             mass_values = range(self.raw_data.mass_range[0], self.raw_data.mass_range[1]+1)
#             keyed_data = [self.raw_data.mass_grid[:, i] for i in range(self.raw_data.mass_grid.shape[1])]
#             self.feature_sources.extend(zip(mass_values, keyed_data))
#
#         if self.config.analyze_peak_areas:
#             get_area = numpy.frompyfunc(sample.Sample.get_area, 1, 1)  # operator.attrgetter("area")
#             all_peak_areas = get_area(self.raw_data.peaks_aligned.values)
#             peak_array = numpy.array(all_peak_areas.T-self.config.filter_minimum_area, dtype=float)  # Assumes areas are positive, with minimum as filter_minimum_area
#             peak_array[peak_array < 0.0] = 0.0  # i.e. don't subtract from 0 values, that would defeat the point
#
#             if self.config.analyze_mass_sums:
#                 feature_data_X = numpy.hstack((feature_data_X, peak_array))
#             else:
#                 feature_data_X = peak_array
#
#             feature_peak_sets = []
#             for peak_index in range(self.raw_data.peaks_aligned.values.shape[0]):
#                 peak_set = self.raw_data.peaks_aligned.values[peak_index, :]
#                 feature_peak_sets.append(peak_set)
#                 # areas = [sample.Sample.get_area(peak) for peak in peak_set]
#             self.feature_sources.extend(feature_peak_sets)
#
#         sample_weights = [x.sample_weight for x in self.samples]
#         normalization_values_by_sample = None
#         if self.config.divide_samples:
#             if self.config.divide_type == 'weights_only':
#                 if not numpy.isnan(sample_weights).any():
#                     normalization_values_by_sample = 1/numpy.array(sample_weights)
#             elif self.config.divide_type == 'sample_total':
#                 normalization_values_by_sample = 1/numpy.sum(feature_data_X, axis=1)
#             elif self.config.divide_type == 'standard_area':
#                 standard_set = collections.defaultdict(list)
#                 for current_sample in self.samples:
#                     standard_data = current_sample.get_standard_peaks()
#                     for key, (peak_data, rt_error) in standard_data.items():
#                         standard_set[key].append(peak_data)
#                 standard_names = []
#                 standard_norm_set = []
#                 for key, peak_list in standard_set.items():
#                     if None in peak_list:
#                         print(f"Missing standard peak for {key}: {peak_list}")
#                     else:
#                         standard_names.append(key)
#                         added_list = [peak.area for peak in peak_list]
#                         if None in added_list:
#                             logging.error(f"No area for normalization standard {key}")
#                         else:
#                             standard_norm_set.append(added_list)
#
#                 if len(standard_names) > 0:
#                     print(f"Normalizing using the following standards: {standard_names}")
#                     twod_standard_areas = numpy.vstack(standard_norm_set)
#                     standard_weights = [sum(samp.get_standard_weights(standard_names)) for samp in self.samples]
#                     standard_area_sum = numpy.sum(twod_standard_areas, axis=0)
#                     normalization_values_by_sample = standard_weights/(standard_area_sum*sample_weights)
#             else:
#                 print("Unknown divide_type in Setting: {Setting.divide_type}")
#             if normalization_values_by_sample is not None:
#                 feature_data_X *= numpy.tile(normalization_values_by_sample, (feature_data_X.shape[1], 1)).T
#             else:
#                 print("Requested Normalization failed, check internal standards")
#         else:
#             print("Weight normalization disabled")
#
#         scaler = sklearn.preprocessing.StandardScaler(copy=False, with_mean=self.config.normalize_mean, with_std=self.config.normalize_std)
#         scaler.fit(feature_data_X)
#         feature_data_X = scaler.transform(feature_data_X)
#
#         if self.config.do_lasso:
#             try:
#                 lasso = sklearn.linear_model.Lasso(alpha=self.config.lasso_alpha)
#                 lasso.fit(feature_data_X, dose_list_y)
#                 self.save_all("_lasso", lasso)
#             except Exception as _e:
#                 logging.exception("LASSO")
#         if self.config.do_lars_lasso:
#             try:
#                 cv_lasso = sklearn.linear_model.LassoLarsCV(verbose=True, normalize=False, n_jobs=2)  # -1 for all processors
#                 cv_lasso.fit(feature_data_X, dose_list_y)
#                 print(cv_lasso.coef_path_)
#                 # TODO: get lars_path
#     #             predicted_values = cv_lasso.predict(Setting.lasso_prediction_target_dose)
#     #             print(predicted_values)
#                 self.save_regression(pathlib.Path(self.config.analysis_directory) / (self.name+"_lars_lasso_regression.pickle"), cv_lasso)
#             except Exception as _e:
#                 logging.exception("LARS LASSO")
#         if self.config.do_linear_svc:
#             try:
#                 clf = sklearn.pipeline.make_pipeline(sklearn.feature_selection.SelectFromModel(estimator=sklearn.linear_model.LogisticRegression()),
#                                                      sklearn.preprocessing.StandardScaler(), sklearn.svm.LinearSVC(random_state=0, tol=1e-5))
#                 clf.fit(feature_data_X, dose_list_y)
#                 self.save_regression(pathlib.Path(self.config.analysis_directory) / (self.name+"_svm_selection.pickle"), clf)
#             except Exception as _e:
#                 logging.exception("Linear SVC")
#         if self.config.do_bayesian_ridge:
#             try:
#                 bay_ridge = sklearn.linear_model.BayesianRidge(verbose=True, normalize=False)  # -1 for all processors
#                 bay_ridge.fit(feature_data_X, dose_list_y)
#                 self.save_regression(pathlib.Path(self.config.analysis_directory) / (self.name+"_bayesian_ridge_regression.pickle"), bay_ridge)
#             except Exception as _e:
#                 logging.exception("Bayesian Ridge")
#
#     def save_all(self, regression_name, regression_obj):
#         self.save_regression(pathlib.Path(self.config.analysis_directory) / (self.name+regression_name+"_regression.pickle"), regression_obj)
#         self.report_regression_results(regression_obj)
#         self.save_results(pathlib.Path(self.config.analysis_directory) / (self.name+regression_name+"_results.pickle"))
#
#     def save_regression(self, target_path: pathlib.Path, regression_obj):
#         with open(target_path, 'wb') as file_out:
#             pickle.dump(regression_obj, file_out)
#
#     def report_regression_results(self, regression_obj):
#         print(f"Results of {regression_obj}:  L1:{regression_obj.l1_ratio}")
#         response_coef = regression_obj.coef_
#         response_sort_keys = numpy.argsort(abs(response_coef))[::-1]
#         sample_numbers = [sample.sample_number for sample in self.samples]
#         print(f"Samples: {sample_numbers}")
#         dose_set = [sample.dose for sample in self.samples]
#         print(f"Doses:   {dose_set}")
#
#         feature_details = numpy.array(self.feature_sources)[response_sort_keys]
#         response_list = response_coef[response_sort_keys]
#         self.plots.load_all_features(self.samples, feature_details, response_list)
#         self.plots.update_feature_set(0)
#
#         with open(self.config.summary_file_path, "a") as file_out:
#             file_out.write(f"{self.name}\n")
#             for feature, response_coeff in zip(feature_details, response_list):
#                 try:
#                     detail_list = [str(x.rt) if x is not None else 'None' for x in feature]
#                 except Exception as _e:
#                     detail_list = [str(x) for x in feature]
#                 file_out.write(f"{response_coeff}, {','.join(detail_list)}\n")
#
#     def save_results(self, target_path: pathlib.Path):
#         with open(target_path, 'wb') as file_out:
#             temp_store = [self.feature_sources, [(x.sample_number, x.dose, x.path) for x in self.samples]]
#             pickle.dump(temp_store, file_out)
#         # Need dose and sample names from self.samples, plus self.feature_sources
