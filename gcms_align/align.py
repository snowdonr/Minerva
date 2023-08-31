'''
Created on Nov. 2, 2021

@author: Ryan
'''
from __future__ import annotations
import numpy
import pickle
import pathlib
import logging
import datetime
import functools
import collections
import copy
import gzip
import csv

import pyms.DPA.Alignment
import pyms.DPA.PairwiseAlignment
import pyms.Peak.List.IO
import pyms.Experiment
# import pandas  # Align object creates pandas DataFrame results

from . import sample
from . import settings
from . import mass_spectrum
from . import sql_interface


class Align(object):
    '''
    Used to align a set of samples
    '''

    def __init__(self, analysis_name: str, sample_set: sample.SampleSet, cache_directory: pathlib.Path, output_directory: pathlib.Path, config: settings.Setting):
        load_from_files = False  # unimplemented
        save_to_file = False
        expr_list = []
        self.config = config
        self.remaining_samples = []
        if len(sample_set) < 1:
            logging.error("No samples to align")

        if load_from_files:
            for file in output_directory.iterdir():
                expr = pyms.Experiment.load_expr(file)
                expr_list.append(expr)
        else:
            largest_min_mass = sample_set[0].min_mass
            smallest_max_mass = sample_set[0].max_mass

            removal_list = []
            for current_sample in sample_set:
                if len(current_sample.peak_list) < 1:
                    removal_list.append(current_sample)
            for sample in removal_list:
                pass

            for current_sample in sample_set:
                if current_sample.min_mass > largest_min_mass:
                    largest_min_mass = current_sample.min_mass
                if current_sample.max_mass < smallest_max_mass:
                    smallest_max_mass = current_sample.max_mass
                # file_name = current_sample.path.stem + ".dat"
                # current_sample.im.export_leco_csv(self._directory_path / file_name)
            self.mass_range = [largest_min_mass, smallest_max_mass]
            if self.config.minimum_mass_limit is not None:
                if largest_min_mass < self.config.minimum_mass_limit:
                    self.mass_range[0] = self.config.minimum_mass_limit
            if self.config.maximum_mass_limit is not None:
                if smallest_max_mass > self.config.maximum_mass_limit:
                    self.mass_range[1] = self.config.maximum_mass_limit

            skip_sample = []
            logging.info(f"Cropping all samples to a mass range of {self.mass_range[0]}-{self.mass_range[1]}")
            for current_sample in sample_set:
                if current_sample.im:
                    current_sample.im.crop_mass(self.mass_range[0], self.mass_range[1])
                removal_list = []
                for target_peak in current_sample.peak_list:
                    try:
                        target_peak.crop_mass(self.mass_range[0], self.mass_range[1])  # !! This change is in place and permanent  !!
                    except IndexError as _e:
                        removal_list.append(target_peak)
                for removed in removal_list:
                    logging.warning(f"No masses left after cropping for {removed.mass_spectrum}")
                    current_sample.peak_list.remove(removed)
                if len(current_sample.peak_list) < 1:
                    skip_sample.append(current_sample)
                    logging.warning(f"Skipping sample {current_sample}, no peaks remaining after mass crop")
            sample_set[:] = [x for x in sample_set if x not in skip_sample]

            logging.info(f"Done cropping mass, setting up experiments and selecting RT range")
            used_names = []
            for count, current_sample in enumerate(sample_set):
                # Create an Experiment

                experiment_peaks = current_sample.get_filtered_peaks(self.peak_filter)
                if not experiment_peaks or len(experiment_peaks) < 1:
                    logging.error(f"No accepted peaks in sample {current_sample}")
                    continue
                else:
                    self.remaining_samples.append(current_sample)
                    logging.info(f"Sample {current_sample} has {len(experiment_peaks)} used for alignment")
                experiment_name_base = analysis_name+"_"+str(current_sample.sample_number)
                experiment_name_candidate = experiment_name_base
                counter = 1
                while experiment_name_candidate in used_names:
                    experiment_name_candidate = experiment_name_base + "_" + str(counter)
                    counter += 1
                used_names.append(experiment_name_candidate)
                expr = pyms.Experiment.Experiment(experiment_name_candidate, experiment_peaks)
                # print(f"\t -> Selecting retention time range between '{lo_rt_limit}' and '{hi_rt_limit}'")
                expr.sele_rt_range([self.config.lo_rt_limit, self.config.hi_rt_limit])
                if save_to_file:
                    expr.dump(output_directory / (analysis_name+"_"+str(count)))
                expr_list.append(expr)

        if len(expr_list) < 1:
            logging.error(f"No samples for analysis from: {sample_set}")
            return

        print(f"Beginning alignment of {len(expr_list)} experiments")
        print(f"{datetime.datetime.now()}")
        temp_path = cache_directory / (str(analysis_name)+'_temp_align.pickle')
        if temp_path.is_file() and self.config.load_align_pickle:
            with open(temp_path, "rb") as pf:
                pair_aligned = pickle.load(pf)
        else:
            # Convert each Experiment object into an Alignment object with the function exprl2alignment()..
            align_sample_experiments = pyms.DPA.Alignment.exprl2alignment(expr_list)
            full_length = len(align_sample_experiments)
            align_sample_experiments = [x for x in align_sample_experiments if len(x.peakalgt) > 0]
            logging.info(f"Experiment list now {full_length} from {len(align_sample_experiments)}")
            pairwise_temp = pathlib.Path(self.config.temp_path) / pathlib.Path((str(analysis_name)+self.config.pairwise_temp_file))
            align_tree = pyms.DPA.PairwiseAlignment.PairwiseAlignment(align_sample_experiments, float(self.config.rt_sensitivity_s), float(self.config.gap_penalty), pairwise_temp, self.config)
            pair_aligned = pyms.DPA.PairwiseAlignment.align_with_tree(align_tree, min_peaks=self.config.alignment_minimum_peaks, config=self.config)
            try:
                with open(temp_path, 'wb') as file_out:
                    pickle.dump(pair_aligned, file_out)
            except Exception as _e:
                logging.exception("aligned pickle")

        try:
            pair_aligned.write_csv(output_directory / (str(analysis_name)+'_rt.csv'),  # for inspection
                                   output_directory / (str(analysis_name)+'_area.csv'))  # data matrix, samples-rows
        except Exception as _e:
            logging.exception("Could not save _rt or _area csv files")

        print(f"Done alignment of {len(expr_list)} experiments, storing details")
        common_ion_list = pair_aligned.common_ion()
        aligned_peaks = pair_aligned.aligned_peaks()
        try:
            pyms.Peak.List.IO.store_peaks(aligned_peaks, output_directory / (str(analysis_name)+'_aligned_peaks.bin'))
            pair_aligned.write_common_ion_csv(output_directory / (str(analysis_name)+'_area_common_ion.csv'), common_ion_list)
        except Exception as _e:
            logging.exception("Could not save _aligned_peaks.bin or _area_common_ion.csv")

        # ms_aligned = pair_aligned.get_ms_alignment(Setting.require_all_peaks)
        self.peaks_aligned = pair_aligned.get_peaks_alignment(self.config.require_all_peaks)

        print(f"Starting main alignment storage of {len(expr_list)} experiments")
        try:
            if self.config.compress_cache_files:
                cache_file_name = cache_directory / (str(analysis_name)+'_align.pickle.cmp')
                with gzip.open(cache_file_name, 'wb') as file_out:
                    pickle.dump(self, file_out)
            else:
                cache_file_name = cache_directory / (str(analysis_name)+'_align.pickle')
                with open(cache_file_name, 'wb') as file_out:
                    # self.peaks_aligned can also have sample references
                    sample_temp = self.remaining_samples
                    self.remaining_samples = None
                    pickle.dump(self, file_out)
                    self.remaining_samples = sample_temp
        except Exception as _e:
            logging.exception("Failed to save checkpoint during alignment")

        try:
            if self.config.comparison_export:
                self.comparison_export(sample_set, self.config.analysis_directory)
        except Exception as _e:
            logging.exception("Failed to create sample comparison")

        print(f"Done alignment storage of {len(expr_list)} experiments")
        print(f"{datetime.datetime.now()}")
        self.pickle_load_finish(sample_set)

    def peak_filter(self, peak_input: pyms.Peak.Peak) -> bool:
        if peak_input.area is None:
            return False
        elif peak_input.rt < self.config.manual_baseline_section_rt*60:
            return peak_input.area > self.config.filter_minimum_area
        else:
            return peak_input.area > 20000  # TODO: Update

    def _sample_set_filter(self, sample_set: sample.SampleSet) -> list:
        # TODO: Remove duplication

        filtered_samples = []
        for current_sample in sample_set:

            experiment_peaks = current_sample.get_filtered_peaks(self.peak_filter)
            if not experiment_peaks or len(experiment_peaks) < 1:
                logging.error(f"No accepted peaks in sample {current_sample}")
                continue
            else:
                filtered_samples.append(current_sample)
                logging.info(f"Sample {current_sample} has {len(experiment_peaks)} used for alignment")
        return filtered_samples

    def pickle_load_finish(self, sample_set: sample.SampleSet):
        self.sample_set = sample_set
        self._add_sample_ref(sample_set)
        if self.remaining_samples is None:
            self.remaining_samples = self._sample_set_filter(sample_set)

        if self.config.comparison_export:
            self.comparison_export(sample_set, self.config.analysis_directory)

    def _add_sample_ref(self, sample_set: sample.SampleSet):
        for aligned_row in self.peaks_aligned.values:
            for index, entry in enumerate(aligned_row):
                if entry is not None:
                    try:
                        entry.source_sample = sample_set[index]
                    except IndexError:
                        logging.warning("Source sample out of range")

    def comparison_export(self, sample_set: sample.SampleSet, output_directory: pathlib.Path):
        ''' Set-wise compare samples and create csv summary files '''
        sample_dict = {}
        aligned_size = None
        for sample, aligned_row in zip(sample_set, self.peaks_aligned.T.values):
            sample_dict[sample.sample_number] = aligned_row
            if aligned_size is None:
                aligned_size = len(aligned_row)
            elif aligned_size != len(aligned_row):
                logging.warning("Unequal comparison length")
        # sample_sets = itertools.combinations(sample_set, 2)
        sample_sets = [[sample] for sample in sample_set]
        pair_results = {}
        for pair in sample_sets:
            in_compounds = [None]*aligned_size
            out_compounds = [None]*aligned_size
            for test_sample in sample_set:
                if test_sample in pair:
                    in_compounds = self.union_peak_lists(sample_dict[test_sample.sample_number], in_compounds)
                else:
                    out_compounds = self.union_peak_lists(sample_dict[test_sample.sample_number], out_compounds)
            pair_name = "_".join([x.sample_number for x in pair])
            pair_results[pair_name] = self.difference_peak_lists(in_compounds, out_compounds)  # set .difference

        manual_set = self.difference_peak_lists(sample_dict["C1"], self.union_peak_lists(sample_dict["C2"], sample_dict["B1"]))
        self.write_csv(output_directory / ("C1-C2_B1"+".csv"), manual_set)

        manual_set = self.difference_peak_lists(sample_dict["C1"], sample_dict["C2"])
        self.write_csv(output_directory / ("C1-C2"+".csv"), manual_set)

        manual_set = self.difference_peak_lists(sample_dict["B1"], sample_dict["B2"])
        self.write_csv(output_directory / ("B1-B2"+".csv"), manual_set)

        for pair_name, pair_item in pair_results.items():
            self.write_csv(output_directory / (pair_name+".csv"), pair_item)

    def write_csv(self, file_name, data_set):
        ''' Create a csv from a list of Peaks '''
        if sum(1 for i in data_set if i is not None) == 0:
            print(f"Skipping {file_name}")
            return
        try:
            with open(file_name, 'w') as file_out:
                writer = csv.writer(file_out, dialect="excel")
                writer.writerow(["Peak UID", "Sample", "Area", "RT(min)", "Mass", "Ion Area", "Mass", "Ion Area"])
                for peak_item in data_set:
                    if self.config.comparison_skip_empty and peak_item is None:
                        writer.writerow("")
                    else:
                        mass_spec_list = []
                        for mass_value, ion_area in peak_item.ion_areas.items():
                            mass_spec_list.extend((mass_value, ion_area))
                        writer.writerow([peak_item.UID, peak_item.source_sample.sample_number, peak_item.area, peak_item.rt/60]+mass_spec_list)
        except Exception as _e:
            logging.exception("aligned pickle")

    def union_peak_lists(self, left_list, right_list):
        ''' Set union equivalent for positional lists of Peak/None '''
        result_list = [a or b for (a, b) in zip(left_list, right_list)]
        return result_list

    def difference_peak_lists(self, first_list, second_list):
        ''' Set difference equivalent for positional lists of Peak/None '''
        result_list = [a if a and not b else None for a, b in zip(first_list, second_list)]
        return result_list


def cache_ignore_args(input_function):
    @functools.wraps(input_function)
    def cached_function(*args, **kwargs):
        this = args[0]
        if input_function in this._cache:
            return this._cache[input_function]
        this._cache[input_function] = input_function(*args, **kwargs)
        return this._cache[input_function]
    return cached_function


class AlignedPeakSet(collections.UserList):
    ''' Static group of peaks created by an alignment
    sqlalchemy version at sql_interface.AlignedPeakSet_SQL '''

    def __init__(self, initList: list, config: settings.Setting, sample_set: sample.SampleSet):
        super().__init__(initList)
        self._cache = {}
        self.compound_matches = []
        self.config = config
        self.samples = sample_set
        self.relations = []
        self.considered = False
        self.sql_ref: sql_interface.AlignedPeakSet_SQL = None  # SQL object with this data

    def add_match(self, new_matcher):
        self.compound_matches.append(new_matcher)
        self.sql_ref.compound_match.append(new_matcher)

    @property
    def count(self) -> int:
        ''' Number of Peaks (not None) in this list '''
        peak_list = [x for x in self.data if x is not None]
        return len(peak_list)

    @property
    def single_peak(self) -> pyms.Peak.PeakClass.Peak:
        peak_list = [x for x in self.data if x is not None]
        if len(peak_list) > 0:
            return peak_list[0]
        else:
            return None

    # @functools.cached_property might be an alternate solution
    @cache_ignore_args
    def highest_peak(self) -> pyms.Peak.PeakClass.Peak:
        intensity_list = [x.area if x is not None else 0.0 for x in self.data]
        highest_peak_index = numpy.argmax(intensity_list)
        return self.data[highest_peak_index]

    @cache_ignore_args
    def average_rt(self) -> float:
        peak_list = self.data
        rt_values = [x.rt for x in peak_list if x is not None]
        if numpy.any(numpy.isclose(rt_values, 0.0)):
            logging.warning(f"Warning: Zero rt value in {peak_list}")
        mean_rt = numpy.mean(rt_values)
        if numpy.any(abs(rt_values-mean_rt)) > float(self.config.rt_sensitivity_s) / 2.0:
            logging.warning(f"Warning: Outlier rt value in {peak_list} Expected {mean_rt}: {rt_values}")
        return mean_rt

    @cache_ignore_args
    def merged_ms(self) -> ('MassSpectrum', numpy.array):
        ''' Mass spectrum average and sample mapping '''
        mass_spec_list = [x.mass_spectrum for x in self.data if x is not None]
        result = mass_spectrum.MassSpectrum.merge_mass_spectrum_list(mass_spec_list)
        return result

    @cache_ignore_args
    def merged_ion_areas(self) -> collections.Counter:
        mass_dict_list = [x.ion_areas for x in self.data if x is not None]
        result_counter = collections.Counter()
        for md in mass_dict_list:
            result_counter.update(md)
        return result_counter

    @cache_ignore_args
    def total_area(self) -> float:
        mass_dict_list = [x.ion_areas for x in self.data if x is not None]
        return sum([sum(x.values()) for x in mass_dict_list])

    def merge_relations(self) -> 'AlignedPeakSet':
        if not self.relations:
            self.considered = True
            return self
        merged_peak_list = self._recursive_merge(0)
        return AlignedPeakSet(merged_peak_list, self.config, self.samples)

    def _recursive_merge(self, count: int) -> list:
        if self.considered:
            return []
        if count >= self.config.max_relation_grouping_depth:
            return []
        result = copy.copy(self.data)
        self.considered = True
        for item in self.relations:
            add = self._check_ion_continuity(item)
            if add:
                result.extend(item._recursive_merge(count+1))
        return result

    def _check_ion_continuity(self, other: AlignedPeakSet, used_mass_count=3, minimum_factor=0.75) -> bool:
        ''' Checks to see if the ion *intensities* are always relatively high between two peaks
            Masses to check are based on the highest *area* ions of this set
            Compares only samples that have a peak in this set  '''
        previous_rt = self.highest_peak().rt
        current_rt = other.highest_peak().rt
        masses = numpy.array(list(self.merged_ion_areas().keys()))
        areas = numpy.array(list(self.merged_ion_areas().values()))
        selected_indices = numpy.argsort(areas)[::-1][:used_mass_count]
        mass_list = masses[selected_indices]
        intensity_mask = numpy.isin(self.highest_peak().mass_spectrum.mass_list, mass_list)
        if sum(intensity_mask) < used_mass_count:
            if sum(intensity_mask) < 1:
                logging.warning("No ions found for continuity check")
                return True
            new_mass_list = []
            for single_mass in mass_list:
                if single_mass in self.highest_peak().mass_spectrum.mass_list:
                    new_mass_list.append(single_mass)
            mass_list = numpy.array(new_mass_list)
        minimum_intensities = numpy.array(self.highest_peak().mass_spectrum.intensity_list)[intensity_mask]*minimum_factor  # consider merged
        sample_mapping = numpy.array(self.data) != None  # Cannot use "is None", we need the numpy != overload or a generator @IgnorePep8
        if len(self.samples.data) != len(sample_mapping):
            logging.error("Loaded inconsistent file counts, change data source or realign")
        samples_to_check = numpy.array(self.samples.data)[sample_mapping]
        result_list = []
        for used_sample in samples_to_check:
            result_list.append(self._check_sample_ions(used_sample, previous_rt, current_rt, mass_list, minimum_intensities))
        return numpy.mean(result_list) > 0.5

    def _check_sample_ions(self, sample_input, low_rt: float, high_rt: float, target_mass_list: list[int], minimum_intensities: list[float]) -> bool:
        rt_index_low = numpy.searchsorted(sample_input.im.time_list, low_rt)
        rt_index_high = numpy.searchsorted(sample_input.im.time_list, high_rt)
        mass_mask = numpy.isin(sample_input.im.mass_list, target_mass_list)
        ms_slice = sample_input.im.intensity_matrix[rt_index_low:rt_index_high, mass_mask]
        if sum(mass_mask) < len(target_mass_list):
            logging.error("Error: Mixed mass ranges csions")
        if len(ms_slice) < 1:
            return True
        return numpy.all(ms_slice.min(axis=0) > minimum_intensities)

    def get_output_dict(self, sample_names: list, output_compounds: bool =False, max_compound_output: int = 5) -> dict:
        ''' Creates a dictionary with the peaks areas per sample and calculated values '''
        result_dict = collections.defaultdict(lambda: 0.0)
        # ["UID", "RT_s", "RT_m", "Mass1", "Portion1", ... "Ion M1", "Ion I1", ... "Total area", "Count", "Related RT"]
        result_dict["UID"] = self.highest_peak().UID
        result_dict["RT_s"] = self.average_rt()
        result_dict["RT_m"] = result_dict["RT_s"]/60.0
        result_dict["Count"] = self.count
        try:
            result_dict["Related RT"] = " ".join([str(x.average_rt()) for x in self.relations])
        except AttributeError:
            pass

        try:
            mass_area_spec = collections.OrderedDict(self.highest_peak().ion_areas)  # Is forcing order sensitivity required?
            total_area = sum(mass_area_spec.values())
            mass_spec_masses = numpy.array(list(mass_area_spec.keys()))
            mass_spec_intensities = numpy.array(list(mass_area_spec.values()))
            intensity_pos = numpy.argsort(mass_spec_intensities)[::-1]
            for i in range(1, min(len(intensity_pos), 4)):
                target_index = intensity_pos[i-1]
                result_dict["Mass"+str(i)] = mass_spec_masses[target_index]
                result_dict["Portion"+str(i)] = mass_spec_intensities[target_index]/total_area

            mass_spec = self.merged_ms()[0]
            intensity_pos = numpy.argsort(mass_spec.intensity_list)[::-1]
            for i in range(1, min(len(intensity_pos), 4)):
                target_index = intensity_pos[i-1]
                result_dict["Ion M"+str(i)] = mass_spec.mass_list[target_index]
                result_dict["Ion I"+str(i)] = mass_spec.intensity_list[target_index]
        except Exception as _e:
            logging.exception("No mass spec data available")
        total_area = 0
        for index, entry in enumerate(self.data):
            if entry is None:
                continue
            if index >= len(sample_names):
                index = index % len(sample_names)
            result_dict[sample_names[index]] += entry.area
            total_area += entry.area
        result_dict["Total area"] = total_area

        if output_compounds and self.sql_ref:
            for index, entry in enumerate(self.sql_ref.compound_match):
                if index >= max_compound_output:
                    break
                result_dict[f"C{index} name"] = entry.compound.name.encode('ascii', 'ignore')
                result_dict[f"C{index} formula"] = entry.compound.formula
                result_dict[f"C{index} mw"] = entry.compound.mw
                result_dict[f"C{index} score"] = entry.score

        return result_dict
