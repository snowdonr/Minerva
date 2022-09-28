'''
Created on Oct 1, 2020

@author: Ryan
'''
from __future__ import annotations
import pyms.IonChromatogram
import pyms.GCMS.IO.ANDI
import pyms.GCMS.IO.JCAMP
import pyms.GCMS.IO.MZML
import pathlib
import pickle
import pickletools
import re
import numpy
import typing
import datetime
import collections
import logging
import openpyxl
import wx
import sqlalchemy.orm

from . import settings
from . import sql_interface
from gcms import peak_detect
from pyms.Utils.Time import time_str_secs


class SampleSet(collections.UserList):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.data is now empty

    def clear(self):
        self.data = []

    def load_directory(self, config: settings.Setting, progress_bar: wx.Gauge=None):
        directory_path = pathlib.Path(config.data_path)
        sample_relation_path = pathlib.Path(config.sample_relation_path)
        internal_standard_path = pathlib.Path(config.internal_standard_path)
        if not internal_standard_path.is_file():
            logging.warning(f"Warning: No weights found at {internal_standard_path}, samples will not be normalized!")
        sample_path_list = []
        self._load_subdirectory_recursive(directory_path, sample_path_list)
        if progress_bar is not None:
            progress_bar.SetRange(len(sample_path_list))
            progress_bar.SetValue(0)
            progress_value = 0

        print(f"STATUS: Finished scanning {directory_path}, checking against excel sheet {sample_relation_path}")
        excel_map = ExcelMap(config)
        excel_map.file_report(sample_path_list)

        print(f"Beginning load of {len(sample_path_list)} files")
        print(f"{datetime.datetime.now()}")
        pickle_samples = []
        pickle_paths = [self._get_pickled_path(path, config) for path in sample_path_list if path.suffix.lower() != ".pickle"]
        if config.read_pickle and len(pickle_paths) > 0:
            for pickle_path in pickle_paths:
                if pickle_path.is_file():
                    print(f"Opening: {pickle_path}")
                    try:
                        with open(pickle_path, "rb") as pf:
                            loaded_sample = pickle.load(pf)
                            if loaded_sample.sample_number is None:
                                loaded_sample.get_sample_number()  # Recheck file name in case regex is updated
                            if config.limit_loading_memory:
                                loaded_sample.memory_truncate(config.lo_rt_limit, config.hi_rt_limit)
                            pickle_samples.append(loaded_sample)
                    except Exception as _e:
                        logging.exception(f"Loading {pickle_path} failed, remove cached version")
                    if progress_bar is not None:
                        progress_value += 1
                        progress_bar.SetValue(progress_value)
                else:
                    logging.debug(f"No pickled entry at {pickle_path}")

        if config.force_load or len(pickle_samples) < len(sample_path_list):
            for sample_path in sample_path_list:
                pickle_path = self._get_pickled_path(sample_path, config)
                if not pickle_path.is_file():
                    new_sample = Sample(sample_path, config)
                    if new_sample.valid:
                        self.data.append(new_sample)  # .data from UserList implementation
                        if config.create_pickle:
                            print("Saving .pickle cache")
                            new_sample.save(pickle_path)
                elif config.overwrite_pickle:  # Assumed to be loaded in earlier loop
                    new_sample = Sample(sample_path, config)
                    if new_sample.valid:
                        self.data.append(new_sample)
                        new_sample.save(pickle_path)

                if progress_bar is not None:
                    progress_value += 1  # Might be able to exceed max if a pickle is reloaded
                    try:
                        progress_bar.SetValue(progress_value)
                    except Exception as _e:
                        pass

        self.data.extend(pickle_samples)

        if len(self.data) < 1:
            print("ERROR: No files found for processing")
            return False

        self._print_debug()

        for sample in self.data:
            sample.load_meta_from_excel(excel_map)

        print(f"STATUS: Finished loading {len(self)} samples")
        print(f"{datetime.datetime.now()}")
        if progress_bar is not None:
            progress_bar.SetRange(1)
            progress_bar.SetValue(1)
        return True

    def _get_pickled_path(self, sample_path: pathlib.Path, config: settings.Setting):
        return pathlib.Path(config.temp_path)/sample_path.with_suffix(".pickle").name

    def _print_debug(self):
        try:
            for sample in self.data:
                print(f"DEBUG: sample {sample.sample_number} mass range= {sample.min_mass}-{sample.max_mass}")
        except Exception as _e:
            logging.exception("DEBUGGING data failed, skipping")

    def _load_subdirectory_recursive(self, base_path: pathlib.Path, result_list: list):
        for item in base_path.iterdir():
            if item.is_dir():
                if item.suffix.lower() == ".d":
                    result_list.append(item)
                else:
                    self._load_subdirectory_recursive(item, result_list)
            else:
                result_list.append(item)

    def export_all_samples(self, directory: pathlib.Path):
        for sample_data in self.data:
            file_name = sample_data.path.stem + ".dat"
            # pyms.Utils.IO.save_data(directory / file_name, sample_data.im.intensity_array)
            sample_data.im.export_leco_csv(directory / file_name)

#     def load_dat_samples(self, directory: pathlib.Path, config: settings.Setting):
#         for sample_data_file in directory.iterdir():
#             if sample_data_file.suffix == ".dat":
#                 iim = pyms.IntensityMatrix.IntensityMatrix([0], [0], [[0]])
#                 iim.import_leco_csv(sample_data_file)
#                 new_sample = Sample(None, config)
#                 new_sample.path = sample_data_file
#                 new_sample._load_from_im(iim)
#                 self.data.append(new_sample)

    def get_tic(self, index: int):
        if index >= len(self):
            logging.warning(f"Could not plot tic, no sample at index {index}")
            return
        return self.data[index].get_tic()

    def remove_zero_peak_samples(self):
        self.data[:] = [x for x in self.data if len(x.peak_list) > 0]

#     def __init__(self, data=[]):
#         self._list = []
#         self.update(data)
#
#     def __getitem__(self, index):
#         return self._list[index]
#
#     def __setitem__(self, index, newValue):
#         self._list[index] = newValue
#
#     def __delitem__(self, index):
#         del self._list[index]
#
#     def __len__(self):
#         return len(self._list)
#
#     def insert(self, index, newItem):
#         self._list.insert(index, newItem)


class ExcelMap(object):
    '''
    Load data from excel file and initialize samples on request
    '''
    def __init__(self, config: settings.Setting):
        self.config = config
        excel_path = pathlib.Path(config.sample_relation_path)
        weight_path = pathlib.Path(config.internal_standard_path)
        self.base_sample_list = {}  # Map of sample numbers to list of derived samples. Not time efficient if # of samples gets large
        self._standard_headers = {}
        if not excel_path.is_file():
            logging.warning(f"Could not open excel file {excel_path}")
            return
        wb = openpyxl.load_workbook(excel_path)
        ws = wb.active  # Assumes only 1 worksheet
        for row in ws.rows:
            self._read_row(row)

        if weight_path.is_file():
            try:
                wbw = openpyxl.load_workbook(weight_path)
                wsw = wbw.worksheets[0]  # Weights
                for row in wsw.rows:
                    self._read_weight_row(row)
                wssw = wbw.worksheets[1]  # Standard weights
                for row in wssw.rows:
                    self._read_standard_weight_row(row)
            except Exception as _e:
                logging.exception('ERROR: Failed to read weight spreadsheet')

    def _read_row(self, row: tuple):
        if isinstance(row[0].value, str):  # row[0].data_type == 's' or 'n'
            return
        new_entry = ExcelEntry()
        new_entry.base_sample = row[0].value
        new_entry.derived_sample = row[1].value
        new_entry.dose = row[2].value
        new_entry.comment = row[3].value
        try:
            new_entry.sample_weight = row[4].value
        except Exception as _e:
            new_entry.sample_weight = None
        if new_entry.base_sample not in self.base_sample_list:
            self.base_sample_list[new_entry.base_sample] = [ExcelEntry(base_sample=new_entry.base_sample, dose=0)]

        self.base_sample_list[new_entry.base_sample].append(new_entry)

    def _read_weight_row(self, row: tuple):
        ''' Read weight spreadsheet and update sample_weight attribute for matching items in base_sample_list '''
        try:
            weight = float(row[2].value)
        except (TypeError, ValueError) as _e:
            logging.info(f"Ignoring weight row: {row[0].value}, {row[1].value}, {row[2].value}")
            return  # Assume header or empty row
        try:
            sample_number = int(row[1].value)
        except ValueError as _e:
            sample_number = int(re.findall(r'\d+', row[1].value)[0])
            repeat = True
        else:
            repeat = False
        existing_entry = self.find_data(sample_number)
        if existing_entry:
            if existing_entry.sample_weight is not None:
                if repeat:
                    logging.warning(f"Weight row repeat data {row[1].value} is currently ignored")
                else:
                    logging.warning(f"Found duplicated weight information at {row[1].value}, {weight} ignored")
            else:
                existing_entry.sample_weight = weight

    def _read_standard_weight_row(self, row: tuple):
        if row[0].value == 'Date':  # TODO: Improve spreadsheet
            self._setup_standard_headers(row)
            return
        elif not self._standard_headers or row[1] is None or not row[1].value:
            return
        base_sample_id = int(row[1].value.split(",")[0])  # TODO: Improve spreadsheet logic
        weights = {}
        for cmpd_name, column_position in self._standard_headers.items():
            weights[cmpd_name] = row[column_position].value
        for entry in self.base_sample_list[base_sample_id]:
            entry.standard_weights = weights

    def _setup_standard_headers(self, row: tuple):
        standard_expected_names = [x[0] for x in self.config.standard_list]
        self._standard_headers = {}
        for index, cell in enumerate(row):
            if cell.value in standard_expected_names:
                self._standard_headers[cell.value] = index

    def full_sample_list(self) -> list:
        # TODO: Calculate and store immediately?
        result = []
        count_none = 0  # Not included in output, report this?
        for data_entry_list in self.base_sample_list.values():
            for data_entry in data_entry_list:
                if data_entry.derived_sample is not None:
                    result.append(data_entry.derived_sample)
                else:
                    count_none += 1
        return result

    def file_report(self, data_file_list: list):
        '''
        Check for the presence of pickle or data files for all sample numbers given in the excel sheet
        '''
        used_file_dict = collections.defaultdict(list)
        unknown_file_list = []
        excel_missing_list = []
        total_found_count = 0
        for data_file_name in data_file_list:
            sample_number = Sample.extract_sample_number(data_file_name.stem)
            if sample_number is None:
                unknown_file_list.append(data_file_name)
            else:
                entry = self.find_data(sample_number)
                if entry is None:
                    print(f"Could not find any data for {sample_number}")
                    excel_missing_list.append(sample_number)
                elif entry.dose is not None:
                    entry.data_file = data_file_name
                    used_file_dict[sample_number].append(entry)
                    total_found_count += 1
                else:
                    excel_missing_list.append(sample_number)
        data_file_missing_list = set(self.full_sample_list()).difference(used_file_dict.keys())

        print("**File Report**")
        print(f"Samples missing excel entries: {excel_missing_list}")
        print(f"Files without sample data: {unknown_file_list}")
        print(f"Excel entries missing files: {data_file_missing_list}")
        print(f"**File matches: {len(used_file_dict)} samples with {total_found_count} files **")
        for key, entry_list in used_file_dict.items():
            print(f"Sample: {key}")
            for entry in entry_list:
                file_size = entry.data_file.stat().st_size
                if file_size < 7000000:
                    file_type = "SIM "
                else:
                    file_type = "SCAN"
                print(f"     {file_type}  {entry.data_file}  {entry.dose}  {entry.base_sample}  {entry.comment}")
        print("**End File Report**")

    def find_data(self, target_sample_number: int) -> typing.Tuple[int, float, str]:
        if target_sample_number in self.base_sample_list.keys():
            return self.base_sample_list[target_sample_number][0]  # The vase sample is always the first entry in the list, with empty values
        else:
            for data_list in self.base_sample_list.values():
                for entry in data_list:
                    if entry.derived_sample == target_sample_number:
                        return entry
        return None


class ExcelEntry(object):
    # TODO: Slots
    def __init__(self, base_sample=None, dose=None):
        self.base_sample = base_sample
        self.derived_sample = None
        self.dose = dose
        self.comment = None
        self.sample_weight = None
        self.standard_weights = {}
        self.repeat_run = False
        self.data_file = None


class Sample(sql_interface.Base):
    '''
    Data store for all information related to a GCMS run, including integrated peaks and meta-data
    '''
    __tablename__ = "Lab_Sample_Analysis"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    gc_peak = sqlalchemy.orm.relationship("Peak", back_populates="source_sample")
    project_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey("Minerva_Project.id"))
    project = sqlalchemy.orm.relationship("Project", back_populates="samples")

    sample_number = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    source_path = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    comment = sqlalchemy.Column(sqlalchemy.String, nullable=True)

    def __init__(self, raw_data_path: pathlib.Path, config: settings.Setting):
        super().__init__()
        self.path = raw_data_path  # pathlib.Path, should be static/update source_path
        self.source_path = str(raw_data_path)  # SQL output corresponding to self.path
        self.config = config
        self._raw_data = None  # pyms.GCMS_data, possibly not available if loaded from prior runs
        self.im = None  # pyms.IonChromatogram
        self.sample_number = None
        self.base_sample_number = None
        self.dose = None
        self.comment = None
        self.sample_weight = None
        self.standard_weights = {}
        self.get_sample_number()
        self.peak_list = []
        self.peak_source = peak_detect.Parse_IM(self)
        self.min_mass = None
        self.max_mass = None
        self._tic = None
        if raw_data_path is not None:
            if True:
                import cProfile
                import pstats
                profileName = 'profile.out'
                cProfile.runctx("self.valid = self._load_file(raw_data_path)", globals(), locals(), profileName)
                try:
                    p = pstats.Stats(profileName)
                    p.strip_dirs().sort_stats('time').print_stats(30)
                    p.strip_dirs().sort_stats('cumulative').print_stats(40)
                    # p.strip_dirs().sort_stats('tottime').print_stats(80)
                except Exception as _e:
                    print("Unable to save debug profile")
            else:
                self.valid = self._load_file(raw_data_path)

    @staticmethod
    def get_area(peak, replacement=0.0):
        try:
            value = peak.area
        except AttributeError as _e:
            value = replacement  # numpy.NaN, If there is no peak for a sample, input to this is: peak==None
        return value

    @classmethod
    def extract_sample_number(cls, file_string: str) -> int:
        file_name_match = re.match(r"(PRG|W\d+)-([^-]+).*", file_string)
        if file_name_match is not None:
            return file_name_match.group(2)
        else:
            return file_string

    def get_sample_number(self):
        self.sample_number = self.extract_sample_number(self.path.stem)
        if self.sample_number is None:
            logging.error(f"Could not parse sample id from {self.path}")

    def _load_file_type(self, path_input: pathlib.Path):
        data = None
        try:
            if path_input.suffix.lower() == ".cdf":
                data = pyms.GCMS.IO.ANDI.ANDI_reader(path_input)
            if path_input.suffix.lower() == ".jdx":
                data = pyms.GCMS.IO.JCAMP.JCAMP_reader(path_input)
            if path_input.suffix.lower() == ".mzml":
                data = pyms.GCMS.IO.MZML.mzML_reader(path_input)
        except Exception as _e:
            logging.exception(f"Failed to read file {path_input}")
        return data

    def _load_file(self, path_input: pathlib.Path) -> bool:
        ''' Load the first file in directory given by path_input (pathlib.Path) that has a known format and updates all object state parameters '''
        data = None
        if path_input.is_dir():
            for item in path_input.iterdir():
                if item.is_file():
                    if self.config.read_pickle and path_input.suffix.lower() == ".pickle":
                        logging.warning(f"Attempting to load pickle file at {path_input}")
                    else:
                        data = self._load_file_type(item)
                        if data is not None:
                            break
        elif path_input.is_file():
            data = self._load_file_type(path_input)
        else:
            logging.error(f"Unexpected sample target type {path_input}")
            return False

        if data is None:
            logging.error(f"Could not load any compatible data at {path_input}")
            return False
        self._raw_data = data
        # data.trim(4101, 4350)
        self.im = pyms.IntensityMatrix.build_intensity_matrix_i(data)  # default bin -0.3 +0.7
        self.min_mass = self.im.min_mass
        self.max_mass = self.im.max_mass
        self.peak_list = self.peak_source.load_from_im(self.im)
        return True

    def get_tic(self) -> pyms.IonChromatogram.IonChromatogram:
        if not hasattr(self, '_tic') or self._tic is None:  # TODO: hasattr can be removed once cache files are updated
            intensity_matrix = self.im
            intensity_array = sum(intensity_matrix.intensity_array.T)
            self._tic = pyms.IonChromatogram.IonChromatogram(intensity_array, intensity_matrix.time_list)
        return self._tic

#     def get_scan_index(self, peak):
#         pass

    def memory_truncate(self, low_limit, high_limit):
        self.min_mass = self.im.min_mass
        self.max_mass = self.im.max_mass
        self.get_tic()
        self.im = None
        self._raw_data = None
        rt_low = time_str_secs(low_limit)
        rt_high = time_str_secs(high_limit)
        new_peak_list = [x for x in self.peak_list if rt_low < x.rt < rt_high]
        self.peak_list = new_peak_list

    def mass_spec(self, scan_index: int):
        return self.im.get_ms_at_index(scan_index)

    def load_meta_from_excel(self, excel_map: ExcelMap):
        ''' Updates self.base_sample_number self.dose and self.comment '''
        if self.sample_number is None:
            return False
        map_result = excel_map.find_data(self.sample_number)
        if map_result is None:
            print(f"No excel data for sample {self.sample_number}")
            return
        self.base_sample_number = map_result.base_sample
        self.dose = map_result.dose
        self.comment = map_result.comment
        self.sample_weight = map_result.sample_weight
        self.standard_weights = map_result.standard_weights

    def get_filtered_peaks(self, peak_filter):
        try:
            return [peak for peak in self.peak_list if peak_filter(peak)]
        except Exception as _e:
            logging.exception("Failed to load peak areas")

    def get_standard_peaks(self):
        result = {}
        for peak_name, fragment_mass, estimated_rt in self.config.standard_list:
            if estimated_rt is None:
                continue  # Skip this compound, the peak isn't reliable enough
            best_peak = None
            best_rt_error = self.config.standard_maximum_rt_error+1
            best_intensity = 0
            for peak in self.peak_list:
                fragment_intensity = peak.mass_spectrum.get_intensity_for_mass(fragment_mass)
                if fragment_intensity > self.config.standard_minimum_intensity:
                    rt_error = abs(estimated_rt-peak.rt)
                    if rt_error < best_rt_error:
                        best_peak = peak
                        best_rt_error = rt_error
                        best_intensity = fragment_intensity
                    if peak.rt > estimated_rt+best_rt_error:  # peaks are in sorted order by rt, and we've gone past
                        break
            result[peak_name] = (best_peak, best_rt_error)  # Peak names must be unique
        return result

    def get_standard_weights(self, standard_name_list):
        return [self.standard_weights[x] for x in standard_name_list]

    def get_mass_sums(self, min_mass, max_mass):
        intensity_matrix = self.im.intensity_matrix
        filtered_range = (min_mass <= numpy.array(self.im.mass_list)) & (numpy.array(self.im.mass_list) <= max_mass)
        return numpy.sum(intensity_matrix[:, filtered_range], axis=0)

    def save(self, target_path: pathlib.Path):
        with open(target_path, 'wb') as file_out:
            _temp_source = self.peak_source
            self.peak_source = None
            opt_pickle = pickletools.optimize(pickle.dumps(self))
            file_out.write(opt_pickle)
            # self.peak_source = _temp_source

    def __eq__(self, other):
        return self.im == other.im and self.peak_list == other.peak_list

    def __str__(self):
        return f"{self.path.stem} with {len(self.peak_list)} peaks"

    def __repr__(self):
        try:
            return f"{self.path}:{len(self.peak_list)} peaks"
        except Exception as _e:
            return f"{self.sample_number}"
