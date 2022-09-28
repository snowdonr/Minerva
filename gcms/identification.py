'''

Created on Sep. 21, 2021

@author: Ryan
'''
import enum
import pyms.Spectrum
# import pyms_nist_search
import openpyxl
import sqlalchemy
import logging
import pathlib
import numpy

from . import sql_interface
from . import mass_spectrum
from . import settings


class Library_Enum(enum.Enum):
    NONE = enum.auto()
    NIST = enum.auto()
    CHEM_SPIDER = enum.auto()
    PUB_CHEM = enum.auto()
    USER = enum.auto()
    MANUAL_CONFIRMED = enum.auto()


class ID_Library(object):
    NORM_MAX = 1

    def __init__(self):
        self.type = Library_Enum.NONE
#         self.mass_plot = None  # Lazy init of output window

#     def plot_data(self, selected_sample, selected_peak, library_results):
#         if self.mass_plot is None:
#             self.mass_plot = plot.MassSpecPlot(None, title="Library Comparison Mass Spectrum")
#         self.mass_plot.show_peak_mass_spectrum(selected_peak, selected_sample, normalize=self.NORM_MAX, first_plot=True)
#         self.mass_plot.show_reference_spectrum(library_results, first_plot=False)

    def search(self, mass_spec: pyms.Spectrum.MassSpectrum, target_rt: float, max_hits: int=5):
        print(f"Search not implemented in ID_Library  {mass_spec}, {max_hits}")  # Make this class abstract?


class NIST(ID_Library):
    NORM_MAX = 999  # This appears to the maximum intensity in the reference data mass spec

    def __init__(self, nist_path: pathlib.Path, working_dir_path: pathlib.Path):
        super().__init__()
        self.type = Library_Enum.NIST
        # self.engine = pyms_nist_search.Engine(nist_path, pyms_nist_search.NISTMS_MAIN_LIB, working_dir_path)  # @UndefinedVariable

    def search(self, mass_spec: pyms.Spectrum.MassSpectrum, target_rt: float, max_hits: int =15) -> list:
        return []
        # return self.engine.full_search_with_ref_data(mass_spec, n_hits=max_hits)


class UserExcel(ID_Library):
    NORM_MAX = 100

    def __init__(self, source_spreadsheet: pathlib.Path):
        super().__init__()
        self.type = Library_Enum.USER
        self._rt_data = []  # collections.defaultdict(list)
        self._ms_data = []
        self.ms_min_correlation = 0.7  # TODO: Parameter
        self.max_rt_difference_s = 30  # TODO: Parameter
        self._source_path = source_spreadsheet
        self._column_names = {
            'compound': (1, 'name'), 'formula': (2, 'formula'), 'mw': (3, 'mw'), 'rt': (4, 'rt'),
            'cas#': (11, 'cas'), 'odor': (12, "odor"), 'synonyms': (13, "Synonyms"),
            'm/z 1': (5, 'mass1'), 'frac 1': (6, 'portion1'),
            'm/z 2': (7, 'mass2'), 'frac 2': (8, 'portion2'),
            'm/z 3': (9, 'mass3'), 'frac 3': (10, 'portion3')
        }
        if source_spreadsheet:  # and source_spreadsheet.is_file():   
            try:
                self.source_wb = openpyxl.load_workbook(source_spreadsheet)
                self.source_ws = self.source_wb.active
                self.active = True
            except FileNotFoundError:
                print(f"Could not load user identification excel {source_spreadsheet}")
                self.active = False
            except PermissionError:
                logging.warning(f"User identification excel {source_spreadsheet} is open or locked")
                self.active = False
            except Exception as _e:
                logging.exception(f"Could not load user identification excel {source_spreadsheet}")
                self.active = False
        else:
            self.active = False

        if self.active:
            self._read_excel()

    def _read_excel(self):
        first_line = next(self.source_ws.rows, None)
        if first_line is None:
            logging.error("Could not read excel user library first row")
        self._read_header(first_line)
        for row in self.source_ws.rows:
            new_entry = UserEntry()
            valid = new_entry.load(row, self._column_names)
            if valid:
                if new_entry.rt is not None:
                    self._rt_data.append(new_entry)
                elif new_entry.mass1 is not None:
                    self._ms_data.append((new_entry.mass1, new_entry))
        self._rt_data.sort(key=lambda x: x.rt if x.rt is not None else -numpy.Inf)
        logging.info(f"Read RT:{len(self._rt_data)} + Mass:{len(self._ms_data)} entries from user library {str(self._source_path)}")

    def _read_header(self, header_tuple: tuple):
        for index, entry in enumerate(header_tuple):  # TODO: Make sure defaults don't duplicate updated values
            if entry.value is not None and entry.value.casefold() in self._column_names:
                _old_index, atr_name = self._column_names[entry.value.casefold()]
                self._column_names[entry.value.casefold()] = (index, atr_name)

    def search(self, mass_spec: pyms.Spectrum.MassSpectrum, target_rt_s: float, max_hits: int=5) -> list:
        if not self.active:
            return None

        results = []
        if target_rt_s is not None:
            rt_list = [x.rt for x in self._rt_data]
            closest_rt_index = numpy.searchsorted(rt_list, target_rt_s, side='left')
            rt_indicies_up = self._add_rt(rt_list, target_rt_s, closest_rt_index, 1)
            rt_indicies_down = self._add_rt(rt_list, target_rt_s, closest_rt_index-1, -1)
            candidates = []
            candidates.extend(numpy.array(self._rt_data)[rt_indicies_up])
            candidates.extend(numpy.array(self._rt_data)[rt_indicies_down])
            if mass_spec is None:
                results.extend(candidates)
            else:
                ms_list = [x.mass_spec for x in candidates]
                score_list = numpy.array(mass_spectrum.MassSpectrum.compare_to_target_ms(mass_spec, ms_list, verify_highest=True))
                selected = score_list > self.ms_min_correlation
                has_no_ms = numpy.array(numpy.isnan([x.min_mass for x in ms_list]), dtype=numpy.bool_)
                results.extend(numpy.array(candidates)[selected | has_no_ms])

        if mass_spec is not None:
            intensity_sort = numpy.argsort(mass_spec.intensity_list)[::-1]
            base_i = mass_spec.intensity_list[intensity_sort[0]]
            mass_select = 1
            while True:
                next_i = mass_spec.intensity_list[intensity_sort[mass_select]]
                if next_i <= self.ms_min_correlation*base_i or mass_select >= len(intensity_sort)-1:
                    break
                mass_select += 1
            target_masses = numpy.array(mass_spec.mass_list)[intensity_sort][:mass_select]  # The largest ion must be included here to be considered
            short_list = [x[1] for x in self._ms_data if x[0] in target_masses]  # and x[1].rt is None]  # Ignore items that were scanned by rt already
            short_ms = [x.mass_spec for x in short_list]
            short_scores = mass_spectrum.MassSpectrum.compare_to_target_ms(mass_spec, short_ms, verify_highest=True)
            selected_ms = numpy.array(short_scores) > 0.7  # TODO: Parameter
            results.extend(numpy.array(short_list)[selected_ms])
        # TODO: Score, sort by priority, truncate list length to max_hits
        results = [(x.name, x) for x in results]  # Match nist result type
        return results

    def _add_rt(self, rt_list: list, target_rt_s: float, start_index: int, direction: int) -> list:
        result = []
        current_index = start_index
        while True:
            if current_index < 0 or current_index >= len(rt_list):
                break
            if abs(rt_list[current_index]-target_rt_s) > self.max_rt_difference_s:
                break
            result.append(current_index)
            current_index += direction
        return result


class UserEntry(sql_interface.Base):
    __tablename__ = "User_ID_Compound"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    formula = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    mw = sqlalchemy.Column(sqlalchemy.Integer, nullable=True)
    cas = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    flavor_list = sqlalchemy.Column(sqlalchemy.String, nullable=True)

    rt = sqlalchemy.Column(sqlalchemy.Float, nullable=True)
    mass1 = sqlalchemy.Column(sqlalchemy.Integer, nullable=True)
    portion1 = sqlalchemy.Column(sqlalchemy.Integer, nullable=True)
    mass2 = sqlalchemy.Column(sqlalchemy.Integer, nullable=True)
    portion2 = sqlalchemy.Column(sqlalchemy.Integer, nullable=True)
    mass3 = sqlalchemy.Column(sqlalchemy.Integer, nullable=True)
    portion3 = sqlalchemy.Column(sqlalchemy.Integer, nullable=True)

    def load(self, row: tuple, column_lookup: dict) -> bool:
        try:
            for _column_name, (column_index, attribute_name) in column_lookup.items():
                setattr(self, attribute_name, row[column_index].value)
        except IndexError:
            return False  # The excel sheet does not have as many columns as the defaults expect and a header is missing
        try:
            self.rt = float(self.rt)*60.0
        except (ValueError, TypeError) as _e:
            no_rt = True
        else:
            no_rt = False
        try:
            int(self.mass1)
        except (ValueError, TypeError) as _e:
            if no_rt:
                return False  # Must have one of rt or mass1

        return True

    @property
    def mass_spec(self):
        mass_list = []
        portion_list = []
        if self.mass1 and self.portion1 and not numpy.isnan(self.mass1):
            mass_list.append(self.mass1)
            portion_list.append(self.portion1)
        if self.mass2 and self.portion2 and not numpy.isnan(self.mass2):
            mass_list.append(self.mass2)
            portion_list.append(self.portion2)
        if self.mass3 and self.portion3 and not numpy.isnan(self.mass3):
            mass_list.append(self.mass3)
            portion_list.append(self.portion3)
        return mass_spectrum.MassSpectrum.convert_from_lists(mass_list, portion_list)

# See sql_interface
# class Match(object):
#     ''' Estimated correlation between a compound and a peak '''
#     def __init__(self):
#         self.compound = None
#         self.peak = None
#         self.score = 0.0
# 
# 
# class Compound(object):
#     '''
#     A compound used for matching during identification.
#     '''
#     def __init__(self):
#         self.name = "N/A"
#         self.fragment_masses = collections.OrderedDict()
#         self.rt_kovat = -1.0
#         self.cas = "N/A"
#         self.reference_data = None
#         self._source_library_type = Library_Enum.NONE
