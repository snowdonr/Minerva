'''

Created on Sep. 21, 2021

@author: Ryan
'''
import enum
import logging
import pathlib
import time
import json
import openpyxl
import numpy
import sqlalchemy.orm
import urllib.request
import urllib.parse

import pyms.Spectrum
from . import sql_base
from . import mass_spectrum
from .mass_spectrum import MassSpectrum_SQL


class Library_Enum(enum.Enum):
    NONE = enum.auto()
    NIST = enum.auto()
    CHEM_SPIDER = enum.auto()
    PUB_CHEM = enum.auto()
    USER = enum.auto()
    MANUAL_CONFIRMED = enum.auto()


class ID_Library(object):
    ''' Base class for peak-> compound identifications '''
    NORM_MAX = 1

    def __init__(self):
        self.type = Library_Enum.NONE

    def search(self, mass_spec: pyms.Spectrum.MassSpectrum, target_rt: float, max_hits: int=5):
        print(f"Search not implemented in ID_Library  {mass_spec}, {max_hits}, {target_rt}")  # Make this class abstract?


class NIST(ID_Library):
    ''' Disabled reference to pyms_nist_search library '''
    NORM_MAX = 999  # This appears to the maximum intensity in the reference data mass spec

    def __init__(self, nist_path: pathlib.Path, working_dir_path: pathlib.Path):
        super().__init__()
        self.type = Library_Enum.NIST
        # self.engine = pyms_nist_search.Engine(nist_path, pyms_nist_search.NISTMS_MAIN_LIB, working_dir_path)  # @UndefinedVariable

    def search(self, mass_spec: pyms.Spectrum.MassSpectrum, target_rt: float, max_hits: int =15) -> list:
        return []
        # return self.engine.full_search_with_ref_data(mass_spec, n_hits=max_hits)


class UserExcel(ID_Library):
    ''' User provided library (spreadsheet) for matching peaks to compounds '''
    NORM_MAX = 100

    def __init__(self, source_spreadsheet: pathlib.Path):
        super().__init__()
        self.type = Library_Enum.USER
        self._rt_data = []  # collections.defaultdict(list)
        self._ms_data = []
        self.ms_min_correlation = 0.7  # TODO: Parameter
        self.max_rt_difference_s = 30  # TODO: Parameter
        self._source_path = source_spreadsheet
        self._column_names = {  # Note: Column name keys should be all lower case to match casefold
            'category': (4, 'category'), 'compound': (5, 'name'), 'nickname': (6, 'nickname'), 'quality': (7, 'quality'),
            'formula': (8, 'formula'), 'mw': (9, 'mw'), 'rt': (10, 'rt'),
            'cas#': (11, 'cas'), 'cid': (17, "cid"), "inchi": (18, "inchi"), "inchikey": (19, "inchikey"),
            'iupac name': (12, 'iupacname'),
            'monoisotopic mass': (13, 'monoisotopicmass'), 'xlogp': (14, 'xlogp'),
            'canonical smiles': (15, 'canonicalsmiles'),
            'odor': (38, "odor"), 'synonyms': (39, "synonyms"),
            'm/z 1': (22, 'mass1'), 'frac 1': (23, 'portion1'),
            'm/z 2': (24, 'mass2'), 'frac 2': (25, 'portion2'),
            'm/z 3': (25, 'mass3'), 'frac 3': (27, 'portion3')
        }
        if source_spreadsheet:
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
        do_write=False
        for row in self.source_ws.rows:
            # Check on header column position
            column_position, _category_name = self._column_names['m/z 1']
            if row[column_position].value == 'm/z 1':
                continue
            new_entry = UserEntry()
            entry_status = new_entry.load(row, self._column_names)
            if entry_status != "error":
                if new_entry.rt is not None:
                    self._rt_data.append(new_entry)
                elif new_entry.mass1 is not None:
                    self._ms_data.append((new_entry.mass1, new_entry))
                if entry_status == "update" or entry_status == "write":
                    do_write = True
                elif entry_status[:7] == "ignore/":
                    match_type = entry_status[7:]
                    output_column = None
                    for _type_name, (column_number, category_name) in self._column_names.items():
                        if category_name == match_type:
                            output_column = column_number
                            break
                    if output_column:
                        row[output_column].value = "(?)"+str(row[output_column].value)
                        do_write = True
        if do_write:
            try:
                self.source_wb.save(self._source_path)
            except IOError:
                logging.error(f"Could not save user library {self._source_path}")
        self._rt_data.sort(key=lambda x: x.rt if x.rt is not None else -numpy.Inf)
        logging.info(f"Read RT:{len(self._rt_data)} + Mass:{len(self._ms_data)} entries from user library {str(self._source_path)}")

    def _read_header(self, header_tuple: tuple):
        for index, entry in enumerate(header_tuple):  # TODO: Make sure defaults don't duplicate updated values
            if entry.value is not None and entry.value.casefold() in self._column_names:
                _old_index, atr_name = self._column_names[entry.value.casefold()]
                self._column_names[entry.value.casefold()] = (index, atr_name)

    def search(self, mass_spec: pyms.Spectrum.MassSpectrum, target_rt_s: float, max_hits: int=5) -> list:
        '''
        Search for a library match from peak data
        '''
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
        results = [(x.name, x) for x in results]  # Match nist library result type
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


class UserEntry(sql_base.Base):
    ''' An SQL entry for compound information '''
    __tablename__ = "User_ID_Compound"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    nickname = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    category = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    formula = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    mw = sqlalchemy.Column(sqlalchemy.Float, nullable=True)
    iupacname = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    inchikey = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    inchi = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    cid = sqlalchemy.Column(sqlalchemy.Integer, nullable=True)
    cas = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    canonicalsmiles = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    molecularformula = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    # molecularweight = sqlalchemy.Column(sqlalchemy.String, nullable=True)  # now mw to match NIST format
    # exactmass = sqlalchemy.Column(sqlalchemy.String, nullable=True)  # see monoisotpic mass
    monoisotopicmass = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    xlogp = sqlalchemy.Column(sqlalchemy.String, nullable=True)

    flavor_list = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    rt = sqlalchemy.Column(sqlalchemy.Float, nullable=True)
    spectrum_id = sqlalchemy.Column(sqlalchemy.Integer, sqlalchemy.ForeignKey('Mass_Spectrum.id'))
    spectrum = sqlalchemy.orm.relationship("MassSpectrum_SQL", back_populates="user_peak")
    aligned_match = sqlalchemy.orm.relationship("CompoundMatch_SQL", back_populates="compound")

    def __init__(self):
        super().__init__()
        if self.spectrum is None:
            self.spectrum = MassSpectrum_SQL()

    def load(self, row: tuple, column_lookup: dict) -> str:
        try:
            for _column_name, (column_index, attribute_name) in column_lookup.items():
                setattr(self, attribute_name, row[column_index].value)
        except IndexError:
            return "error"  # The excel sheet does not have as many columns as the defaults expect and a header is missing
        requires_write = "exists"
        if not self.inchi or not self.inchikey:
            try:
                requires_write = self._check_and_update(row, column_lookup)
            except Exception as _e:
                logging.exception("Excel update failed")
        try:
            self.rt = float(self.rt)*60.0  # Assumes user rt are in minutes
        except (ValueError, TypeError) as _e:
            self.rt = None
            no_rt = True
        else:
            no_rt = False
        try:
            int(self.mass1)
        except (ValueError, TypeError) as _e:
            self.mass1 = None
            if no_rt:
                logging.info(f"No identification in User Library at {row}")  # Must have one of rt or mass1
        return requires_write

    @property
    def molecularweight(self):
        return self.mw

    @molecularweight.setter
    def molecularweight(self, value):
        self.mw = value

    @property
    def mass1(self):
        return self.spectrum.values[0].mass

    @mass1.setter
    def mass1(self, value):
        if len(self.spectrum.values) < 1:
            self.spectrum.values.append(mass_spectrum.MassSpectrumEntry())
        self.spectrum.values[0].mass = value

    @property
    def portion1(self):
        result = self.spectrum.values[0].area
        if not result:
            result = self.spectrum.values[0].intensity  # Could cause issues if area/intensity is not consistently set
        return result

    @portion1.setter
    def portion1(self, value):
        if len(self.spectrum.values) < 1:
            self.spectrum.values.append(mass_spectrum.MassSpectrumEntry())
        self.spectrum.values[0].area = value

    @property
    def mass2(self):
        return self.spectrum.values[1].mass

    @mass2.setter
    def mass2(self, value):
        if len(self.spectrum.values) < 2:
            self.spectrum.values.append(mass_spectrum.MassSpectrumEntry())
        self.spectrum.values[1].mass = value

    @property
    def portion2(self):
        result = self.spectrum.values[1].area
        if not result:
            result = self.spectrum.values[1].intensity
        return result

    @portion2.setter
    def portion2(self, value):
        if len(self.spectrum.values) < 2:
            self.spectrum.values.append(mass_spectrum.MassSpectrumEntry())
        self.spectrum.values[1].area = value

    @property
    def mass3(self):
        return self.spectrum.values[2].mass

    @mass3.setter
    def mass3(self, value):
        if len(self.spectrum.values) < 3:
            self.spectrum.values.append(mass_spectrum.MassSpectrumEntry())
        self.spectrum.values[2].mass = value

    @property
    def portion3(self):
        result = self.spectrum.values[2].area
        if not result:
            result = self.spectrum.values[2].intensity
        return result

    @portion3.setter
    def portion3(self, value):
        if len(self.spectrum.values) < 3:
            self.spectrum.values.append(mass_spectrum.MassSpectrumEntry())
        self.spectrum.values[2].area = value

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

    def _check_and_update(self, row, column_lookup) -> str:
        update_excel = "exists"
        if self.inchi:
            id_type = "inchi"
            id_details = self.inchi
        elif self.inchikey:
            id_type = "inchikey"
            id_details = self.inchikey
        elif self.cid:
            id_type = "cid"
            id_details = self.cid
        elif self.iupacname:
            id_type = "name"
            id_details = self.iupacname
        elif self.name:
            id_type = "name"
            id_details = self.name
        else:
            try:
                row_number = row[0].row
            except Exception as _e:
                row_number = "?"
            logging.info(f"Unknown use row at {row_number}")
            return "error"
        if isinstance(id_details, str) and id_details.startswith("(?)"):
            logging.info(f"Skipping {id_details}")
            return "skip"
        lookup_results = self.scan_pubchem(id_type, id_details)
        if len(lookup_results) > 1:
            logging.info(f"Multiple results for {id_details}")
        elif len(lookup_results) < 1:
            return f"ignore/{id_type}"
        selected_results = lookup_results[0]

        for check_property in self.property_list:
            if check_property in selected_results and selected_results[check_property] is not None:
                if getattr(self, check_property.casefold()) and getattr(self, check_property.casefold()) != selected_results[check_property]:
                    logging.warning(f"Updating property {check_property} {getattr(self, check_property.casefold())} -> {selected_results[check_property]}")
                setattr(self, check_property.casefold(), selected_results[check_property])
                update_excel = "update"
        self.cid = selected_results["CID"]
        # update row
        try:
            for _column_name, (column_index, attribute_name) in column_lookup.items():
                if getattr(self, attribute_name) is None:
                    continue
                if attribute_name == "rt":
                    new_value = getattr(self, attribute_name)/60.0  # Userlibrary uses minutes
                else:
                    new_value = getattr(self, attribute_name)
                row[column_index].value = new_value
        except Exception as _e:
            logging.exception(f"Update of excel row {row} failed.")
        return update_excel

    property_list = "IUPACName", "InChI", "InChIKey", "MolecularFormula", "MolecularWeight", "MonoisotopicMass", "XLogP", "CanonicalSMILES"  # "cid" always included

    @classmethod
    def scan_pubchem(cls, id_type: str, id_details: str):
        '''
        Given a type as one of "name", "inchi", "inchikey", "cid" do a lookup of compound data.
        Applies a delay to reduce requests per second if this function is run across a list of targets
        type cid allows a comma separated list of details for multiple results
        '''
        property_str = ",".join(cls.property_list)
        target_string = urllib.parse.quote(str(id_details).strip().replace("/", "."))
        if id_type.casefold() == "inchi":
            try:
                inchi_target = f"inchi={str(id_details).strip()}".encode('utf-8')
                pubchem_response = urllib.request.urlopen(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{id_type}/property/{property_str}/JSON", data=inchi_target).read()
            except Exception as _e:
                logging.exception("Encode failed")
                return []
        else:
            try:  # with urllib... as pubchem_respone:
                pubchem_response = urllib.request.urlopen(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{id_type}/{target_string}/property/{property_str}/JSON").read()
            except urllib.error.HTTPError:
                logging.error(f"User library query failed on {id_type}: {str(id_details)}")
                return []
            except urllib.error.URLError:
                logging.error(f"User library url failed on {id_type}: {str(id_details)}")
                return []
        # xml propertytable, properties CID InChI InChIKey
        json_result = json.loads(pubchem_response)
        time.sleep(0.25)  # time limit for below 5 requests/sec
        return json_result["PropertyTable"]["Properties"]
