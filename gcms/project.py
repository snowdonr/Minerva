'''
Created on Sep. 16, 2021

@author: Ryan
'''
from __future__ import annotations
import logging
import pathlib
import wx
import sqlalchemy.orm
import sqlalchemy.dialects.mysql

from . import sample
from . import analysis
from . import plot
from . import identification
from . import sql_interface
from . import settings
from . import main_ui


class Project(sql_interface.Base):
    '''
    Holding class for interface selections and data structures
    '''
    __tablename__ = "Minerva_Project"

    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    experiment_name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    pr = sqlalchemy.Column(sqlalchemy.String, nullable=True)  # TODO: Refer to table of people
    process_version = sqlalchemy.Column(sqlalchemy.String, nullable=True)
    text_type = sqlalchemy.Text().with_variant(sqlalchemy.dialects.mysql.LONGTEXT(), "mysql")
    comments = sqlalchemy.Column(text_type, nullable=True)
    create = sqlalchemy.Column(sqlalchemy.DateTime, default=sqlalchemy.func.now())
    update = sqlalchemy.Column(sqlalchemy.DateTime, onupdate=sqlalchemy.func.utc_timestamp())

    aligned_peaks = sqlalchemy.orm.relationship("AlignedPeakSet_SQL", back_populates="project")
    analysis_settings = sqlalchemy.orm.relationship("AnalysisSettings_SQL", back_populates="project")
    samples = sqlalchemy.orm.relationship("Sample", back_populates="project")

    def __init__(self, sql_input: sql_interface.SQL_Interface, config_input: settings.Setting, version: str, main_ref: main_ui.MainUI):
        self.mainUI = main_ref
        self.sample_set = sample.SampleSet()
        self.sql = sql_input
        self.config = config_input
        self.process_version = version
        self.user_library = identification.UserExcel(None)  # .active will be False
        self.nist_library = None
        self._analysis_set = None
        self.chrom_plot = None

    def startup(self):
        if self.config.output_sqlite:
            sql_path = pathlib.Path(self.config.analysis_directory) / (self.config.experiment_name+"_Output.db")
            try:
                self.sql.link_sqlite(sql_path, clear_file=not sql_path.is_file())
            except sqlalchemy.exc.OperationalError as _e:
                logging.exception(f"Failed to open db file {sql_path}")
        else:
            self.sql.link_cloud_postgre(self)
        try:
            self.nist_library = identification.NIST(self.config.main_nist_path, self.config.nist_working_dir_path)
        except Exception as _e:
            logging.warning(f"Failed to start NIST libarary at {self.config.main_nist_path} {self.config.nist_working_dir_path}")

    def load_samples(self, input_dir: pathlib.Path=None, progress_bar: wx.Gauge=None):
        self.sample_set.clear()
        if input_dir is None:
            input_dir = pathlib.Path(self.config.data_path)
        if input_dir and input_dir.is_dir():  # Override setting if an argument is provided
            success = self.sample_set.load_directory(self.config, progress_bar)
        else:
            logging.error("Input location is not a directory")
            success = False
        if success:
            for samp in self.sample_set:
                self.samples.append(samp)
        else:
            logging.error("ERROR: No files found")  # Program is expected to crash, but is allowed to run for debugging
        self._update_duplicate_sample_names()

    def _update_duplicate_sample_names(self):
        sample_taken_names = set()
        for current_sample in self.sample_set:
            if current_sample.sample_number in sample_taken_names:
                next_name = str(current_sample.sample_number) + "_1"
                i = 1
                while next_name in sample_taken_names:
                    i += 1
                    next_name = next_name[:-1]+str(i)
                current_sample.sample_number = next_name
            sample_taken_names.add(current_sample.sample_number)

        # test_sample_set[0].get_mass_sums(50, 550)

    def merged_sample_plot(self):
        pass

    def identify_compounds(self):
        test_sample_data = self.sample_set[0]
        peaks = test_sample_data.get_filtered_peaks(lambda x: x.area>30000)

        identification_results = []
        for current_peak in peaks:
            test_result = self.nist_library.search(current_peak.mass_spectrum)
            identification_results.append((current_peak, test_result))

        return identification_results

    def show_debug_plots(self):
        if self.chrom_plot is None:
            self.chrom_plot = plot.SingleChromatogramPlot(None)
            self.mass_spec_plot = plot.MassSpecPlot(self.chrom_plot)
            self.chrom_plot._mass_plot = self.mass_spec_plot
            self.debug_plot = plot.LoadedDataPlot(self.mass_spec_plot, self.sample_set)
            self.debug_plot.change_plot(0, maintain_range=False)
        else:
            self.chrom_plot.Show()
            self.mass_spec_plot.Show()
            self.debug_plot.Show()

    def process_samples(self) -> bool:
        ''' Create analysis set and return status '''
        self.name = self.config.experiment_name
        self.experiment_name = self.config.method_analysis_name
        self.user_library = identification.UserExcel(self.config.user_library_location)
        sample_names = [x.sample_number for x in self.sample_set.data]
        if len(set(sample_names)) != len(sample_names):
            print("Duplicated sample name")
        try:
            # test_sample_set.remove_zero_peak_samples()  Not here: Tends to be an issue AFTER mass cropping
            self._analysis_set = analysis.Analysis_Set(self.sample_set, self.config.experiment_name, self.config, self.sql)  # [13490, 15196, 15205]
        except Exception as _e:
            logging.exception("Analysis")
            return False
        if self.config.merge_aligned_peaks:
            self._analysis_set.merge_aligned_peaks()
            try:
                self._analysis_set.write_aligned_csv(self._analysis_set.alignment_data_groups, pathlib.Path(self.config.analysis_directory) / (self.config.experiment_name+"_grouping.csv"))
            except PermissionError:
                print("Analysis Relationship could not write to 'grouping' Excel file")
            try:
                self._analysis_set.write_aligned_csv(self._analysis_set.merged_data, pathlib.Path(self.config.analysis_directory) / (self.config.experiment_name+"_grouped.csv"))
            except PermissionError:
                print("Analysis Relationship could not write to 'grouped' Excel file")

        self.mainUI.full_progress_gauge.SetValue(main_ui.DONE_ALIGN)
        self._analysis_set.add_library_matches(self.user_library, self.sql)
        self._analysis_set.add_library_matches(self.nist_library, self.sql)

        if self.config.merge_aligned_peaks:
            output_data = self._analysis_set.merged_data
        else:
            output_data = self._analysis_set.alignment_data_groups
        try:
            self._analysis_set.write_aligned_csv(output_data, pathlib.Path(self.config.analysis_directory) / (self.config.experiment_name+"_results.csv"), output_compounds=True)
        except PermissionError:
            print("Analysis Relationship could not write to 'results' Excel file")

        try:
            self.sql.update_project_table(self, self._analysis_set)
        except sqlalchemy.exc.OperationalError as _e:
            logging.exception(f"Could not store data")
        return True

    def get_aligned_set(self) -> list:
        # result = self._analysis_set._alignment_data.peaks_aligned  # TODO: use alignment_data_groups
        # TODO: Filter when samples are eliminated
        return self._analysis_set.alignment_data_groups

    def __getstate__(self):
        return None  # used to locate pickle attempts on this object
