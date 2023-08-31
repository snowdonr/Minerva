'''
Created on Sept 16, 2021

@author: Ryan Snowdon
'''
import matplotlib
matplotlib.interactive(True)  # show non-blocking
matplotlib.use('WXAgg')
import wx
from wx.lib.delayedresult import startWorker  # @UnresolvedImport
import pathlib
import logging
import functools
import itertools
import typing
import datetime
import sys

from . import project
from . import plot
from . import sql_interface
from . import settings

PROGRESS_MAX = 10
DONE_LOAD = 1
DONE_ALIGN = 2
DONE_COMPOUND_LOOKUP = 3
FULL_PROGRESS_ITEMS = 3


class MainProcessing(object):
    '''
    Non-interactive processing using configuration file and producing a full set of results
    plus caches for plotting using the UI version
    '''
    def __init__(self, version: str, config_file: typing.Optional[pathlib.Path], data_input_dir: typing.Optional[pathlib.Path]):
        sql = sql_interface.SQL_Interface()
        if config_file is not None:
            config = settings.Setting(config_file)
        else:
            config = settings.Setting(pathlib.Path.home() / "Minerva.ini")
        self.config = config
        if data_input_dir is not None:
            self.config.data_path = data_input_dir
        self._project = project.Project(sql, config, version, self)

        # self.config.load_align_pickle = False  # But does not save the change to the ini file
        logging.info(f"Loading {str(datetime.datetime.now().time())}")
        self._project.startup()
        self._project.load_samples(progress_bar=None)
        self._project.process_samples()
        logging.info(f"Finished alignment")


class MainUI(wx.Frame):
    '''
    UI related functions, resulting data is stored in _project
    '''
    def __init__(self, version: str, default_input_dir: pathlib.Path):
        super().__init__(None, wx.ID_TOP, f"Minerva {version}", size=(490, 380))
        # Setup main objects
        self.version = version
        config = settings.Setting(pathlib.Path.home() / "Minerva.ini")
        self.sql = sql_interface.SQL_Interface()
        self.config = config
        self._project = None  # set during start_processing

        self._status_loading = False  # Not set atomically, but help to indicate what operations make sense
        self._status_loaded = False
        self._status_aligning = False
        self._status_aligned = False

        self._merged_sample_plot = None

        if default_input_dir is not None:
            self.config.data_path = default_input_dir

        main_sizer = wx.GridBagSizer(vgap=1, hgap=1)
        self.SetSizer(main_sizer)

        tab_notebook = wx.Notebook(self, id=wx.ID_ANY, style=wx.BK_DEFAULT)
        experiment_panel = wx.Panel(tab_notebook)
        experiment_sizer = wx.GridBagSizer(vgap=1, hgap=1)
        experiment_panel.SetSizer(experiment_sizer)
        tab_notebook.AddPage(experiment_panel, "Experiment")
        method_panel = wx.Panel(tab_notebook)
        method_sizer = wx.GridBagSizer(vgap=1, hgap=1)
        method_panel.SetSizer(method_sizer)
        tab_notebook.AddPage(method_panel, "Method")
        directory_panel = wx.Panel(tab_notebook)
        directory_sizer = wx.GridBagSizer(vgap=1, hgap=1)
        directory_panel.SetSizer(directory_sizer)
        tab_notebook.AddPage(directory_panel, "File Locations")

        grid_row = itertools.count(0)
        main_sizer.Add(tab_notebook, pos=(next(grid_row), 0), span=(1, 3), flag=wx.EXPAND)
        self._ui_section_list = []

        experiment_row = itertools.count(0)
        self._ui_section_list.append(self._create_configuration_text(experiment_panel, "Experiment Name", "experiment_name", str, next(experiment_row)))
        self._ui_section_list.append(self._create_target_path(experiment_panel, "Source Directory", "data_path", next(experiment_row), select_directory=True))
        # self._ui_section_list.append(self._create_target_path(experiment_panel, "Sample Metadata", "sample_relation_path", next(experiment_row)))
        self._ui_section_list.append(self._create_target_path(experiment_panel, "Destination Directory", "analysis_directory", next(experiment_row), select_directory=True))

        directory_row = itertools.count(0)
        # self._ui_section_list.append(self._create_target_path(directory_panel, "Nist Library Directory", "main_nist_path", next(directory_row), select_directory=True))
        # self._ui_section_list.append(self._create_target_path(directory_panel, "Nist Working Dir", "nist_working_dir_path", next(directory_row), select_directory=True))
        self._ui_section_list.append(self._create_target_path(directory_panel, "Excel Library", "user_library_location", next(directory_row)))
        self._ui_section_list.append(self._create_target_path(directory_panel, "Temporary Directory", "temp_path", next(directory_row), select_directory=True))

        method_row = itertools.count(0)
        self._ui_section_list.append(self._create_configuration_text(method_panel, "Method Name", "method_analysis_name", str, next(method_row)))

        self._lo_rt_limit_box = self._create_configuration_text(method_panel, "Minimum RT", "lo_rt_limit", str, next(method_row))
        self._lo_rt_limit_box.SetToolTip("Must end with unit: s or m")
        self._ui_section_list.append(self._lo_rt_limit_box)
        self._hi_rt_limit_box = self._create_configuration_text(method_panel, "Maximum RT", "hi_rt_limit", str, next(method_row))
        self._hi_rt_limit_box.SetToolTip("Must end with unit: s or m")
        self._ui_section_list.append(self._hi_rt_limit_box)
        self._ui_section_list.append(self._create_configuration_text(method_panel, "Minimum Area", "filter_minimum_area", float, next(method_row)))

        self._ui_section_list.append(self._create_configuration_text(method_panel, "RT sensitivity (sec)", "rt_sensitivity_s", float, next(method_row)))
        self._gap_penalty_box = self._create_configuration_text(method_panel, "Gap Penalty", "gap_penalty", float, next(method_row))
        self._gap_penalty_box.SetToolTip("[0-1] Lower values favors insertion of gaps, higher values align less related peaks")
        self._ui_section_list.append(self._gap_penalty_box)

        self._ui_section_list.append(self._create_configuration_text(method_panel, "Min Peak % (Normalized)", "peak_trim_percent_of_max", float, next(method_row)))
        self._ui_section_list.append(self._create_configuration_text(method_panel, "Minimum Ion count", "min_ions", int, next(method_row)))

        self._ignore_cache = wx.CheckBox(self, wx.ID_ANY, "Recalculate Alignment Cache")  # State is not saved, default off
        main_sizer.Add(self._ignore_cache, pos=(next(grid_row), 0))

        progress_row = next(grid_row)
        label = wx.StaticText(self, wx.ID_ANY, "Item Progress")
        main_sizer.Add(label, pos=(progress_row, 0), flag=wx.ALIGN_LEFT)
        self.item_progress_gauge = wx.Gauge(self, id=wx.ID_ANY, style=wx.GA_HORIZONTAL, name="item_progress_gauge")
        self.item_progress_gauge.SetRange(PROGRESS_MAX)
        self.item_progress_gauge.SetValue(0)
        main_sizer.Add(self.item_progress_gauge, pos=(progress_row, 1), span=(1, 2), flag=wx.EXPAND)

        progress_row = next(grid_row)
        label = wx.StaticText(self, wx.ID_ANY, "Overall Progress")
        main_sizer.Add(label, pos=(progress_row, 0), flag=wx.ALIGN_LEFT)
        self.full_progress_gauge = wx.Gauge(self, id=wx.ID_ANY, style=wx.GA_HORIZONTAL, name="full_progress_gauge")
        self.full_progress_gauge.SetRange(FULL_PROGRESS_ITEMS)
        self.full_progress_gauge.SetValue(0)
        main_sizer.Add(self.full_progress_gauge, pos=(progress_row, 1), span=(1, 2), flag=wx.EXPAND)

        self._load_button = wx.Button(self, wx.ID_ANY, "Load and Align")
        self.Bind(wx.EVT_BUTTON, self.start_processing, self._load_button)
        main_sizer.Add(self._load_button, pos=(next(grid_row), 0), span=(1, 2), flag=wx.ALIGN_LEFT)

        self._raw_plots_button = wx.Button(self, wx.ID_ANY, "Raw Data Plots")
        self.Bind(wx.EVT_BUTTON, self.show_raw_plots, self._raw_plots_button)
        main_sizer.Add(self._raw_plots_button, pos=(next(grid_row), 0), span=(1, 2), flag=wx.ALIGN_LEFT)
        self._raw_plots_button.Enable(enable=False)

        self._plot_aligned_button = wx.Button(self, wx.ID_ANY, "Plot Aligned Data")
        self.Bind(wx.EVT_BUTTON, self.merged_sample_plot, self._plot_aligned_button)
        main_sizer.Add(self._plot_aligned_button, pos=(next(grid_row), 0), span=(1, 2), flag=wx.ALIGN_LEFT)
        self._plot_aligned_button.Enable(enable=False)

        self.Bind(wx.EVT_CLOSE, self.on_close)
        main_sizer.AddGrowableCol(1, proportion=1)
        experiment_sizer.AddGrowableCol(1, proportion=1)
        method_sizer.AddGrowableCol(1, proportion=1)
        directory_sizer.AddGrowableCol(1, proportion=1)
        self.Layout()
        self.Show()
        if default_input_dir is not None:
            self.start_processing(None)

    def _create_target_path(self, parent: wx.Window, label_text: str, settings_item: str, current_row: int, select_directory: bool=False):
        ''' Create UI elements to set a path and save it to the settings '''
        parent_sizer = parent.GetSizer()
        label = wx.StaticText(parent, wx.ID_ANY, label_text)
        default_value = getattr(self.config, settings_item)
        path_box = wx.TextCtrl(parent, value=str(default_value), size=(160, -1), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.path_changed, path_box)

        selectButton = wx.Button(parent, wx.ID_ANY, "...")
        if select_directory:
            self.Bind(wx.EVT_BUTTON, functools.partial(self.update_directory_path, text_ctrl=path_box), selectButton)
        else:
            self.Bind(wx.EVT_BUTTON, functools.partial(self.update_path, text_ctrl=path_box), selectButton)
        parent_sizer.Add(label, pos=(current_row, 0), span=wx.DefaultSpan, flag=wx.ALIGN_CENTER_VERTICAL)
        parent_sizer.Add(path_box, pos=(current_row, 1), span=wx.DefaultSpan, flag=wx.EXPAND)
        parent_sizer.Add(selectButton, pos=(current_row, 2), span=wx.DefaultSpan, flag=wx.ALIGN_CENTER)
        path_box.type_caster = str
        path_box.settings_item = settings_item  # TODO: Inherit a new TextCtrl class for this?
        return path_box

    def _create_configuration_text(self, parent: wx.Window, label_text: str, settings_item: str, value_type_cast: typing.Callable,
                                   current_row: int) -> wx.TextCtrl:
        parent_sizer = parent.GetSizer()
        label = wx.StaticText(parent, wx.ID_ANY, label_text)
        default_value = getattr(self.config, settings_item)
        item_box = wx.TextCtrl(parent, value=str(default_value), size=(160, -1), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, functools.partial(self.item_changed, settings_item, value_type_cast), item_box)

        parent_sizer.Add(label, pos=(current_row, 0), span=wx.DefaultSpan, flag=wx.ALIGN_CENTER_VERTICAL)
        parent_sizer.Add(item_box, pos=(current_row, 1), span=wx.DefaultSpan, flag=wx.EXPAND)
        item_box.type_caster = value_type_cast
        item_box.settings_item = settings_item
        return item_box

    def path_changed(self, event: wx.Event):
        pass

    def item_changed(self, settings_name: str, value_type_caster: typing.Callable, source_event: wx.Event):
        ''' Save setting to configuration file '''
        new_value = source_event.GetEventObject().Value
        self.config.save_value(settings_name, value_type_caster(new_value))  # Note: Not committed!

    def update_path(self, _event, text_ctrl: wx.TextCtrl):
        ''' Update a text box value with the result of a user file selection dialog '''
        file_menu = wx.FileDialog(self, message="Choose path", wildcard="*", style=wx.FD_OPEN)  # , style=)
        if file_menu.ShowModal() == wx.ID_OK:
            path_value = pathlib.Path(file_menu.Directory)/ pathlib.Path(file_menu.Filename)
            text_ctrl.SetValue(str(path_value))

    def update_directory_path(self, _event, text_ctrl: wx.TextCtrl):
        ''' Update a text box value with the result of a user directory selection dialog '''
        file_menu = wx.DirDialog(self, message="Choose path")
        if file_menu.ShowModal() == wx.ID_OK:
            path_value = pathlib.Path(file_menu.Path)
            if path_value.is_dir():
                text_ctrl.SetValue(str(path_value))
            else:
                print(f"Could not find directory at {file_menu.Path}")

    def start_processing(self, _event):
        ''' Update UI state and begin processing samples '''
        # This is assumed to be fairly atomic prior to activating the startWorker, so UI bugs can be triggered here depending on timing
        self._load_button.Enable(enable=False)  # Not fully thread-safe
        self._project = project.Project(self.sql, self.config, self.version, self)

        for item in self._ui_section_list:
            self._check_and_update(item)

        if self._ignore_cache.IsChecked():
            self.config.load_align_pickle = False  # ** Deliberately does not save this change to the ini file!

        self.config.commit()

        self._status_loading = True
        self._status_aligned = False
        startWorker(consumer=self._finished_loading, workerFn=self._load_file_worker, jobID=f"loadf_{str(datetime.datetime.now().time())}")  # wargs=

    def _load_file_worker(self):
        ''' Threaded load operation '''
        self._project.startup()
        self._project.load_samples(progress_bar=self.item_progress_gauge)
        return True

    def _finished_loading(self, delayed_result):
        ''' End threaded load operation '''
        self._status_loading = False
        # self._load_button.Enable(enable=True)  # TODO: Fix sqlalchemy errors on second load
        print(f"Finished loading {delayed_result.getJobID()}")
        if delayed_result.get():
            self.full_progress_gauge.SetValue(DONE_LOAD)
            self._status_loaded = True
            self._raw_plots_button.Enable(enable=True)
            self.start_alignment_worker()

    def show_raw_plots(self, _event):
        ''' Show plots for the integrated but not aligned data (one sample visible at a time) '''
        if self._status_loading or not self._status_loaded:
            print("DEBUG: Can't identify until loading is complete")
            return  # Button should be disabled
        self._project.show_debug_plots()

    def merged_sample_plot(self, _event):
        ''' Show plots for the data post alignment '''
        if self._merged_sample_plot is None:
            self._merged_sample_plot = plot.SampleSetPlot(self._project)
        self._merged_sample_plot.Show()
        aligned_set = self._project.get_aligned_set()
        self._merged_sample_plot.show_plot(aligned_set)

    def start_alignment_worker(self):
        ''' Trigger a thread to process the alignment of the selected sample set '''
        startWorker(consumer=self._finished_alignment, workerFn=self._project.process_samples, jobID=f"align_{len(self._project.sample_set)}")

    def _finished_alignment(self, delayed_result):
        ''' Clean up after threaded alignment '''
        self._status_aligning = False
        print(f"Finished alignment {delayed_result.getJobID()}")
        if delayed_result.get():
            self.full_progress_gauge.SetValue(DONE_COMPOUND_LOOKUP)
            self._status_aligned = True
            self._plot_aligned_button.Enable(enable=True)

    def _check_and_update(self, ui_text_box: wx.TextCtrl):
        config_item_name = ui_text_box.settings_item
        current_config_value = getattr(self.config, config_item_name)
        if current_config_value != ui_text_box.Value:
            try:
                new_value = ui_text_box.type_caster(ui_text_box.Value)
            except ValueError as _e:
                print(f"Could not parse configuration setting {ui_text_box.Value}")
                return False
            self.config.save_value(config_item_name, new_value)
            return True
        return False

    def on_close(self, _event):
        ''' Force exit of all threads '''
        self.Destroy()
        sys.exit(0)  # TODO: figure out why processing threads are otherwise preventing a full close

    def __getstate__(self):
        ''' Force pickling of the UI to fail '''
        logging.debug("*Attempting to pickle main UI*")
        return None  # used to locate pickle attempts on this object
