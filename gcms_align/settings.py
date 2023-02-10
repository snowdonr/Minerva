'''
Created on Oct 16, 2020

@author: Ryan
'''
import configparser
import pathlib
import logging


class Setting(object):
    '''
    Central location for default values of user configuration options
    These can be overridden by the active configparser values of by call arguments
    Some parameters can be updated via UI and will be saved to configparser, not here
    '''

    settings_version = "0.8"
    method_analysis_name = "Requires Unique Analysis Name"
    experiment_name = "Experiment X"
    # **NOTE! Changes here may require deletion of invalid pickle temporary files**
    data_path = r"C:/GCMS/Minerva/SourceFolder.AIA"
    temp_path = r"C:/GCMS/Minerva/Data_Cache"
    pairwise_temp_file = r"_pairwise_relations.tmp"
    sample_relation_path = r"C:/GCMS/Minerva/GCMS_Sample_information.xlsx"
    internal_standard_path = r"C:/GCMS/Minerva/"
    user_library_location = "C:/GCMS/Minerva/Libraries/UserLibrary.xlsx"
    output_sqlite = True
    remote_db_location = "?.ca-central-1.rds.amazonaws.com"
    db_username = "postgres"
    db_password = ""
    db_name = "minerva"

    main_nist_path = r"C:/NIST11/MSSEARCH/mainlib"
    nist_working_dir_path = r"C:/GCMS/Minerva/Data_Cache"

    read_pickle = True
    create_pickle = True
    overwrite_pickle = False
    compress_cache_files = True
    force_load = False  # Even if pickles are present, load and process CDR files

    do_ion_calc = True  # Integrate areas for each ion individually
    do_peak_bounds_calc = True  # Store the beginning and end points of each peak area

    baseline_moving_average_width = 30
    baseline_multiplier = 1.0

    analysis_directory = r"C:/GCMS/Minerva/ExperimentData"
    # result_folder = r"E:/GCMS/Minerva"
    summary_file_path = r"C:/GCMS/Minerva/Summary.csv"
    tophat_baseline_noise_smooth = "1.5m"
    load_align_pickle = True
    bb_peak_points = 5
    bb_scan_merge = 2
    peak_trim_percent_of_max = 2
    min_ions = 3
    top_ion_areas = 10  # Number of top ions to record in peak
    # Use the same retention time range for all experiments
    lo_rt_limit = "80s"
    hi_rt_limit = "1500s"
    minimum_mass_limit = 31
    maximum_mass_limit = None
    alignment_minimum_peaks = 3
    peak_width_st_dev = 4  # "Peak width standard deviation in seconds. Determines search window width."
    # Define the within-state alignment parameters.
    rt_sensitivity_s = 2.5  # rt modulation [s]
    gap_penalty = 0.40

    align_multiprocess = False
    align_sparse_mode = True
    align_diagonal_width = 4
    align_full_compare = (0, 5, 10, 15, 20)
    align_score_full = False
    align_score_backtrack = 4  # Reverse check of the ~diagonal score matrix
    align_end_row_check = 80

    library_rt_scale = 2.5  # multiplier applied to rt_sensitivity for library scoring
    minimum_compound_score = 0.001  # score to list compound option in result
    spectrum_merge_required_match_portion = 0.3

    limit_loading_memory = False
    meta_data_excel = False
    merge_aligned_peaks = True
    merge_compare_normalize = True
    consider_rt_range = 15
    max_relation_grouping_depth = 3
    require_all_peaks = False  # Align result must have a peak in each sample
    divide_samples = False  # Enables normalization by one of the following methods:
    # divide_type = 'weights_only'  # mutually exclusive with other divide types
    # divide_type = 'sample_total'
    divide_type = 'standard_area'  # Assuming internal standards are working, this is what you should use
    detail_plot_peaks = 100

    # Name, Fragment masses and estimated retention times (in seconds) for the standards that are used in normalization
    # *Names must be unique*
    standard_list = ["Naphthalene-D8", ]
    # ("d10 Phenanthrene", 178+10, 2058),  # mass 97 cluttered
    # ("1,1 Binaphthyl", 254, 2550),  #
    # ("Squalane", 85, None),  # n-alkanes
    # ("d4 Cholestane", 217+4, 3853),  # mass 376
    # ("d16 Adamantane", 102, 1152),  # 136+16 rt 1152*, 1243, 1280  base fragments 79, 93
    # ("d30 Phenyldodacane", 138, 2642),  # mass 100 246+30 rt 2643*, 2878  Missing
    # ("d32 Pentadecane", 212+32, 2008),
    # ("d52 Pentacosane", 138, 3435)  # 352+52+1 rt 3435, very small  Missing

    standard_minimum_intensity = 1e4  # for selected fragment
    standard_maximum_rt_error = 120

    filter_minimum_area = 5e5
    manual_baseline_section_rt = 25  # minutes
    analyze_mass_sums = True
    analyze_peak_areas = True
    restrict_to_mass_scan = True  # Default plot zoom during identification does not show library below scan

    do_lasso = True
    do_lars_lasso = True  # Default fold requires >5 samples
    do_linear_svc = True
    do_bayesian_ridge = True

    lasso_alpha = 0.8  # Higher -> fewer features
    normalize_mean = True
    normalize_std = False
    profile_load = False
    blank_name_search = "qc-etoh"
    standard_name_search = "std-"

    def __init__(self, filename: pathlib.Path):
        '''
        Check configuration file to update values
        '''
        self._store = Settings_Store(filename)
        # check through attributes of the class, not instance
        for key, value in Setting.__dict__.items():
            # Ignore anything not of type str, float, int, or bool
            if isinstance(value, str):  # TODO: match value:   case str():  ...  case _: "Unknown type"
                read_value = self._store.read(key, None)
            elif isinstance(value, float):
                read_value = self._store.read_float(key, None)
            elif isinstance(value, bool):
                read_value = self._store.read_boolean(key, None)
            elif isinstance(value, int):  # must be after bool, as int is a supertype
                read_value = self._store.read_int(key, None)
            elif isinstance(value, list):
                read_value = self._store.read_list(key, None)
            else:
                logging.info(f"Unknown setting type {value}, {type(value)}")
                read_value = None

            if read_value is not None:
                setattr(self, key, read_value)  # set at the instance level, not the class
            elif not key.startswith("__"):
                logging.info(f"No value set for setting key {key}")

    def save_value(self, value_name: str, new_value: object):
        self._store.write(value_name, str(new_value))
        setattr(self, value_name, new_value)

    def commit(self):
        self._store.commit_file()


class Settings_Store(object):
    '''
    Interface for configparser storage of settings values
    '''
    def __init__(self, filename="Minerva.ini"):
        self._filename = filename
        self._config=configparser.ConfigParser()
        self.section_name = "Settings"
        success = self._config.read(self._filename)  # configparser read
        if not success:
            print(f"Failed to open settings at {filename}")

    def read(self, item_name, default_value):
        try:
            return self._config.get(self.section_name, item_name)
        except configparser.NoOptionError:
            return default_value
        except configparser.NoSectionError:
            return default_value
        except AttributeError:
            return default_value

    def write(self, item_name, new_value):
        self._changed = True
        try:
            return self._config.set(self.section_name, item_name, str(new_value))
        except configparser.NoSectionError:
            self._config.add_section(self.section_name)
            return self._config.set(self.section_name, item_name, str(new_value))

    def read_float(self, item_name, default_value):
        try:
            return self._config.getfloat(self.section_name, item_name)
        except (configparser.NoOptionError, configparser.NoSectionError, ValueError, AttributeError):
            return default_value

    def read_int(self, item_name, default_value):
        try:
            return self._config.getint(self.section_name, item_name)
        except (configparser.NoOptionError, configparser.NoSectionError, ValueError, AttributeError):
            return default_value

    def read_boolean(self, item_name, default_value):
        try:
            return self._config.getboolean(self.section_name, item_name)
        except (configparser.NoOptionError, configparser.NoSectionError, ValueError, AttributeError):
            return default_value

    def read_list(self, item_name, default_value):
        try:
            result = self._config.get(self.section_name, item_name, fallback=default_value)
            return result.split(",")
        except (configparser.NoOptionError, configparser.NoSectionError, ValueError, AttributeError):
            return default_value

    def commit_file(self):
        if not self._changed:
            return
        with open(self._filename, 'w') as configFile:
            self._config.write(configFile)
        self._changed = False
