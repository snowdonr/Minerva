'''
Created on Oct 2, 2020

@author: Ryan
'''
from __future__ import annotations
import wx
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot
import matplotlib.axes
import numpy
import operator
import logging
from . import samplefrom . import projectfrom . import mass_spectrum

USE_BOXPLOT_ONLY = False
QC_PLOT_MODE = False


class GeneralPlotFrame(wx.Frame):
    '''
    Base class for all plots
    '''
    TIME_NORM = 60.0

    def __init__(self, title):
        super().__init__(None, wx.ID_ANY, title=title)
        self.create_main_panel()
        self.Layout()
        self.Show()

    def create_main_panel(self, skip_buttons=False):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        # Create the mpl Figure and FigCanvas objects.
        self.dpi = 72
        self.fig = Figure((14.0, 12.0), dpi=self.dpi)
        self.canvas = FigureCanvasWxAgg(self, -1, self.fig)

        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.plots = []
        self.axes = self.fig.add_subplot(1, 1, 1)
        self.plots.append(self.axes)
        self.current_plot = self.plots[0]
        self.fig.tight_layout()

        # Bind the 'pick' event for clicking on one of the bars
        self.canvas.mpl_connect('pick_event', self.on_pick)

        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)
        self.toolbar.Realize()

        self.button_sizer = wx.BoxSizer(wx.HORIZONTAL)  # This sizer is also used by inherited classes
        if not skip_buttons:
            previous_button = wx.Button(self, wx.ID_ANY, "Previous")
            self.Bind(wx.EVT_BUTTON, self.on_previous, previous_button)
            self.button_sizer.Add(previous_button, 0, flag=wx.ALIGN_LEFT)
            next_button = wx.Button(self, wx.ID_ANY, "Next")
            self.Bind(wx.EVT_BUTTON, self.on_next, next_button)
            self.button_sizer.Add(next_button, 0, flag=wx.ALIGN_LEFT)

        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.Add(self.button_sizer, 0, wx.EXPAND)
        self.SetSizer(self.vbox)
        self.vbox.Fit(self)

    def on_next(self, _event):
        print("No next button implementation")

    def on_previous(self, _event):
        print("No previous button implementation")

    def on_pick(self, event):
        print(f"No pick implemented for {event}")

    def update_plot(self):
        self.Show(show=True)
        self.canvas.draw()

    def on_close(self, _event):
        self.Show(show=False)


class SampleSetPlot(GeneralPlotFrame):
    '''
    Summary plot of select data from all samples
    '''
    def __init__(self, main_project: project.Project):
        super().__init__("Sample Summary")
        self._project = main_project
        self._detail_plot = PeakDetailPlots(main_project)
        self._last_source_data = None
        self._last_data_filtered = None
        self._last_x_values = None

    def create_main_panel(self):
        super().create_main_panel(skip_buttons=True)
        sizer = self.button_sizer
        label_text = wx.StaticText(self, -1, "     Get details from RT:")
        sizer.Add(label_text, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        self.select_rt = wx.TextCtrl(self, value="0.0", size=(80, -1), style=wx.TE_PROCESS_ENTER)
        sizer.Add(self.select_rt, 0, flag=wx.ALIGN_LEFT)
        detail_button = wx.Button(self, wx.ID_ANY, "Update Details")
        self.Bind(wx.EVT_BUTTON, self.on_details, detail_button)
        sizer.Add(detail_button, 0, flag=wx.ALIGN_LEFT)
        self.Layout()
        self.Show()

    def on_details(self, _event):
        if self._last_data_filtered is None:
            return
        try:
            user_rt_value = float(self.select_rt.Value)*self.TIME_NORM
        except Exception as _e:
            print("Could not interpret rt plot selection")
            return
        insertion_index = numpy.searchsorted(self._last_x_values, user_rt_value)
        if insertion_index == 0:
            user_index = insertion_index
        elif insertion_index == len(self._last_x_values):
            user_index = insertion_index-1
        elif insertion_index > len(self._last_x_values):
            logging.warning(f"Details read out of bounds {insertion_index} > {len(self._last_x_values)}")
            user_index = len(self._last_x_values)-1
        else:
            lower_point = self._last_x_values[insertion_index-1]
            higher_point = self._last_x_values[insertion_index]
            user_index = (insertion_index-1) if abs(user_rt_value-lower_point) < abs(user_rt_value-higher_point) else insertion_index
        self._detail_plot.update_data(self._last_x_values, self._last_data_filtered)
        self._detail_plot.show_data(user_index)

    def on_next(self, event):
        self._detail_plot.update_data(self._last_x_values, self._last_data_filtered)
        self._detail_plot.on_next(event)

    def on_previous(self, event):
        self._detail_plot.update_data(self._last_x_values, self._last_data_filtered)
        self._detail_plot.on_previous(event)

    def show_plot(self, aligned_data: list):
        max_peaks_to_show = self._project.config.detail_plot_peaks
        log_scale = True
        self.current_plot.clear()
        self._last_source_data = aligned_data
        # alternatively transpose and use pandas box plot
        area_lists = [[0 if x is None else x.area for x in aligned_row] for aligned_row in aligned_data]

        area_grid = numpy.array(area_lists)
        area_max_set = area_grid.max(axis=1).argsort()
        max_filter = area_max_set[-max_peaks_to_show:]  # Ordered smallest to largest

        increasing_filtered_areas = area_grid[max_filter]
        rt_lists = [[numpy.NaN if x is None else x.rt for x in aligned_row] for aligned_row in aligned_data]
        rt_grid = numpy.array(rt_lists)
        rt_list = numpy.nanmean(rt_grid, axis=1).T
        rt_filter = rt_list[max_filter]
        rt_sorting = rt_filter.argsort()
        filtered_areas = increasing_filtered_areas.T[:, rt_sorting]
        if log_scale:
            filtered_areas = numpy.log10(filtered_areas+1)

        self._last_x_values = rt_filter[rt_sorting]
        self._last_data_filtered = numpy.zeros(max_peaks_to_show, dtype=object)
        for i, sorted_index in enumerate(max_filter[rt_sorting]):
            self._last_data_filtered[i] = aligned_data[sorted_index]
        if USE_BOXPLOT_ONLY:
            self.current_plot.boxplot(filtered_areas)
            tick_labels = [f"{(x/self.TIME_NORM):.2f}" for x in self._last_x_values]
            xtick_names = matplotlib.pyplot.setp(self.current_plot, xticklabels=tick_labels)
            matplotlib.pyplot.setp(xtick_names, rotation=45, fontsize=8)
        else:
            sample_labels = [x.sample_number for x in self._project.sample_set]
            for row in range(filtered_areas.shape[0]):
                self.current_plot.scatter(self._last_x_values/self.TIME_NORM, filtered_areas[row, :], label=sample_labels[row], picker=5)
            # self.current_plot.errorbar(self._last_x_values, y, yerr=bar_widths)
            x_list = []
            y_list = []
            for sample_index in range(len(self._project.sample_set)):
                input_ic = self._project.sample_set.get_tic(sample_index)
                x_list.append(numpy.array(input_ic.time_list)/self.TIME_NORM)
                y_list.append(numpy.array(input_ic.intensity_array))

            new_size = int(numpy.median([x.shape[0] for x in x_list]))
            [x.resize(new_size, refcheck=False) for x in x_list]  # In place truncate, affects all views
            [x.resize(new_size, refcheck=False) for x in y_list]
            average_x = numpy.mean(x_list, axis=0)[:-3]  # TODO: Sort out edge effects instead of truncating
            average_y = numpy.mean(y_list, axis=0)[:-3]
            if log_scale:
                average_y = numpy.log10(average_y+1)

            self.current_plot.plot(average_x, average_y)
            self.current_plot.boxplot(filtered_areas, positions=self._last_x_values/self.TIME_NORM, notch=False, widths=0.03, manage_ticks=False)
        try:
            self.current_plot.set_title(f"Summarized peak areas for {len(aligned_data[0])} samples")
        except Exception as _e:
            pass
        self.current_plot.set_xlabel("RT (min)")
        self.current_plot.set_ylabel("Log$_{10}$(Intensity) or Log$_{10}$(Area)")
        self.update_plot()
        self.on_details(None)

    def on_pick(self, event):
        ''' Selection of a point on the plot '''
        if event.mouseevent.button == "up" or event.mouseevent.button == "down":
            return False  # No pick on scroll wheel
        targetIndices = event.ind
        if len(targetIndices) < 1:
            return True
        self._detail_plot.update(event)


class PeakDetailPlots(GeneralPlotFrame):
    '''
    A set of plots related to a single aligned peak set
    '''
    def __init__(self, main_project: project.Project):
        super().__init__("Aligned Peak sets")
        self._project = main_project
        self._sample_set = main_project.sample_set
        self._current_rt_index = 0
        self.plot_count = 4
        self._source_data = None
        self._rt_list = []
        self._select_sample_index = 0
        self._override_compound_update = False
        del self.current_plot  # TODO: Clean up inheritance
        self.fig.clear()
        count_width = 2
        self.normalize_area = True
        axes1 = self.fig.add_subplot(self.plot_count//count_width, count_width, 1)
        axes2 = self.fig.add_subplot(self.plot_count//count_width, count_width, 2)
        if QC_PLOT_MODE:
            axes3 = self.fig.add_subplot(self.plot_count//count_width, count_width, 3, sharex=axes1)
        else:
            axes3 = self.fig.add_subplot(self.plot_count//count_width, count_width, 3)
        axes4 = self.fig.add_subplot(self.plot_count//count_width, count_width, 4, sharex=axes2)
        self.plots = [axes1, axes2, axes3, axes4]

    def create_main_panel(self):
        super().create_main_panel()
        sizer = self.button_sizer
        label_text = wx.StaticText(self, -1, "     Select compound")
        sizer.Add(label_text, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        self.select_compound = wx.ComboBox(self, id=wx.ID_ANY, value="Select_Compund_Value_Test", size=(240, -1), choices=["All"], name="select_compound")  # style=wx.CB_SIMPLE, 
        self.select_compound.Select(0)
        sizer.Add(self.select_compound, 0, flag=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_COMBOBOX, self.on_change_compound, self.select_compound)

        # self.confirm_compound = wx.Button(self, wx.ID_ANY, "Confirm compound id")
        # sizer.Add(self.confirm_compound, 0, flag=wx.ALIGN_LEFT)
        # self.Bind(wx.EVT_BUTTON, self.on_confirm_compound, self.confirm_compound)
        self.Layout()
        self.Show()

    def update_data(self, rt_list: list, data: numpy.ndarray):
        ''' Change source data for all potential peak sets '''
        self.Show()  # Display this window, if it has not been used yet this session
        self._rt_list = rt_list
        self._source_data = data

    def show_data(self, target_rt_index: int):
        ''' Set a new aligned peak set, and update all subplots from the largest area sample in the set '''
        if target_rt_index < 0:
            target_rt_index = len(self._rt_list)-1
        elif target_rt_index >= len(self._rt_list):
            target_rt_index = 0
        self._current_rt_index = target_rt_index
        data = self._source_data[target_rt_index]

        self._param_plot(self.plots[0], data, target_rt_index, "area", callback_target=self.plots[1])

        max_sample_index = numpy.argmax([0 if x is None else x.area for x in data])
        self._select_sample_index = max_sample_index
        self._source_peaks = data
        self.update_sample_subplots(max_sample_index, target_rt_index)
        self.update_plot()

    def update_sample_subplots(self, sample_index, target_rt_index):
        ''' Update 3 sub-plots with data from a single sample '''
        frag_mass_range = self._mass_spec_plot(self.plots[1], sample_index)
        library_mass_range = self._library_plot(self.plots[3], self._source_peaks, sample_index)
        if self._project.config.restrict_to_mass_scan:
            min_mass = frag_mass_range[0]
        else:
            min_mass = frag_mass_range[0] if frag_mass_range[0] < library_mass_range[0] else library_mass_range[0]
        max_mass = frag_mass_range[1] if frag_mass_range[1] > library_mass_range[1] else library_mass_range[1]
        self.plots[1].set_xlim(min_mass, max_mass)
        data = self._source_data[target_rt_index]
        if QC_PLOT_MODE:
            self._param_plot(self.plots[2], data, target_rt_index, "rt")
        else:
            self._rt_plot(self.plots[2], target_rt_index)

    def show_single_data(self, current_plot, source_reference):
        current_plot.clear()
        library_ms = source_reference.mass_spec
        # compare_scores = mass_spectrum.MassSpectrum.compare_to_target_ms(peak_data.mass_spectrum, library_ms)

        legend_str = f"{source_reference.name}({source_reference.formula}) mw:{source_reference.mw}"  # CAS {reference.cas}
        intensities = library_ms.intensity_list
        masses = library_ms.mass_list
        _plt = current_plot.stem(masses, intensities, label=legend_str)

        self._annotate_largest(current_plot, masses, intensities)
        current_plot.legend()
        self.update_plot()

    def remove_common_string_parts(self, string_list):
        pre_count=0
        first_string = str(string_list[0])
        for i in range(1, len(first_string)):
            pre = first_string[:i]
            if False in [str(path_string).startswith(pre) for path_string in string_list]:
                pre_count = i-1
                break
        post_count = 0
        for i in range(1, len(first_string)):
            post = first_string[-i:]
            if False in [str(path_string).endswith(post) for path_string in string_list]:
                post_count = i-1
                break
        pathValues = [str(path)[pre_count:len(str(path))-post_count] for path in string_list]
        return pathValues

    def _param_plot(self, current_plot, data, target_rt_index, parameter_name, callback_target=None):
        current_plot.clear()
        target_rt = self._rt_list[target_rt_index]
        sample_labels = [f"{x.sample_number}" for x in self._project._analysis_set._alignment_data.remaining_samples]
        sample_labels = self.remove_common_string_parts(sample_labels)

        param_get = operator.attrgetter(parameter_name)
        y_values = numpy.array([numpy.NaN if x is None else param_get(x) for x in data])
        if self.normalize_area:
            sample_areas = numpy.array([x.total_area for x in self._project._analysis_set._alignment_data.remaining_samples])
            if all(sample_areas > 0):
                y_values = 1000*y_values/sample_areas
            else:
                logging.error(f"Total area {sample_areas}")

        plt = current_plot.plot(sample_labels, y_values, marker='x', picker=5)
        plt[0].callback_plot = callback_target

        tick_locations = current_plot.get_xticks()
        current_plot.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(tick_locations))
        xtick_names = matplotlib.pyplot.setp(current_plot, xticklabels=sample_labels)
        matplotlib.pyplot.setp(xtick_names, rotation=45, fontsize=8)
        current_plot.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(sample_labels))  # xtick_names

        sample_count = numpy.count_nonzero(data)
        title = f"{parameter_name} at RT {target_rt/self.TIME_NORM:0.2f} with {sample_count}/{len(data)} peaks"
        current_plot.grid(which='major', axis="x", alpha=0.5, linewidth=0.5)
        current_plot.set_xlabel("Samples")
        current_plot.set_ylabel(f"{parameter_name}")
        current_plot.set_title(title)

    def _rt_plot(self, current_plot, target_rt_index: int):
        current_plot.clear()
        input_ic = self._project.sample_set.get_tic(self._select_sample_index)
        x_values = numpy.array(input_ic.time_list)
        y_values = numpy.array(input_ic.intensity_array)
        current_plot.plot(x_values/self.TIME_NORM, y_values)

        peak_data = self._source_data[target_rt_index]
        target_rt = peak_data.sql_ref.average_rt
        target_intensity = y_values[numpy.searchsorted(x_values, target_rt)]
        target_rt_norm = target_rt/self.TIME_NORM
        current_plot.annotate(f"RT {target_rt_norm:0.3f}",
                              xy=(target_rt_norm, target_intensity), xycoords='data',
                              xytext=(target_rt_norm*1.01, target_intensity*1.1), textcoords='data',
                              arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
                              )
        sample_id = self._sample_set[self._select_sample_index].sample_number
        current_plot.set_title(f"TIC for {sample_id}")
        current_plot.set_xlabel("RT (min)")
        current_plot.set_ylim(0, min(target_intensity*2.0, max(y_values)*1.05))

    def _mass_spec_plot(self, current_plot: matplotlib.axes.SubplotBase, sample_index: int):
        current_plot.clear()
        sample_id = self._sample_set[sample_index].sample_number
        peak_data = self._source_peaks[sample_index]

        _stem_plot = current_plot.stem(peak_data.mass_spectrum.mass_list, peak_data.mass_spectrum.intensity_list)
        current_plot.plot(peak_data.ion_areas.keys(), numpy.array(list(peak_data.ion_areas.values()))/8.0, 'x')
        
        self._annotate_largest(current_plot, peak_data.mass_spectrum.mass_list, peak_data.mass_spectrum.intensity_list)
        current_plot.set_title(f"Mass Spectrum for {sample_id}")
        current_plot.set_xlabel("Mass")
        current_plot.set_ylabel("Intensity")
        valid_points = numpy.nonzero(peak_data.mass_spectrum.intensity_list)[0]
        mass_range = [peak_data.mass_spectrum.mass_list[0], peak_data.mass_spectrum.mass_list[valid_points[-1]]]  # Always show the bottom limit
        return mass_range

    def _library_plot(self, current_plot, peak_data, selected_sample_index):
        current_plot.clear()
        # found_reference_data = self._project.nist_library.search(peak_data.mass_spectrum, peak_data.rt, max_hits=5)
        # user_entries = self._project.user_library.search(peak_data.mass_spectrum, peak_data.rt, max_hits=5)
        # if user_entries is not None:
        #    found_reference_data.extend(user_entries)
        # self.mass_plot.show_peak_mass_spectrum(selected_peak, selected_sample, normalize=self.NORM_MAX, first_plot=True)
        found_reference_data = peak_data.compound_matches
        if len(found_reference_data) < 1:
            return [0, 9]
        library_ms = [x.compound.mass_spec for x in found_reference_data]
        sample_mass_spec = peak_data[selected_sample_index].mass_spectrum
        compare_scores = mass_spectrum.MassSpectrum.compare_to_target_ms(sample_mass_spec, library_ms)
        plot_list = []
        all_masses = []
        all_intensities = []
        marker_list = ["x", "+", "*", "v", "^", "o"]
        for index, reference in enumerate(found_reference_data):
            ms = reference.compound.mass_spec

            legend_str = f"{100*compare_scores[index]:.1f} {reference.compound.name}({reference.compound.formula}) mw:{reference.compound.mw}"  # CAS {reference.cas}
            intensities = ms.intensity_list
            current_marker = marker_list[index % len(marker_list)]
            plt = current_plot.scatter(ms.mass_list, intensities, marker=current_marker, s=25, label=legend_str)
            if numpy.any(numpy.isnan(intensities)):
                logging.info(f"NaN: Unknown mass spectrum for {reference.compound.name}")
            all_masses.extend(ms.mass_list)
            all_intensities.extend(intensities)
            plot_list.append(plt)
        self._select_compound_data = [x.compound for x in found_reference_data]
        string_list = ["All", ]
        string_list.extend([x.compound.name if x is not None else "?" for x in found_reference_data])
        self.select_compound.SetItems(string_list)
        # self._override_compound_update not required
        self.select_compound.Select(0)
        self._annotate_largest(current_plot, all_masses, all_intensities)
        current_plot.legend()
        if len(all_masses) < 1:
            return [0, 10]
        else:
            return [min(all_masses), max(all_masses)]

    def _annotate_largest(self, current_plot, x_data: list, y_data: list, max_count: int=10):
        index_largest_y = numpy.argsort(y_data)[::-1]
        annotate_count = 0  # Tried to make this with numpy.unique on x; non-trivial
        shown_x = []
        for index in index_largest_y:
            try:
                current_x = x_data[index]
                if numpy.isnan(y_data[index]) or numpy.isnan(current_x):
                    print(f"Mass unknown in annotation for x,y {current_x}, {y_data[index]}")
                    continue
                if current_x not in shown_x:
                    current_plot.annotate(f"{x_data[index]}", (current_x, y_data[index]))
                    shown_x.append(current_x)
                    annotate_count += 1
                if annotate_count >= max_count:
                    break
            except Exception as _e:
                pass

    def update(self, pick_event: wx.Event):
        pass

    def on_change_compound(self, event: wx.Event):
        if self._override_compound_update:
            return
        if event.GetSelection() == 0:
            _library_mass_range = self._library_plot(self.plots[3], self._source_peaks, self._select_sample_index)
            self.show_data(self._current_rt_index)
            return
        source_reference = self._select_compound_data[event.GetSelection()-1]
        self.show_single_data(self.plots[3], source_reference)

    def on_confirm_compound(self, _event: wx.Event):
        try:
            selected_index = self.select_compound.GetSelection()-1
            if selected_index < 0:
                print("Could not confirm: No selection")
                return
            if selected_index >= len(self._select_compound_data):
                print("DEBUG: Selection out of range, could not confirm")
                return
            selected_raw_data = self._select_compound_data[selected_index]

            aligned_peak = self._project.aligned_peaks[self._current_rt_index]

            # self._select_sample_index
            print("Unimplemented: Adding confirmation for compound")
            aligned_peak.update_compound(selected_raw_data, "UserConfirmed-Minerva")
        except Exception as _e:
            logging.exception("Confirm function failed")

    def on_pick(self, event: wx.Event):
        event.artist.callback_plot
        self._override_compound_update = True
        self._select_sample_index = event.ind[0]
        self.update_sample_subplots(self._select_sample_index, self._current_rt_index)
        self.update_plot()
        self._override_compound_update = False

    def on_next(self, _event: wx.Event):
        self.show_data(self._current_rt_index+1)

    def on_previous(self, _event: wx.Event):
        self.show_data(self._current_rt_index-1)


# class FeaturePlot(GeneralPlotFrame):
#     ''' Compare meta-data for multiple samples '''
#     def __init__(self):
#         super().__init__("Feature Detail Plot")
#         self.pair_length = 3
#         del self.current_plot  # Change inheritance
#         self.fig.clear()
#         self.plots = []
#         for index in range(self.pair_length*2):
#             axes = self.fig.add_subplot(self.pair_length, 2, index+1)
#             self.plots.append(axes)
#         self.select_subplot(0)
#         self.first_feature_display = 0
#
#     def select_subplot(self, plot_number):
#         if plot_number < self.pair_length and plot_number >= 0:
#             self.selected_plot_index = plot_number
#             self.current_plot_left = self.plots[plot_number*2]
#             self.current_plot_right = self.plots[plot_number*2+1]
#         else:
#             print(f"No plots available for display at {plot_number}")
#
#     def clear_feature(self):
#         self.current_plot_left.clear()
#         self.current_plot_right.clear()
#
#     def load_all_features(self, sample_set, feature_sources, regression_response_list):
#         self.first_feature_display = 0
#         self.sample_set = sample_set
#         self.dose_set = [current_sample.dose for current_sample in sample_set]
#         self.feature_sources = feature_sources
#         self.feature_regression_response = regression_response_list
#
#     def update_feature_set(self, feature_index):
#         self.first_feature_display = feature_index
#         for i in range(self.pair_length):
#             if feature_index+i < 0 or feature_index+i >= len(self.feature_sources):
#                 print(f"No data to display at feature {feature_index+i}")
#                 continue
#             self.select_subplot(i)
#             self.clear_feature()
#             source_data = self.feature_sources[feature_index+i]
#             try:
#                 mass_data_set = isinstance(source_data[0], int)
#             except Exception as _e:
#                 mass_data_set = False
#
#             if mass_data_set:
#                 self.current_plot_left.plot(self.dose_set, source_data[1], "x-", label="Total Ion Area")
#                 self.current_plot_left.set_title(f"Feature {feature_index+i}: Total mass at {source_data[0]}")
#                 self.current_plot_left.set_xlabel("Dose")
#                 self.current_plot_left.set_ylabel("Total Fragment Area")
#                 self.current_plot_right.clear()  # Redundant, replace with some kind of plot?
#                 self.current_plot_right.set_title(f"Feature coef={self.feature_regression_response[feature_index+i]:.5f}")
#             else:
#                 peak_set = source_data
#                 areas = [sample.Sample.get_area(peak) for peak in peak_set]
#                 self.current_plot_left.plot(self.dose_set, areas, "x-", label="Ion Area")
#                 try:
#                     rt = numpy.average([peak.rt for peak in peak_set if peak is not None])  # TODO: Check deviation
#                     self.current_plot_left.set_title(f"Feature {feature_index+i}: Peak at time {rt/self.TIME_NORM:.3f}")
#                 except Exception as _e:
#                     print("WARNING: Could not add RT to plot")
#                 self.current_plot_left.set_xlabel("Dose")
#                 self.current_plot_left.set_ylabel("Ion Area")
#
#                 self._used_annotations = []  # Avoid labelling the same masses multiple times?
#                 for current_sample, peak in zip(self.sample_set, peak_set):
#                     if peak is not None:
#                         # self.current_plot_right.plot(peak.mass_spectrum.mass_list, peak.mass_spectrum.intensity_list)
#                         x_set, y_set = numpy.array(list(peak.ion_areas.keys())), numpy.array(list(peak.ion_areas.values()))
#                         self.current_plot_right.scatter(x_set, y_set, label=f"{current_sample.sample_number}-{current_sample.dose}")
#                         annotate_index = numpy.argsort(y_set)[-3:][::-1]
#                         self._annotate_x(x_set[annotate_index], y_set[annotate_index])
#                 self.current_plot_right.legend()
#                 self.current_plot_right.set_title(f"Feature coef={self.feature_regression_response[feature_index+i]:.5f}. Corresponding mass spectrum:")
#                 self.current_plot_right.set_xlabel("Mass")
#                 self.current_plot_right.set_ylabel("Area")  # Intensity
#
#     def _annotate_x(self, x_array, y_array):
#         for x, y in zip(x_array, y_array):
#             if x not in self._used_annotations:
#                 self._used_annotations.append(x)  # TODO: Swap to the most intense peak
#                 self.current_plot_right.annotate(x, (x, y), xycoords='data')
#
#     def on_next(self, _event):
#         self.update_feature_set(self.first_feature_display+self.pair_length)
#         self.update_plot()
#
#     def on_previous(self, _event):
#         self.update_feature_set(self.first_feature_display-self.pair_length)
#         self.update_plot()


class LoadedDataPlot(GeneralPlotFrame):
    '''
    Plots raw GCMS TIC data for verification and testing
    '''
    def __init__(self, mass_spec_frame, sample_set=None):
        super().__init__(title='GCMS TIC plot')
        self._mass_spec_frame = mass_spec_frame
        self.sample_set = sample_set
        self._current_sample_count = -1
        self._source_peak_list = []  # Used for indexing data for sub plots

    def on_next(self, _event):
        self.change_plot(1)

    def on_previous(self, _event):
        self.change_plot(-1)

    def change_plot(self, change_index: int, maintain_range: bool=True):
        if self.sample_set is None:
            print("No loaded samples for display")
            return
        self._current_sample_count += change_index
        if self._current_sample_count >= len(self.sample_set):
            self._current_sample_count = 0
        elif self._current_sample_count < 0:
            self._current_sample_count = len(self.sample_set)-1
        x_range = self.current_plot.get_xlim()
        self.current_plot.clear()
        self.show_ion_chromatogram(self.sample_set.get_tic(self._current_sample_count))
        self.show_peak_total_intensity(self.sample_set[self._current_sample_count])
        if maintain_range:
            self.current_plot.set_xlim(x_range)
        self.update_plot()

    def on_pick(self, event):
        ''' Selection of a point on the plot '''
        if event.mouseevent.button == "up" or event.mouseevent.button == "down":
            return False  # No pick on scroll wheel
        targetIndices = event.ind
        if len(targetIndices) < 1:
            return True
        peak_index = event.ind[0]
        self._mass_spec_frame.load_peak_list(event.artist.active_sample, self._source_peak_list)
        self._mass_spec_frame.show_peak_mass_spectrum(peak_index)
        print(f"{event.mouseevent.xdata}  {self._source_peak_list[peak_index]}")
        return True

    def show_ion_chromatogram(self, input_ic):
        ''' Plot time list and intensity array of input_ic pyms chomatogram '''
        self.current_plot.plot(numpy.array(input_ic.time_list)/self.TIME_NORM, input_ic.intensity_array)
        x_data = numpy.array(input_ic.time_list)/self.TIME_NORM
        y_data = [self.sample_set.baseline_at(x) for x in numpy.array(input_ic.time_list)]
        raw_y_data = self.sample_set.full_baseline

        self.current_plot.plot(x_data, y_data, "r")
        if not raw_y_data.any():
            # no baseline
            return
        if len(x_data) > len(raw_y_data):
            raw_x_data = numpy.resize(x_data, len(raw_y_data))
        else:
            raw_y_data = numpy.resize(raw_y_data, len(x_data))
            raw_x_data = x_data
        self.current_plot.plot(raw_x_data, raw_y_data, "g")
        # self.update_plot()  # Required if not using show_peak_total_intensity

    def show_peak_total_intensity(self, sample):
        ''' Scatter plot: sum of ion intensity for each peak '''
        peak_list = sample.peak_list
        plot_x = []
        plot_y = []
        for peak in peak_list:
            plot_x.append(peak.rt/self.TIME_NORM)
            plot_y.append(sum(peak.mass_spectrum.intensity_list))
        self._source_peak_list = peak_list
        plt = self.current_plot.scatter(plot_x, plot_y, marker="+", c="blue", zorder=20, s=30, picker=5)
        plt.active_sample = sample
        plt.x_data = plot_x
        self.current_plot.set_title(f"GCMS Load {sample}")
        self.current_plot.set_xlabel("RT (min)")
        self.current_plot.set_ylabel("Intensity")


class MassSpecPlot(GeneralPlotFrame):
    def __init__(self, chromatogram_detail, title='Mass Spectrum'):
        super().__init__(title=title)
        self._chromatogram_plot = chromatogram_detail
        self._last_index = 0
        self._peak_list = []
        self._source_sample = None

    def load_peak_list(self, source_sample, source_peak_list):
        self._source_sample = source_sample
        self._peak_list = source_peak_list

    def get_index_by_rt(self, rt_value: float) -> int:
        rt_list = numpy.array([x.rt for x in self._peak_list])
        insert_index = numpy.searchsorted(rt_list, rt_value, side='left')
        dif_right = abs(rt_list[insert_index] - rt_value)
        try:
            dif_left = abs(rt_list[insert_index-1] - rt_value)
            if dif_left < dif_right:
                return insert_index-1
        except IndexError:
            pass
        return insert_index

    def show_peak_mass_spectrum(self, source_peak_index, normalize: float=False, first_plot: bool=True):
        '''
        source_peak_index
        active_sample: sample associated with this source_peak for title and picker linking
        normalize: scale the largest source_peak to this value
        first_plot: reset figure before plotting
        '''
        if first_plot:
            self.current_plot.clear()
        if self._source_sample is None:
            print("No data for MassSpecPlot")
            return
        if source_peak_index >= len(self._peak_list):
            source_peak_index = 0
        elif source_peak_index < 0:
            source_peak_index = len(self._peak_list)-1
        self._last_index = source_peak_index
        source_peak = self._peak_list[source_peak_index]

        ms = source_peak.mass_spectrum
        intensities = ms.intensity_list
        if normalize:
            intensities = normalize*numpy.array(intensities)/max(intensities)
        plt = self.current_plot.scatter(ms.mass_list, intensities, marker="*", s=20, picker=5)
        self.current_plot.stem(ms.mass_list, intensities)
        try:
            point_count = len(source_peak.ion_areas)  # TODO: Filter by intensity minimum
        except Exception as _e:
            point_count = "Unknown"
        if first_plot:
            self.current_plot.set_title(f"{self._source_sample} RT {source_peak.rt/self.TIME_NORM:.3f} with {point_count}")
        else:
            old_title = self.current_plot.get_title()
            self.current_plot.set_title(old_title+f" and {source_peak.rt/self.TIME_NORM:.3f} with {point_count}")
        plt.active_sample = self._source_sample
        plt.x_data = ms.mass_list
        self.current_plot.set_xlabel("Mass")
        self.current_plot.set_ylabel("Intensity")
        self.update_plot()
        return source_peak

#     def show_reference_spectrum(self, found_reference_data, normalize: float=False, first_plot=False):
#         if first_plot:
#             self.current_plot.clear()
#         plot_list = []
#         for _compound_name, reference in found_reference_data:
#             ms = reference.mass_spec
#             legend_str = f"{reference.name}({reference.formula}) CAS {reference.cas}"
#             intensities = ms.intensity_list
#             if normalize:
#                 intensities = normalize*numpy.array(intensities)/max(intensities)
#             plt = self.current_plot.scatter(ms.mass_list, intensities, marker="x", s=25, label=legend_str)
#             plot_list.append(plt)
#
#         self.current_plot.legend()
#         self.current_plot.set_xlabel("Mass")
#         self.current_plot.set_ylabel("Intensity")
#         self.update_plot()

    def on_pick(self, event):
        selected_sample = event.artist.active_sample
        try:
            selection_indices = event.ind
            mass_set = numpy.take(event.artist.x_data, selection_indices)
            selected_mass = mass_set[0]
        except Exception as _e:
            selected_mass = event.mouseevent.xdata
        if self._chromatogram_plot:
            print(f"Showing chromatogram at {selected_mass} for sample {selected_sample}")
            self._chromatogram_plot.show_mass(selected_sample, selected_mass)
        else:
            print(f"Selected {selected_mass} for sample {selected_sample}")

    def on_next(self, _event):
        self.show_peak_mass_spectrum(self._last_index+1)

    def on_previous(self, _event):
        self.show_peak_mass_spectrum(self._last_index-1)


class SingleChromatogramPlot(GeneralPlotFrame):
    ''' '''
    def __init__(self, parent_plot):
        super().__init__(title='Single Chromatogram')
        self._mass_plot = parent_plot
        self._last_sample = None
        self._last_mass = 0

    def show_mass(self, selected_sample, selected_mass):
        self.current_plot.clear()
        try:
            if selected_mass > selected_sample.max_mass:
                selected_mass = selected_sample.min_mass
            elif selected_mass < selected_sample.min_mass:
                selected_mass = selected_sample.max_mass
        except AttributeError:
            pass  # Older versions do not have this attribute
        self._last_sample = selected_sample
        self._last_mass = selected_mass
        try:
            ion_chromatogram = selected_sample.im.get_ic_at_mass(selected_mass)
        except AttributeError:
            logging.warning("Detailed data has been removed, ion chromatogram is not possible. Disable 'limit_loading_memory'")
            return
        self.current_plot.plot(numpy.array(ion_chromatogram.time_list)/self.TIME_NORM, ion_chromatogram.intensity_array)
        # list_mass = round(selected_mass)
        mass_index = selected_sample.im.get_index_of_mass(selected_mass)
        single_mass_peaks = [peak for peak in selected_sample.peak_list if peak.mass_spectrum.intensity_list[mass_index] > 0.0]  # .get_intensity_for_mass(list_mass)]
        y_data = [peak.mass_spectrum.intensity_list[mass_index] for peak in single_mass_peaks]
        x_data = [peak.rt/self.TIME_NORM for peak in single_mass_peaks]
        plt = self.current_plot.scatter(x_data, y_data, marker="x", s=20, picker=5)
        plt.x_data = x_data
        plt.profile_y_data = ion_chromatogram.intensity_array
        self.current_plot.set_title(f"{selected_sample} Mass {selected_mass:.3f} with {len(single_mass_peaks)}")
        self.current_plot.set_xlabel("RT (min)")
        self.current_plot.set_ylabel("Intensity")
        self.update_plot()

    def on_pick(self, event):
        try:
            selection_indices = event.ind
            rt_set = numpy.take(event.artist.x_data, selection_indices)
            selected_rt = rt_set[0]
        except Exception as _e:
            selected_rt = event.mouseevent.xdata
        if not self._mass_plot:
            return
        peak_index = self._mass_plot.get_index_by_rt(selected_rt*self.TIME_NORM)
        selected_peak = self._mass_plot.show_peak_mass_spectrum(peak_index)
        # TODO: Display integration area for peak, mass
        try:
            for ic_peak in selected_peak._ion_peaks:
                if ic_peak.mass == self._last_mass:
                    left_b = ic_peak.left
                    right_b = ic_peak.right
            center_b = selected_peak.bounds[1]
            left_t = self._last_sample.im.get_time_at_index(center_b-left_b)/self.TIME_NORM
            center_t = self._last_sample.im.get_time_at_index(center_b)/self.TIME_NORM
            right_t = self._last_sample.im.get_time_at_index(center_b+right_b)/self.TIME_NORM
            try:
                left_i = event.artist.profile_y_data[center_b-left_b]
                right_i = event.artist.profile_y_data[center_b+right_b]
            except Exception as _e:
                left_i = 0
                right_i = 0
            self.current_plot.plot([left_t, left_t, center_t, right_t, right_t], [left_i, 0, 0, 0, right_i], marker="v")
            self.update_plot()
        except Exception as _e:
            print("Could not update peak bounds")

    def on_next(self, _event):
        self.show_mass(self._last_sample, self._last_mass+1)

    def on_previous(self, _event):
        self.show_mass(self._last_sample, self._last_mass-1)
