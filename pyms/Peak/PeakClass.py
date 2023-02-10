"""
Provides a class to model signal peak.
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                              #
#    Copyright (C) 2023 Ryan Snowdon                                           #
#                                                                              #
#    This program is free software; you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License version 2 as         #
#    published by the Free Software Foundation.                                #
#                                                                              #
#    This program is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with this program; if not, write to the Free Software               #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                              #
################################################################################

# stdlib
import copy
import collections
from typing import List, Optional, Sequence, Tuple, Union, cast
from warnings import warn
import math

import numpy

# this package
from pyms.Base import pymsBaseClass
from pyms.IntensityMatrix import BaseIntensityMatrix
from pyms.Spectrum import MassSpectrum
from pyms.Utils.Utils import is_number, is_sequence

__all__ = ["AbstractPeak", "Peak", "ICPeak"]


class AbstractPeak(pymsBaseClass):
    """
    Models a signal peak.

    :param rt: Retention time.
    :param minutes: Retention time units flag. If :py:obj:`True`, retention time
        is in minutes; if :py:obj:`False` retention time is in seconds.
    :param outlier: Whether the peak is an outlier.

    :authors: Vladimir Likic, Andrew Isaac,
        Dominic Davis-Foster (type assertions and properties), David Kainer (outlier flag)

    .. versionadded:: 2.3.0
    """

    def __init__(
            self,
            rt: Union[int, float] = 0.0,
            minutes: bool = False,
            outlier: bool = False,
    ):

        if not is_number(rt):
            raise TypeError("'rt' must be a number")

        if minutes:
            rt = rt * 60.0

        # basic peak attributes
        self.is_outlier = outlier
        self._rt = float(rt)
        self._pt_bounds: Optional[Tuple[int, int, int]] = None
        self._area: Optional[float] = None

        self._UID = None

    def __eq__(self, other) -> bool:
        """
        Return whether this Peak object is equal to another object.

        :param other: The other object to test equality with.
        """

        if isinstance(other, self.__class__):
            return (
                self.UID == other.UID and self.bounds == other.bounds and self.rt == other.rt and
                self.area == other.area
            )

        return NotImplemented

    @property
    def area(self):
        return self._area

    @property
    def UID(self):
        if self._UID is not None:
            return self._UID
        else:
            self._make_UID()
            return self._UID

    def _make_UID(self):
        self._UID = f"{self._rt:.2f}"

    @property
    def bounds(self) -> Optional[Tuple[int, int, int]]:
        """
        The peak boundaries in points.

        :return: A 3-element tuple containing the left, apex, and right
            peak boundaries in points. Left and right are offsets.

        :author: Andrew Isaac
        """

        return self._pt_bounds

    @bounds.setter
    def bounds(self, value: Sequence[int]):
        """
        Sets peak boundaries in points.

        :param value: A 3-element tuple containing the left, apex, and right
            peak boundaries in points. Left and right are offsets.
        """

        if not is_sequence(value):
            raise TypeError("'Peak.bounds' must be a Sequence")

        if len(value) != 3:
            raise ValueError("'Peak.bounds' must have exactly 3 elements")

        for index, item in enumerate(value):
            if not isinstance(item, int):
                raise TypeError(f"'Peak.bounds' element #{index} must be an integer")

        self._pt_bounds = cast(Tuple[int, int, int], tuple(value[:3]))

    @property
    def rt(self) -> float:
        """
        The retention time of the peak, in seconds.
        """

        return self._rt

    def set_bounds(self, left: int, apex: int, right: int):
        """
        Sets peak boundaries in points.

        :param left: Left peak boundary, in points offset from apex
        :param apex: Apex of the peak, in points
        :param right: Right peak boundary, in points offset from apex
        """

        self.bounds = (left, apex, right)


class CompositePeak(AbstractPeak):
    """
    Peak created by combining a list of Peaks
    """
    def __init__(self, peak_list: list, minutes: bool = False):

        # if ignore_outliers:
        #     rts = []
        # if len(peak_list) > 3:
        #     for peak in peak_list:
        #         rts.append(peak.rt)
        #
        #     is_outlier = median_outliers(rts)
        #
        #     for i, val in enumerate(is_outlier):
        #         if val:
        #             peak_list[i].is_outlier = True

        rt_average = numpy.mean([x.rt for x in peak_list])
        if minutes:
            rt_average /= 60
        outlier = numpy.any([x.is_outlier for x in peak_list])
        self.is_outlier = outlier
        super().__init__(rt_average, minutes, outlier)
        combined_ions = collections.defaultdict(list)
        for current_peak in peak_list:
            for ion_peak in current_peak._ion_peaks:
                combined_ions[ion_peak.mass].append(ion_peak)

        self._ion_peaks = []
        self.ion_areas = {}
        for target_mass, ion_peak_list in combined_ions.items():
            new_ic_peak = ICPeak(rt_average, target_mass, minutes, outlier)
            self._combine_peaks(new_ic_peak, ion_peak_list)
            self._ion_peaks.append(new_ic_peak)
            self.ion_areas[target_mass] = new_ic_peak.area

        self._combine_peaks(self, self._ion_peaks)

    def _combine_peaks(self, new_ic_peak: AbstractPeak, peak_list: list):
        new_ic_peak._area = numpy.mean([x.area for x in peak_list])
        new_ic_peak.left = numpy.mean([x.left for x in peak_list])
        new_ic_peak.right = numpy.mean([x.right for x in peak_list])
        new_ic_peak.l_share = numpy.count_nonzero([x.l_share for x in peak_list]) > len(peak_list)//2
        new_ic_peak.r_share = numpy.count_nonzero([x.r_share for x in peak_list]) > len(peak_list)//2

    def crop_mass(self, mass_min: float, mass_max: float):
        pass


class Peak(AbstractPeak):
    """
    Subclass of :class:`~.Peak` representing a peak in a mass spectrum.

    :param rt: Retention time.
    :param ms: The mass spectrum at the apex of the peak.
    :param minutes: Retention time units flag. If :py:obj:`True`, retention time
        is in minutes; if :py:obj:`False` retention time is in seconds.
    :param outlier: Whether the peak is an outlier.

    :authors: Vladimir Likic, Andrew Isaac,
        Dominic Davis-Foster (type assertions and properties), David Kainer (outlier flag)

    .. versionchanged:: 2.3.0

        Functionality related to single ion peaks has moved to the :class:`~.ICPeak` class.
        The two classes share a common base class, :class:`AbstractPeak`, which can be used
        in type checks for functions that accept either type of peak.

    .. versionchanged:: 2.3.0

        If the ``ms`` argument is unset an empty mass spectrum is used,
        rather than :py:obj:`None` in previous versions.

    .. TODO:: Change type hint of ``ms`` to Optional[MassSpectrum] once __new__ removed.
    """

    _mass_spectrum: MassSpectrum

    def __init__(
        self,
        rt: Union[int, float] = 0.0,
        ms: Union[float, MassSpectrum, None] = None,
        im: BaseIntensityMatrix = None,
        minutes: bool = False,
        outlier: bool = False,
    ):

        if ms is None:
            ms = MassSpectrum([], [])
        elif not isinstance(ms, MassSpectrum):
            raise TypeError("'ms' must be a MassSpectrum object")

        self._mass_spectrum = ms
        super().__init__(rt, minutes, outlier)
        # self._im = im  space hog when pickled
        self._area = None
        self._ion_peaks = None  # Not initialised, and empty list is used for peaks with no ions
        self.ion_areas = {}

        if im is not None:
            self._setup_ions(im)
        else:
            pass

    def __eq__(self, other) -> bool:
        """
        Return whether this Peak object is equal to another object.

        :param other: The other object to test equality with.
        """

        if isinstance(other, self.__class__):
            return (
                self.UID == other.UID and self.bounds == other.bounds and self.rt == other.rt and
                self.mass_spectrum == other.mass_spectrum and self.area == other.area
            )

        return NotImplemented

    @property
    def area(self):
        return self._area

    def get_ion_area(self, ion: float) -> Optional[float]:
        """
        Returns the area of a single ion chromatogram under the peak.

        :param ion: The ion to calculate the area for.

        :return: The area of the ion under this peak.
        """

        try:
            return self.ion_areas[ion]
        except KeyError:
            return None

#     @property
#     def ion_areas(self) -> Dict:
#         """
#         Returns a copy of the ion areas dict.
#
#         :return: The dictionary of ``ion: ion area`` pairs
#         """  # noqa: D400
#
#         if len(self._ion_areas) == 0:
#             raise ValueError("no ion areas set")
#
#         return copy.deepcopy(self._ion_areas)
#
#     @ion_areas.setter
#     def ion_areas(self, value: Dict):
#         raise ValueError
#         """
#         Sets the ``ion: ion area`` pairs dictionary.
#
#         :param value: The dictionary of ion:ion_area pairs
#         """
#
#         if not isinstance(value, dict) or not is_number(list(value.keys())[0]):
#             raise TypeError("'Peak.ion_areas' must be a dictionary of ion:ion_area pairs")
#
#         self._ion_areas = value

    def crop_mass(self, mass_min: float, mass_max: float):
        """
        Crops mass spectrum.

        :param mass_min: Minimum mass value.
        :param mass_max: Maximum mass value.

        :author: Andrew Isaac
        """
        if not self._mass_spectrum:
            raise ValueError("Mass spectrum is unset.")
        if not is_number(mass_min) or not is_number(mass_max):
            raise TypeError("'mass_min' and 'mass_max' must be numbers")
        if mass_min >= mass_max:
            raise ValueError("'mass_min' must be less than 'mass_max'")

        mass_list = self._mass_spectrum.mass_list

        if mass_min < min(mass_list):
            raise ValueError(f"'mass_min' is less than the smallest mass: {min(mass_list)}")
        if mass_max > max(mass_list):
            raise ValueError(f"'mass_max' is greater than the largest mass: {max(mass_list)}")

        # pre build mass_list and list of indices
        mass_mapping = numpy.logical_and(mass_min <= numpy.array(mass_list), numpy.array(mass_list) <= mass_max)

        self._mass_spectrum.mass_list = numpy.array(mass_list)[mass_mapping]
        self._mass_spectrum.mass_spec = numpy.array(self._mass_spectrum.mass_spec)[mass_mapping]

        if len(self._mass_spectrum.mass_list) == 0:
            raise ValueError("mass spectrum is now empty")
        elif len(self._mass_spectrum.mass_list) < 10:
            warn("peak mass spectrum contains < 10 points", Warning)

        # update UID
        self._make_UID()

    def get_int_of_ion(self, ion: int) -> float:
        """
        Returns the intensity of a given ion in this peak.

        :param ion: The m/z value of the ion of interest
        """

        try:
            index = self._mass_spectrum.mass_list.index(ion)
            return self._mass_spectrum.mass_spec[index]
        except (ValueError, IndexError):
            raise IndexError(
                f"'ion' out of range of mass spectrum (range "
                f"{min(self._mass_spectrum.mass_list)} to "
                f"{max(self._mass_spectrum.mass_list)})"
            )

    def get_third_highest_mz(self) -> int:
        """
        Returns the *m/z* value with the third highest intensity.
        """

        if not self._mass_spectrum:
            raise ValueError("Mass spectrum is unset.")

        mass_list = self._mass_spectrum.mass_list
        mass_spec = self._mass_spectrum.mass_spec
        # find top two masses
        best = 0
        best_ii = 0
        best2_ii = 0
        best3_ii = 0
        for ii, intensity in enumerate(mass_spec):
            if intensity > best:
                best = intensity
                best3_ii = best2_ii
                best2_ii = best_ii
                best_ii = ii

        return int(mass_list[best3_ii])

    def _make_UID(self) -> None:
        """
        Create a unique peak ID (UID):

        - Integer masses of top two intensities and their ratio (as ``Mass1-Mass2-Ratio*100``); or

        :author: Andrew Isaac
        """  # noqa: D400

        if self._mass_spectrum:
            try:
                mass_list = self._mass_spectrum.mass_list
                mass_spec = self._mass_spectrum.mass_spec
                mass_sort_indicies = numpy.argsort(self._mass_spectrum.mass_spec)[::-1]
                best_ii = mass_spec[mass_sort_indicies[0]]
                best2_ii = mass_spec[mass_sort_indicies[1]]
                try:
                    ratio = int(100 * best2_ii / best_ii)
                except ValueError as _e:
                    ratio = -1  # Expected to be 0/0 (or a NaN?)
                self._UID = f"{int(mass_list[mass_sort_indicies[0]]):d}-{int(mass_list[mass_sort_indicies[1]]):d}-{ratio:d}-{self._rt:.2f}"
            except IndexError as _e:
                super()._make_UID()
        else:
            super()._make_UID()

    @property
    def mass_spectrum(self) -> MassSpectrum:
        """
        The mass spectrum at the apex of the peak.
        """
        return copy.copy(self._mass_spectrum)

    @mass_spectrum.setter
    def mass_spectrum(self, value: MassSpectrum):
        """
        Sets the mass spectrum for the apex of the peak.
        """
        if not isinstance(value, MassSpectrum):
            raise TypeError("'Peak.mass_spectrum' must be a MassSpectrum object")
        self._mass_spectrum = value

    def null_mass(self, mass: float):
        """
        Ignore given mass in spectra.
        :param mass: Mass value to remove
        :author: Andrew Isaac
        """

        if not self._mass_spectrum:
            raise ValueError("Mass spectrum is unset.")

        if not is_number(mass):
            raise TypeError("'mass' must be a number")

        mass_list = self._mass_spectrum.mass_list

        if mass < min(mass_list) or mass > max(mass_list):
            raise IndexError("'mass' not in mass range:", min(mass_list), "to", max(mass_list))

        best = max(mass_list)
        ix = 0
        for ii in range(len(mass_list)):
            tmp = abs(mass_list[ii] - mass)
            if tmp < best:
                best = tmp
                ix = ii

        self._mass_spectrum.mass_spec[ix] = 0

    def find_mass_spectrum(self, data: BaseIntensityMatrix, from_bounds: float = False):
        """
        .. TODO:: What does this function do?

        Sets the peak's mass spectrum from the data.

        Clears the single ion chromatogram mass.

        :param data:
        :param from_bounds: Whether to use the attribute :attr:`pyms.Peak.Class.Peak.pt_bounds`
            or to find the peak apex from the peak retention time.
        """  # noqa: D400

        if not isinstance(data, BaseIntensityMatrix):
            raise TypeError("'data' must be an IntensityMatrix")

        if from_bounds:
            if self._pt_bounds is None:
                raise NameError("pt_bounds not set for this peak")
            else:
                pt_apex = self._pt_bounds[1]
        else:
            # get the index of peak apex from peak retention time
            pt_apex = data.get_index_at_time(self._rt)

        # set the mass spectrum
        self._mass_spectrum = data.get_ms_at_index(pt_apex)

    def _top_ions(self, num_ions: int = 5) -> List["ICPeak"]:
        """
        Computes the highest #num_ions intensity ions.

        :param num_ions: The number of ions to be recorded.

        :return: A list of the ions with the highest intensity.

        :authors: Sean O'Callaghan, Dominic Davis-Foster (type assertions)
        """
        if not self._mass_spectrum:
            raise ValueError("Mass spectrum is unset.")
        if not isinstance(num_ions, int):
            raise TypeError("'n_top_ions' must be an integer")

        sorted_ic = sorted(self._ion_peaks, key=lambda x: x.area)
        top_ic = sorted_ic[-num_ions:]

        return top_ic

    def _setup_ions(self, im: BaseIntensityMatrix, max_bound: int = 0):
        """
        """
        if self._ion_peaks is not None:
            raise RuntimeError(f"Ion setup changed for peak {self}")
        if im is None:
            raise RuntimeError(f"Cannot calculate ion values for composite peak")
        self._ion_peaks = []
        mat = im.intensity_array
        ms = self.mass_spectrum

        if ms is None:
            raise ValueError("The peak has no mass spectrum.")

        rt = self.rt
        apex = im.get_index_at_time(rt)

        # get peak masses with non-zero intensity
        mass_ii = [ii for ii in range(len(ms.mass_list)) if ms.mass_spec[ii] > 0]

        # get stats on boundaries
        for current_mass_index in mass_ii:
            # get ion chromatogram as list
            ia = mat[:, current_mass_index]  # ion_chrom = im.get_ic_internal(ion)
            actual_mass = ms.mass_list[current_mass_index]
            new_ion_peak = ICPeak(rt, actual_mass, minutes=False, outlier=self.is_outlier)
            new_ion_peak.setup_ion_area(ia, apex, max_bound)
            self._ion_peaks.append(new_ion_peak)

        self._area = max(sum([x.area for x in self._ion_peaks]), 0.0)

    def peak_pt_bounds(self) -> Tuple[int, int]:
        """
        Approximate the peak bounds (left and right offsets from apex).
        Calculated with the 95th percentile of all ion bounds

        :return: Left and right bounds of most ions in this peak

        :authors: Andrew Isaac, Sean O'Callaghan, Dominic Davis-Foster
        """
        # get peak masses with non-zero intensity
        left_list = []
        right_list = []

        # get stats on boundaries
        for ion_peak in self._ion_peaks:
            left_list.append(ion_peak.left)
            right_list.append(ion_peak.right)

        left_list.sort()
        right_list.sort()

        return int(math.ceil(numpy.percentile(left_list, 95))), int(math.ceil(numpy.percentile(right_list, 95)))

    def peak_top_ion_areas(
            self,
            n_top_ions: int = 5,
            max_bound: int = 0,
    ):
        """
        Calculate and return the ion areas of the five most abundant ions in the peak.

        :param n_top_ions: Number of top ions to return areas for.
        :param max_bound: Optional value to limit size of detected bound.

        :return: Dictionary of ``ion : ion_area pairs``.

        :authors: Sean O'Callaghan,  Dominic Davis-Foster (type assertions)
        """
        if not isinstance(n_top_ions, int):
            raise TypeError("'n_top_ions' must be an integer")

        if not isinstance(max_bound, int):
            raise TypeError("'max_bound' must be an integer")

        ion_areas = {}  # Dictionary to store ion:ion_area pairs
        top_ions = self._top_ions(n_top_ions)
        for ion in top_ions:
            ion_areas[ion.mass] = ion.area

        self.ion_areas = ion_areas

    def median_bounds(self, shared: bool = True) -> Tuple[float, float]:
        """
        Calculates the median of the left and right bounds found for each apexing peak mass.

        :param im: The originating IntensityMatrix object.
        :param peak:
        :param shared: Include shared ions shared with neighbouring peak.

        :return: Median left and right boundary offset in points.

        :authors: Andrew Isaac, Dominic Davis-Foster
        """
        if not isinstance(shared, bool):
            raise TypeError("'shared' must be a boolean")

#         # check if RT based index is similar to stored index
#         if is_sequence(peak.bounds):
#             bounds = cast(Sequence, peak.bounds)
#             if apex - 1 < bounds[1] < apex + 1:
#                 apex = bounds[1]

        # get stats on boundaries
        left_list = []
        right_list = []

        for ion_peak in self._ion_peaks:
            if shared or not ion_peak.l_share:
                left_list.append(ion_peak.left)
            if shared or not ion_peak.r_share:
                right_list.append(ion_peak.right)

        # return medians
        # NB if shared=True, lists maybe empty
        l_med = 0.0
        r_med = 0.0
        if len(left_list) > 0:
            l_med = numpy.median(left_list)
        if len(right_list) > 0:
            r_med = numpy.median(right_list)

        return l_med, r_med


class ICPeak(AbstractPeak):
    """
    Subclass of :class:`~.Peak` representing a peak in an ion chromatogram for a single mass.

    :param rt: Retention time.
    :param mass: The mass of the ion.
    :param minutes: Retention time units flag. If :py:obj:`True`, retention time
        is in minutes; if :py:obj:`False` retention time is in seconds.
    :param outlier: Whether the peak is an outlier.

    :authors: Vladimir Likic, Andrew Isaac,
        Dominic Davis-Foster (type assertions and properties), David Kainer (outlier flag)

    .. versionadded:: 2.3.0
    """

    def __init__(
            self,
            rt: Union[int, float] = 0.0,
            mass: float = None,
            minutes: bool = False,
            outlier: bool = False,
            bounds: Optional[list] = [],
    ):

        if mass and not is_number(mass):
            raise TypeError("'ms' must be a number")

        self.mass = mass
        self._bounds = bounds

        super().__init__(rt, minutes, outlier)
        self._make_UID()
        # Requires setup_ion_area(ia, apex) to populate

    def __eq__(self, other) -> bool:
        """
        Return whether this Peak object is equal to another object.

        :param other: The other object to test equality with.
        """

        if isinstance(other, self.__class__):
            return (
                self.UID == other.UID and self.bounds == other.bounds and self.rt == other.rt and
                self.area == other.area and self.mass == other.mass
            )
        return NotImplemented

    def _make_UID(self) -> None:
        """
        Create a unique peak ID (UID):

        - the single mass as an integer and the retention time.

        :author: Andrew Isaac
        """  # noqa: D400

        if self.mass is not None:
            self._UID = f"{int(self.mass):d}-{self._rt:.2f}"
        else:
            super()._make_UID()

    def setup_ion_area(
            self,
            ia: List,
            apex: int,
            max_bound: int = 0,
            tol: float = 0.5,
    ) -> Tuple[float, float, float, float, float]:
        """
        Find bounds of peak by summing intensities until change in sum is less than
        ``tol`` percent of the current area.
    
        :param ia: List of intensities for a given mass.
        :param apex: Index of the peak apex.
        :param max_bound: Optional value to limit size of detected bound.
        :param tol: Percentage tolerance of added area to current area.
    
        :return: Area, left and right boundary offset, shared left, shared right.
    
        :authors: Andrew Isaac, Dominic Davis-Foster (type assertions)
        """  # noqa: D400

        # Left area
        lhs = ia[:apex + 1]
        # lhs.reverse()  # reverse, as search to right is bounds safe
        l_area, left, l_share = self._half_area(lhs[::-1], max_bound, tol)

        # Right area
        rhs = ia[apex:]
        r_area, right, r_share = self._half_area(rhs, max_bound, tol)
        r_area -= ia[apex]  # Counted apex twice for tolerance now ignore

        # Put it all together
        self._area = l_area + r_area
        self.left = left
        self.right = right
        self.l_share = l_share
        self.r_share = r_share

    def _half_area(
            self,
            ia: List,
            max_bound: int = 0,
            tol: float = 0.5,
    ) -> Tuple[float, float, float]:
        """
        Find bound of peak by summing intensities until change in sum is less than
        ``tol`` percent of the current area.
    
        :param ia: List of intensities from Peak apex for a given mass.
        :param max_bound: Optional value to limit size of detected bound.
        :param tol: Percentage tolerance of added area to current area.
    
        :return: Half peak area, boundary offset, shared (True if shared ion).
    
        :authors: Andrew Isaac, Dominic Davis-Foster (type assertions)
        """  # noqa: D400

        tol = tol / 200.0  # halve and convert from percent

        # Default number of points to sum new area across, for smoothing
        wide = 3

        # start at 0, compare average value of 'wide' points to the right,
        # centre 'wide' points on edge point, and keep moving right until:
        # i) tolerance reached
        # ii) edge area starts increasing
        # iii) bound reached

        # initialise areas and bounds
        shared = False
        area = ia[0]
        edge = float(sum(ia[0:wide])) / wide
        old_edge = 2 * edge  # bigger than expected edge
        index = 1
        if max_bound < 1:
            limit = len(ia)
        else:
            limit = min(max_bound + 1, len(ia))
        # while edge > area * tol and edge < old_edge and index < limit:
        while area * tol < edge < old_edge and index < limit:
            old_edge = edge
            area += ia[index]
            edge = float(sum(ia[index:index + wide])) / wide  # bounds safe
            index += 1
        if edge >= old_edge:
            shared = True

        return area, index, shared
