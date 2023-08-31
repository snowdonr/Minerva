'''
Created on Nov. 2, 2021
From pymassspec documentation

@author: Ryan
'''
import logging
import operator
import time

import pyms.Noise.SavitzkyGolay
import pyms.Noise.Analysis
import pyms.BillerBiemann
import pyms.TopHat
import sqlalchemy.orm
from sqlalchemy import Column, Integer, Float, ForeignKey

from . import sql_base
from . import mass_spectrum


class Parse_IM():
    ''' Wrap and process the Intensity matrix from pymassspec '''

    def __init__(self, input_sample: 'sample.Sample'):
        self.sample = input_sample
        self.config = input_sample.config

    def load_from_im(self, intensity_matrix: pyms.IntensityMatrix.IntensityMatrix):
        '''
        Apply pyms functions to extract peaks, see pymassspec documentation
        '''
        n_scan, n_mz = intensity_matrix.size
        print(f"File has {n_scan} scans and {n_mz} masses")

        # perform necessary pre filtering
        status_update_count = n_mz / 10
        for ii in range(n_mz):
            if ii % status_update_count<1:
                print(f"Working on IC#{ii+1}/{n_mz}")
            ic = intensity_matrix.get_ic_at_index(ii)
            ic_smooth = pyms.Noise.SavitzkyGolay.savitzky_golay(ic)
            ic_bc = pyms.TopHat.tophat(ic_smooth, struct=self.config.tophat_baseline_noise_smooth)  # Noise removal (smoothing) filter for baseline correction
            intensity_matrix.set_ic_at_index(ii, ic_bc)
        print(f"STATUS: Done IC scan, Detecting peaks")

        # Detect Peaks
        peak_list = pyms.BillerBiemann.BillerBiemann(intensity_matrix,
                                                     points=self.config.bb_peak_points,
                                                     scans=self.config.bb_scan_merge)  # Can resolve retention time peaks <0.5scan based on all mass fragments
        print(f"Number of detected peaks: {len(peak_list)}")

        # ######## Filter peaks###############
        # Filter the peak list,
        # first by removing all intensities in a peak less than a given relative threshold,
        # then by removing all peaks that have less than a given number of ions above a given value
        r = self.config.peak_trim_percent_of_max  # % -> percentage ratio of ion intensity to max ion intensity
        peak_list_r = pyms.BillerBiemann.rel_threshold(peak_list, r)  # trim by relative intensity
        print(f"Peaks above {r}% of max: {len(peak_list_r)}")

        noise_level = pyms.Noise.Analysis.window_analyzer(self.sample.get_tic())
        pyms_peak_list = pyms.BillerBiemann.num_ions_threshold(peak_list_r, self.config.min_ions, cutoff=noise_level)  # trim by threshold
        print(f"Remaining with {self.config.min_ions} ions at >{noise_level:.3g} intensity: {len(pyms_peak_list)}")

        self.sample.fragment_count = len(pyms_peak_list)
        # Set the peak areas
        print(f"STATUS: Calculating ion areas")
        try:
            status_update_count = len(pyms_peak_list) / 10
            for count, peak in enumerate(pyms_peak_list):
                if count % status_update_count<1:
                    print(f"Finished {count}/{len(pyms_peak_list)} peaks")
                peak.peak_top_ion_areas(n_top_ions=self.config.top_ion_areas)
        except Exception as _e:
            logging.exception("Could not set peak ion areas")
        print("STATUS: Done calculating peak areas")

        return pyms_peak_list

    def status_update(self, start_time, current_amount: float, total_amount: float):
        current_time = time.process_time()
        frac_complete = current_amount/total_amount
        time_remain = (1.0-frac_complete)*(start_time-current_time)
        print(f"Finished {current_amount}/{total_amount} ({100*frac_complete:.1f}) processing. Section Remaining: {time_remain/60:.1f}")


class Peak(sql_base.Base):
    __tablename__ = "GC_Peak"

    id = Column(Integer, primary_key=True)
    area = Column(Float, nullable=False)
    rt = Column(Float, nullable=False)

    aligned_id = Column(Integer, ForeignKey("Aligned_Peak_Set.id"))
    aligned_set = sqlalchemy.orm.relationship("AlignedPeakSet_SQL", back_populates="gc_peaks")
    sample_id = Column(Integer, ForeignKey("Lab_Sample_Analysis.id"))
    source_sample = sqlalchemy.orm.relationship("Sample", back_populates="gc_peak")
    mass_spectrum_id = Column(Integer, ForeignKey("Mass_Spectrum.id"))
    mass_spectrum = sqlalchemy.orm.relationship("MassSpectrum_SQL", back_populates="gc_peak")

    def __init__(self, pyms_peak):
        super().__init__()
        self.area = pyms_peak.area
        self.rt = pyms_peak.rt
        self.source_sample = pyms_peak.source_sample
        self.mass_spectrum = mass_spectrum.MassSpectrum_SQL()

    def set_ions(self, ion_dict: dict, max_ion_entries=3):
        ion_tuple_list = sorted(list(ion_dict.items()), key=operator.itemgetter(1))[::-1]
        for i, ion_tuple in enumerate(ion_tuple_list):
            new_spectrum_entry = mass_spectrum.MassSpectrumEntry(mass=int(ion_tuple[0]), area=ion_tuple[1])
            self.mass_spectrum.values.append(new_spectrum_entry)
            if i >= max_ion_entries-1:
                break
