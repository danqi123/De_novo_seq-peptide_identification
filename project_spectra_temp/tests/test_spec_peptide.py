"""Unit test for Spec_peptide."""

from project_spectra.Spec_peptide import Spectrum_processing, Peptide_identification
from project_spectra.startup import data_path
import os
import pandas as pd
from collections import Counter

class TestSpecPeptide:
    """Test code for project_spectra.Spec_peptide."""

    # global variable for testing
    sp = Spectrum_processing(ms2_file='data/test_ms2.mzML', ms2_scan_number=0)
    pi = Peptide_identification(ms2_file='data/test_ms2.mzML', scan_number=0)

    def test_store_one_scan_and_get_deisotoped(self):
        """Test function store_one_scan_and_get_deisotoped"""
        self.sp.store_one_scan_and_get_deisotoped()
        assert os.path.exists(data_path+'/0_deisotoped.mzML')

    def test_normalization(self):
        """Test function normalization"""
        dic = {'test': 100}
        result = self.sp.normalization(dic)
        assert result['test'] == 100

    def test_parse_dict(self):
        """Test function parse_dict"""
        global dic
        dic = self.sp.parse_dict()
        assert isinstance(dic, dict)

    def test_precursor_mz(self):
        """Test function precursor_mz"""
        assert isinstance(self.sp.precursor_mz(), float)

    def test_precursor_charge(self):
        """Test function precursor_charge"""
        assert isinstance(self.sp.precursor_charge(), int)

    def test_precursor_ms(self):
        """Test function precursor_ms"""
        global mass
        mass = self.sp.precursor_ms()
        assert isinstance(mass, float)

    def test_plot_spectrum(self):
        """Test function plot_spectrum"""
        self.sp.plot_spectrum(False)
        assert os.path.exists(data_path+"/0_deisotoped.jpg")

    def test_parse_scan_data(self):
        """Test function parse_scan_data"""
        df, mz_dict = self.sp.parse_scan_data(data_path+"/0_deisotoped.mzML")
        assert isinstance(df, pd.DataFrame)
        assert isinstance(mz_dict, dict)

    def test_identify_C_term(self):
        """Test function identify_C_term"""
        global norm_dic
        norm_dic = self.sp.normalization(dic)
        result = self.pi.identify_C_term(precursor_ms=mass, dict_norm=norm_dic)
        assert result == True

    def test_compile(self):
        """Test function compile. Execute this first to generate attributes in class."""
        assert isinstance(self.pi.compile(False), list)

    def test_get_up_downstream_peak_dict(self):
        """Test function get_up_downstream_peak_dict"""
        result = self.pi.get_up_downstream_peak_dict(dict_norm=norm_dic)
        assert isinstance(result, dict)
        assert mass in result.keys()

    def test_get_peak_pair(self):
        """Test function get_peak_pair"""
        global pair_list
        pair_list = self.pi.get_peak_pair()
        assert isinstance(pair_list, list)
        for elem in pair_list:
            assert len(elem) == 2
            assert elem[0] > elem[1]

    def test_end_to_end_concate(self):
        """Test function end_to_end_concate"""
        global merged_list
        merged_list = self.pi.end_to_end_concate(input_list=pair_list)
        assert isinstance(merged_list, list)
        for elem in merged_list:
            n = len(elem)
            for i in range(n-1):
                assert elem[i] > elem[i+1]

    def test_get_peak_list(self):
        """Test function get_peak_list"""
        global new_list
        global merged_list
        new_list = self.pi.get_peak_list(merge_list=merged_list)
        assert isinstance(new_list, list)
        for i in range(3):
            assert len(merged_list[i]) + 1 == len(new_list[i])

    def test_get_peptide_list(self):
        """Test function get_peptide_list"""
        seq_list = self.pi.get_peptide_list()
        for elem in seq_list:
            assert isinstance(elem, list)

