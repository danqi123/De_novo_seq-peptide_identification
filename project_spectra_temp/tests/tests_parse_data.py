"""Unit test for Parse_data_MS2"""

import os
import pandas as pd
from project_spectra.Parse_data_MS2 import get_ms2, parse_data_from_MS2


class TestParseData:
    """Test code for project_spectra.Parse_data_MS2"""

    def test_get_ms2(self):
        """Test function get_ms2"""
        file_path = "data/test_data.mzML"
        output_path = 'data/test_ms2.mzML'
        get_ms2(file_path, output_path)
        assert os.path.exists(output_path)

    def test_parse_data_from_MS2(self):
        """Test function parse_data_from_MS2"""
        input_file = "data/test_ms2.mzML"
        scan_number = 0
        df, mz_int, pre_mz, pre_charge = parse_data_from_MS2(input_file, scan_number)
        assert isinstance(df, pd.DataFrame)
        assert isinstance(mz_int, dict)
        assert isinstance(pre_mz, float)
        assert isinstance(pre_charge, int)



