
"""This module is used to parse data from raw .mzML data
   Transform into MS2 file and get m/z and intensity from one scan
"""

from pyopenms import *
import pandas as pd

def get_ms2(ms_file: str, output_file: str) -> None:
    """
    Store MS2 file from raw LC-MS/MS mzXML file.
    Parameters
    ----------
    ms_file: str
            The .mzXML file
    Returns
    -------
    None
    """
    exp = MSExperiment()
    MzMLFile().load(ms_file, exp)
    spec = []
    for s in exp.getSpectra():
        if s.getMSLevel() == 2:
            spec.append(s)
    exp.setSpectra(spec)
    MzMLFile().store(output_file, exp)
    return


def parse_data_from_MS2(input_file: str, ms2_scan_number: int) -> tuple:
    """Method used to parse data from one specific scan.
       RAW DATA."""
    exp = MSExperiment()
    MzMLFile().load(input_file, exp)
    mz = []
    intensity = []
    for peak in exp[ms2_scan_number]:
        mz.append(peak.getMZ())
        intensity.append(peak.getIntensity())

    # raw MS2 dataFrame (m/z and intensity):
    df = pd.DataFrame(columns=['M/Z', 'Intensity'])
    df['M/Z'] = mz
    df['Intensity'] = intensity

    mz_int = {x: y for x, y in zip(mz, intensity)}

    pre_mz = exp[ms2_scan_number].getPrecursors()[0].getMZ()
    pre_charge = exp[ms2_scan_number].getPrecursors()[0].getCharge()
    return df, mz_int, pre_mz, pre_charge

