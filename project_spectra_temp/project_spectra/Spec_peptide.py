
""" This module is used to find peptides.
    The output is a list contains possible peptides.
"""

import pandas as pd
import matplotlib.pyplot as plt
from pyopenms import *
from project_spectra.Parse_data_MS2 import parse_data_from_MS2
from project_spectra.startup import data_path
from collections import defaultdict
import itertools


class Spectrum_processing():
    def __init__(self, ms2_file: str, ms2_scan_number: int):
        self.ms2_file = ms2_file
        self.s_number = ms2_scan_number
        self.raw_file = data_path + f"/{self.s_number}_full.mzML"
        self.deisotoped_file = data_path + f"/{self.s_number}_deisotoped.mzML"
        self.plot_file = data_path + f"/{self.s_number}_deisotoped.jpg"

    def store_one_scan_and_get_deisotoped(self) -> None:
        """This method is used to deisotope peaks for plotting and store specific (scan number) raw spectrum to .ML file"""
        min_isotopes = 2
        max_isotopes = 10
        use_decreasing_model = True
        start_intensity_check = 3

        e = MSExperiment()
        MzMLFile().load(self.ms2_file, e)
        s = e[self.s_number]
        s.setFloatDataArrays([])
        Deisotoper.deisotopeAndSingleCharge(s, 0.1, False, 1, 3, True,
                                            min_isotopes, max_isotopes,
                                            True, True, True,
                                            use_decreasing_model, start_intensity_check, False)

        # use two variables to store the size of each scan
        self.full_scan_size = e[self.s_number].size()
        self.de_scan_size = s.size()

        # store specific scan to new .ML files.
        e2 = MSExperiment()
        e2.addSpectrum(e[self.s_number])
        MzMLFile().store(self.raw_file, e2)

        e2 = MSExperiment()
        e2.addSpectrum(s)
        MzMLFile().store(self.deisotoped_file, e2)
        return

    def normalization(self, peak_dict: dict) -> dict:
        """Normalize the intensity as relative abundance"""
        max_value = max(list(peak_dict.values()))
        self.norm_dict = {x: (y/max_value)*100 for x, y in peak_dict.items()}
        return self.norm_dict

    def parse_dict(self) -> dict:
        """get dictionary: m/z is key, and intensity is value."""
        return parse_data_from_MS2(self.ms2_file, self.s_number)[1]

    def precursor_mz(self) -> float:
        """get m/z of precursor"""
        return parse_data_from_MS2(self.ms2_file, self.s_number)[2]

    def precursor_charge(self) -> int:
        """get charge of precursor"""
        return parse_data_from_MS2(self.ms2_file, self.s_number)[3]

    def precursor_ms(self) -> float:
        """Based on the m/z and charge of precursor, calculate the mass of precursor"""
        m_z = self.precursor_mz()
        charge = self.precursor_charge()
        self.mass = m_z * charge
        return self.mass

    def plot_spectrum(self, print_option: False) -> None:
        """ method used to generate spectrum plot of deisotoped spectrum"""
        exp = MSExperiment()
        MzMLFile().load(self.deisotoped_file, exp)
        for spec in exp:
            for mz, i in zip(*spec.get_peaks()):
                plt.plot([mz, mz], [0, i], color = 'black')
                plt.text(mz, i, str(mz))

            # for the title add RT and Precursor m/z and charge info if available
            title = ''
            if spec.getRT() >= 0:
                title += 'RT: ' + str(spec.getRT())
            if len(spec.getPrecursors()) >= 1:
                title += '   Precursor m/z: ' + str(spec.getPrecursors()[0].getMZ()) + '   Charge: '+ str(spec.getPrecursors()[0].getCharge())
            plt.title(title)
            plt.ylabel('intensity')
            plt.xlabel('m/z')
            plt.ylim(bottom=0)
            plt.savefig(self.plot_file)

            if print_option:
                plt.show()

    @staticmethod
    def parse_scan_data(specific_scan_MLfile: str) -> tuple:
        """parse specific scan into dictionary and dataframe"""
        exp = MSExperiment()
        MzMLFile().load(specific_scan_MLfile, exp)
        for spec in exp:
            mz, intensity = spec.get_peaks()
        mz_int_dict = {x: y for x, y in zip(mz, intensity)}

        df = pd.DataFrame(columns=['M/Z', 'Intensity'])
        df['M/Z'] = mz
        df['Intensity'] = intensity
        return df, mz_int_dict


class Peptide_identification(Spectrum_processing):
    def __init__(self, ms2_file: str, scan_number: int):
        super().__init__(ms2_file, scan_number)
        self.aa_code = self.AA_code()

    def AA_code(self) -> dict:
        """Amino Acid -> Mass dict"""
        aa_code = {
        'A': 71, 'C': 103, 'D': 115, 'E': 129, 'F': 147,
        'G': 57, 'H': 137, 'I': 113, 'K': 128, 'L': 113,
        'M': 131, 'N': 114, 'P': 97, 'Q': 128, 'R': 156,
        'S': 87, 'T': 101, 'V': 99, 'W': 186, 'Y': 163,
        }
        self.aa = list(aa_code.keys())
        self.aa_mass = list(aa_code.values())
        return aa_code

    def identify_C_term(self, precursor_ms: float, dict_norm: dict, tolerance_error = 1) -> bool:
        """This function is used to get the C-terminal residue.(b-ion)
        It must be R or K, since Trypsin is used for digestion.
        If it has R or K: return True,
        else: False
        """
        mass_set = []
        m_z = list(dict_norm.keys())
        b_ion = []

        for elem in m_z:
            C_term_1 = int(precursor_ms - elem - 18)
            C_term_2 = int(precursor_ms - elem - 18 - 17)
            C_term_3 = int(precursor_ms - elem - 18 - 18)
            if self.aa_code["R"] in range(C_term_1 - tolerance_error, C_term_1 + 1 + tolerance_error):
                # only select the peaks with higher intensity as candidate
                b_ion.append("R")
                mass_set.append(elem)
            elif self.aa_code["K"] in range(C_term_1 - tolerance_error, C_term_1 + 1 + tolerance_error):
                b_ion.append("K")
                mass_set.append(elem)
            elif self.aa_code["R"] in range(C_term_2 - tolerance_error, C_term_2 + 1 + tolerance_error):
                b_ion.append("R")
                mass_set.append(elem)
            elif self.aa_code["K"] in range(C_term_2 - tolerance_error, C_term_2 + 1 + tolerance_error):
                b_ion.append("K")
                mass_set.append(elem)
            elif self.aa_code["R"] in range(C_term_3 - tolerance_error, C_term_3 + 1 + tolerance_error):
                b_ion.append("R")
                mass_set.append(elem)
            elif self.aa_code["K"] in range(C_term_3 - tolerance_error, C_term_3 + 1 + tolerance_error):
                b_ion.append("K")
                mass_set.append(elem)

        intensity_set = {elem: dict_norm[elem] for elem in mass_set}

        try:
            ret = max(intensity_set, key = lambda x: intensity_set[x])

            self.C_term_peak = [ret]
            self.C_term_aa = [b_ion[mass_set.index(ret)]]
            return True
        except:
            self.C_term_peak = None
            self.C_term_aa = None
            self.info_C = "cannot find b_ion at C terminal, no peptide generated here."
            return False


    def have_N_term(self) -> bool:
        """ This function is used to if we can find candidate N terminal AA.
        If true: continue for the algorithm;
        else: return False"""
        T_F_list = []
        for elem in list(self.norm_dict.keys()):
            mass_diff = int(self.mass) - int(elem)
            if mass_diff in self.aa_mass:
                T_F_list.append(True)
        if True in T_F_list:
            return True
        else:
            self.info_N = "cannot find y_ion at N terminal, no peptide generated here."
            return False


    def get_up_downstream_peak_dict(self, dict_norm: dict, min_AA = 57, max_AA = 187) -> dict:
        """
        DFS algorithm. search from N-terminal -> C terminal, largest m/z y-ion -> smallest
        Store up parent peak as key, and its corresponding child peak as list of values.
        Parameters
        ----------
        dict_norm: dict (m/z and normalized intensity)
        min_AA: the minimum of AA (which is ued as the searching window)
        max_AA: the maximum of AA

        Returns
        -------
        route_dict: dict
        """

        stack = []
        stack.append(self.mass)
        self.route_dict = defaultdict(list)
        while len(stack) > 0:
            init = stack.pop()
            searching_window = [elem for elem in list(dict_norm.keys()) if int(init)-max_AA < int(elem) < int(init)-min_AA]
            for candidate_peak in searching_window:
                mass_diff = int(init) - int(candidate_peak)
                # mass_difference matches to mass of AA
                if mass_diff in self.aa_mass:
                    # add down stream peaks here.
                    self.route_dict[init].append(candidate_peak)
                    stack.append(candidate_peak)

        return self.route_dict

    def get_peak_pair(self) -> list:
        """collection of all possible peak linkage."""
        # e.g. [[707.369183213676, 560.1217041015625], [560.1217041015625, 432.19842529296875], [560.1217041015625, 431.2114562988281]]
        up_stream_peaks = self.route_dict.keys()
        cc = list(itertools.combinations(up_stream_peaks, 2))
        self.pair_list = []
        for x, y in cc:
            if int(x) - int(y) in self.aa_mass:
                self.pair_list.append([x, y])
        return self.pair_list

    def end_to_end_concate(self, input_list: list) -> list:
        """This algorithm is used to generate list of peaks."""
        # e.g. [[707.369183213676, 560.1217041015625, 432.19842529296875, 331.0882568359375],
        #       [707.369183213676, 560.1217041015625, 432.19842529296875, 301.0992431640625]]
        l_N = len(input_list)
        sequences = []
        concated = []
        sequences.append(input_list[0])
        concated.append(False)
        for i in range(1, l_N):
            test = input_list[i]
            flag = False
            for j in range(len(sequences)):
                if sequences[j][-1] == test[0]:
                    concated[j] = True
                    sequences.append(sequences[j] + test[1:])
                    concated.append(False)
                    flag = True
            if not flag:
                sequences.append(test)
                concated.append(False)
        self.merge_list = []
        for i in range(len(sequences)):
            if not concated[i]:
                self.merge_list.append(sequences[i])
        return self.merge_list

    def get_peak_list(self, merge_list: list):
        """Add the child peak for the last y-ion parent peak."""
        self.new_list = []
        try:
            for elem in merge_list:
                for r in self.route_dict[elem[-1]]:
                    self.new_list.append(elem + [r])
            return self.new_list
        except:
            # merge_list is not available, not enough peaks there.
            self.merge_info = "Not enough AA generated, NO peptide sequence."
            return False

    def get_peptide_list(self, peptide_identification_tolerance = 1) -> list:
        """From the peak list to generate corresponding peptide lists.
           Here we should add the missing C terminal AA."""
        mass_list = []
        seq_list = []

        for elem in self.new_list:
            total_mass = 0
            seq = ""
            for i in range(0, len(elem) - 1):
                mass_difference = int(elem[i]) - int(elem[i + 1])
                seq += self.aa[self.aa_mass.index(mass_difference)]
                total_mass += mass_difference

            # Add corresponding C-terminal AA
            total_mass += self.AA_code()[self.C_term_aa[0]]
            seq += self.C_term_aa[0]
            mass_list.append(total_mass)
            seq_list.append(seq)

        self.filter_seq = []
        self.filter_peak = []
        for i, elem in enumerate(mass_list):
            if abs(elem - self.mass) < peptide_identification_tolerance:
                self.filter_seq.append(seq_list[i])
                self.filter_peak.append(self.new_list[i])
        # here we should see the total relative abundance of peak to make a score,
        # but not take the precursor (first peak) into account.
        self.filter_peak_intensity = [sum([self.norm_dict[i] for i in x[1:]]) for x in self.filter_peak]
        combine_info = [self.filter_seq, self.filter_peak, self.filter_peak_intensity]
        return combine_info

    def compile(self, show_C_terminal: False):
        """This is a compile version of getting final peptide list."""
        self.AA_code()
        self.store_one_scan_and_get_deisotoped()
        dict_p = self.parse_scan_data(self.raw_file)[1]
        self.normalization(dict_p)
        self.precursor_ms()

        # If N- and C- terminal were identified, we can continue the algorithm...
        if self.identify_C_term(self.mass, self.norm_dict):
            if self.have_N_term():
                # here we can see C terminal (whether R or k)
                if show_C_terminal:
                    print(f"The C_terminal is:{self.C_term_aa, self.C_term_peak}")
                self.get_up_downstream_peak_dict(self.norm_dict)
                self.get_peak_pair()
                merge_list = self.end_to_end_concate(self.pair_list)
                if self.get_peak_list(merge_list):
                    return self.get_peptide_list()

    @staticmethod
    def peptide_info(peptide: str):
        seq = AASequence.fromString(peptide)
        print("The peptide", str(seq), "consists of the following amino acids:")
        for aa in seq:
            print(aa.getName(), ":", aa.getMonoWeight())

