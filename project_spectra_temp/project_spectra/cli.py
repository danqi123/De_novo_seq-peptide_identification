#!/usr/bin/env python
"""CLI module"""

import click
from pyopenms import *
from project_spectra.Parse_data_MS2 import get_ms2, parse_data_from_MS2
from project_spectra.Spec_peptide import Peptide_identification
from project_spectra.identifier import Identifier


@click.group(help = f'The Command Line Utilities of generating peptides and candidate proteins.')
def main():
	"""Entry method."""
	pass

# The first CLI for generate MS2
@main.command(name = "get_ms2")
@click.argument('input_path')
@click.argument('output_path')
def ms2_cli(input_path: str, output_path: str) -> None:

	"""Generates MS2 .ML file from raw data."""
	print("File extracting...")
	get_ms2(input_path, output_path)
	return

# The second CLI for parsing
@main.command(name = "check_scan_size")
@click.argument('input_path')
def check_scan_size(input_path: str) -> int:

	"""Get the total scan number of MS2 file."""
	exp = MSExperiment()
	MzMLFile().load(input_path, exp)
	size = len(exp.getSpectra())
	print(size)



# The third CLI for parsing
@main.command(name = "parse")
@click.argument('ms2_path')
@click.argument('scan_number')
def parse_cli(ms2_path: str, scan_number: int) -> None:

	"""Parse M/Z, intensity, precursor info of specific scan from MS2 file."""

	info = parse_data_from_MS2(ms2_path, int(scan_number))
	print("Here comes the information of PRECURSOR...:")
	print(f"The M/Z is: {info[2]}")
	print(f"The Charge is: {info[3]}")
	print("------------------------------------------")
	print("The M/Z and intensity info of all fragments from this precursor:")
	print(info[0])


# The fourth CLI for generating plot
@main.command(name = 'plot')
@click.argument('ms2_raw')#, help="the MS2 raw file"
@click.argument('scan_number')#, help="The specific scan number."
@click.option('-v', '--verbose', default=False, is_flag=True, help="When used, will print the spectrum to STDOUT.")
def get_plot_spectrum(ms2_raw: str, scan_number: int, verbose: bool):

	"""Generate Plot of de-isotoped spectrum of one MS2 scan."""
	p = Peptide_identification(ms2_raw, int(scan_number))
	p.store_one_scan_and_get_deisotoped()
	p.plot_spectrum(print_option = verbose)


# The fifth CLI for checking C terminal
@main.command(name = 'check_C')
@click.argument('ms2_raw')
@click.argument('scan_number')
def identify_C_terminal(ms2_raw: str, scan_number: int) -> bool:

	"""check C terminal (b-ion) R or K is identified or not."""
	p = Peptide_identification(ms2_raw, int(scan_number))
	p.AA_code()
	p.store_one_scan_and_get_deisotoped()
	dict_p = p.parse_scan_data(p.raw_file)[1]
	p.normalization(dict_p)
	p.precursor_ms()
	print("Can we identify b-ion at C-terminal with R or K....?")
	print(p.identify_C_term(p.mass, p.norm_dict))


# The sixth CLI for checking N terminal
@main.command(name = 'check_N')
@click.argument('ms2_raw')
@click.argument('scan_number')
def identify_N_terminal(ms2_raw: str, scan_number: int) -> bool:

	"""check N terminal (y-ion) is identified or not."""
	p = Peptide_identification(ms2_raw, int(scan_number))
	p.AA_code()
	p.store_one_scan_and_get_deisotoped()
	dict_p = p.parse_scan_data(p.raw_file)[1]
	p.normalization(dict_p)
	p.precursor_ms()
	print("Can we identify the first y-ion at N-terminal using our algorithm....?")
	print(p.have_N_term())


# The seventh CLI for generating peptide list
@main.command(name = 'peptide_list')
@click.argument('ms2_raw')
@click.argument('scan_number')
@click.option('-v', '--verbose', default=False, is_flag=True, help="When used, will print the C_terminal to STDOUT.")
def get_peptide(ms2_raw: str, scan_number: int, verbose:bool) -> list:

	"""get peptide seq list from one scan and the detailed info about the peptide seq"""
	import pandas as pd

	p = Peptide_identification(ms2_raw, int(scan_number))
	try:
		list_peptide = p.compile(show_C_terminal=verbose)[0]
		print("Here comes the peptide lists....: ")
		# show peptide list in STDOUT
		print(list_peptide)
		print("-----------------------------------")
		seq_info = [p.peptide_info(x) for x in list(set(p.compile(show_C_terminal=verbose)[0]))]
		seq = p.filter_seq
		sum_intensity = p.filter_peak_intensity
		lst = zip(seq, sum_intensity)
		result_p = pd.DataFrame(lst, columns=["Sequence", "Sum of Peak Intensity/ Score"])
		result_p = result_p.drop_duplicates()
		# show table of sequence and its total peak intensity in STDOUT
		print(result_p)
		# return the details info of peptide
		return seq_info
	except:
		p.precursor_charge()
		p.precursor_ms()
		if not p.identify_C_term(p.mass, p.norm_dict) or not p.have_N_term():
			print("N or C terminal not found by the algorithm.")
		else:
			print(p.merge_info)


# The eighth CLI for identifying candidate protein
@main.command(name = "protein")
@click.argument('ms2_raw')
@click.argument('scan_number')
@click.option('-v', '--verbose', default=False, is_flag=True, help="When used, will print the C_terminal to STDOUT.")
def get_protein(ms2_raw: str, scan_number: int, verbose:bool) -> None:
	""" Search for candidate proteins from API request.

        Parameters
        ----------
        ms2_raw: str
        scan_number: int
        verbose: bool

        Returns
        -------
        STDOUT: Peptide list
                Protein info: Accession number, Full name and Organism info
	"""
	import pandas as pd
	p = Peptide_identification(ms2_raw, int(scan_number))
	list_ = list(set(p.compile(show_C_terminal=verbose)[0]))
	print("--------------------------------------------------------------------------")
	print("Peptides combination with their corresponding protein info.....")
	print("###############################################################")
	peptide_match = Identifier(sequence=list_)
	# Print list of proteins
	print(f"{list} candidate protein info: {pd.DataFrame(peptide_match.stout.items())}")


if __name__ == "__main__":
	main()
