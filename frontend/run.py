# -*- coding: utf-8 -*-
"""This module contains the Spectra Flask Application application."""

import logging
import time
from pathlib import Path
import pandas as pd
from pyopenms import *
import os

from flask_sqlalchemy import SQLAlchemy
from flask import Flask, render_template, request, redirect, url_for, flash, session
from werkzeug.utils import secure_filename

from project_spectra.Parse_data_MS2 import parse_data_from_MS2
from project_spectra.Spec_peptide import Peptide_identification
from project_spectra.identifier import Identifier

log = logging.getLogger(__name__)

# Constants
ALLOWED_EXTENSIONS = {'mzml', 'fasta', 'fa'}

# Create upload folder
upload_folder = Path(Path.home(), '.projectSpectra', 'uploads')
Path.mkdir(upload_folder, exist_ok=True)
UPLOAD_FOLDER = str(upload_folder)


"""Create the Flask application"""

t = time.time()

app = Flask(__name__)

FLASK_PORT = os.environ.get('FLASK_PORT', default=5000)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

app.config['SECRET_KEY'] = "1P313P4OO138O4UQRP9343P4AQEKRFLKEQRAS230"
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///spectra.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# Initialize the database
db = SQLAlchemy(app)
log.info('Done building %s in %.2f seconds', app, time.time() - t)


"""Build URLs"""

@app.route("/")
@app.route("/start")
def start():
	return render_template("start.html")

# Home
@app.route("/home", methods=['GET', 'POST'])
def home():
	ms2_file = os.path.join(app.config['UPLOAD_FOLDER'], "ms2.mzML")
	filenames = next(os.walk(UPLOAD_FOLDER), (None, None, []))[2]  # Get files
	filenames = [file for file in filenames]
	# Check files are uploaded
	if filenames:
		exp = MSExperiment()
		MzMLFile().load(ms2_file, exp)
		size = len(exp.getSpectra())
	return render_template('home.html', size = size, my_info ="", data_table ="",
						   p_mz ='', p_c ="", number_of_peak = "-", min_m_z = "-",
							max_m_z = "-", max_m_z_peak = "-")

# Generate MS summary
@app.route('/parse', methods=['POST'])
def parse_spec():
	ms2_file = os.path.join(app.config['UPLOAD_FOLDER'], "ms2.mzML")
	filenames = next(os.walk(UPLOAD_FOLDER), (None, None, []))[2]  # Get files
	filenames = [file for file in filenames]
	# Check files are uploaded
	if filenames:
		scan_number = request.form['textbox']
		exp = MSExperiment()
		MzMLFile().load(ms2_file, exp)
		m_z = list(exp.getSpectrum(int(scan_number)).get_peaks()[0])
		intensity = list(exp.getSpectrum(int(scan_number)).get_peaks()[1])

		number_of_peak = len(m_z)
		min_m_z = round(min(m_z), 2)
		max_m_z_peak = round(m_z[intensity.index(max(intensity))], 2)
		max_m_z = round(max(m_z), 2)

		precursor_mz = parse_data_from_MS2(ms2_file, int(scan_number))[2]
		precursor_charge = parse_data_from_MS2(ms2_file, int(scan_number))[3]

		return render_template('home.html', p_mz = precursor_mz, p_c = precursor_charge,
							   number_of_peak = number_of_peak, min_m_z = min_m_z,
							   max_m_z = max_m_z, max_m_z_peak = max_m_z_peak)
	else:
		upload_info = "Please check, no file is uploaded!"
		return render_template('home.html', upload_info = upload_info)

# Plot figure for specific scan
@app.route('/plot.png', methods=['GET', 'POST'])
def plot_png():
	"""Generate the de-isotoped plot of MS2 scan."""
	from project_spectra.Spec_peptide import Peptide_identification
	import matplotlib.pyplot as plt
	import io
	import base64

	ms2_file = os.path.join(app.config['UPLOAD_FOLDER'], "ms2.mzML")
	scan_number = request.form['textbox']
	p = Peptide_identification(ms2_file, int(scan_number))
	p.store_one_scan_and_get_deisotoped()
	exp = MSExperiment()
	MzMLFile().load(p.deisotoped_file, exp)
	for spec in exp:
		for mz, i in zip(*spec.get_peaks()):
			plt.plot([mz, mz], [0, i], color='black')
			plt.text(mz, i, str(mz))

		# for the title add RT and Precursor m/z and charge info if available
		title = ''
		if spec.getRT() >= 0:
			title += 'RT: ' + str(spec.getRT())
		if len(spec.getPrecursors()) >= 1:
			title += '   Precursor m/z: ' + str(spec.getPrecursors()[0].getMZ()) + '   Charge: ' + str(spec.getPrecursors()[0].getCharge())
		plt.title(title)
		plt.ylabel('intensity')
		plt.xlabel('m/z')
		plt.ylim(bottom=0)

	img = io.BytesIO()
	plt.savefig(img, format='png')
	plt.close()
	img.seek(0)
	plot_url = base64.b64encode(img.getvalue()).decode('utf8')
	return render_template('plot.html', plot_url=plot_url)

# Determine peptides possible
@app.route('/peptide', methods=['POST'])
def get_peptide_list():
	"""get peptide list"""
	ms2_file = os.path.join(app.config['UPLOAD_FOLDER'], "ms2.mzML")
	scan_number = request.form['textbox']
	p = Peptide_identification(ms2_file, int(scan_number))

	try:
		vals = p.compile(show_C_terminal=False)[0]
		peaks = p.compile(show_C_terminal=False)[1]
		sum_intensity = p.compile(show_C_terminal=False)[2]
		lst = zip(vals, sum_intensity)
		show_seq_peak = zip(vals,peaks)
		result_p = pd.DataFrame(lst, columns=["Sequence", "Sum of Relative Peak Intensity"], dtype=float)
		result_delete_duplicate = result_p.drop_duplicates()
		data_table_html = result_delete_duplicate.to_html(header="true", table_id="table", index=False, justify="justify", )
		session['peptides'] = list(set(vals))
		return render_template('home.html', data_table_html=data_table_html, show_seq_peak=show_seq_peak, error_info="")
	except:
		error_info = "N/C terminal are not found by the algorithm or no enough AA to match, try another scan number :) "
		return render_template('home.html', error_info=error_info)

# About
@app.route("/about")
def about():
	return render_template("about.html")

# Upload
@app.route("/upload", methods=['GET', 'POST'])
def upload():
	if request.method == 'POST':
		file = request.files.get('file')
		clearAllFlag = request.values.get("clearall")

		# if the clear all button is pressed
		if clearAllFlag:
			clear_all()
			flash('Uploaded files successfully removed', "info")
			return redirect(url_for('upload'))

		# if a file was uploaded
		elif file:
			# if user does not select file, browser also
			# submit an empty part without filename
			if file.filename == '':
				flash('No file selected', "warning")
				return redirect(url_for('home'))

			if file and allowed_file(file.filename):
				filename = secure_filename(file.filename)
				file.save(Path(app.config['UPLOAD_FOLDER'], filename))
				print(file)
				flash(f"File uploaded successfully: '{file.filename}'", "info")
				return redirect(url_for('home'))
			else:
				flash(f"File could not be saved: '{file.filename}'", "warning")
				return redirect(url_for('upload'))

		# if no file was selected
		else:
			flash("Error! No file selected.", "warning")
			return redirect(url_for('upload'))
	else:
		return render_template("upload.html")


def allowed_file(filename):
	return '.' in filename and \
		filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def clear_all():
	for root, dirs, files in os.walk(UPLOAD_FOLDER):
		for file in files:
			os.remove(os.path.join(root, file))
	return

# Compute Protein Matches
@app.route("/proteins")
def protein_matches():
	# Run protein identification
	if session['peptides']:
		i = Identifier(sequence=session['peptides'])
		result = i.stout

		name = []
		full_name = []
		organism = []
		for x, y in result.items():
			name.append(x)
			full_name.append(y[0])
			organism.append(y[1])
		data = {'Accession Number': name, 'Fullname': full_name, 'Organism': organism}
		df = pd.DataFrame(data)
		data_table_html = df.to_html(header="true", table_id="table", index=False, justify="justify")
		return render_template("proteins.html", results = data_table_html)
	else:
		error_info = "No protein identified."
		return render_template("proteins.html", info = error_info)


'''
    Run app
'''
if __name__ == "__main__":
	flask_port = int(os.environ.get('FLASK_PORT', '5005'))
	app.run(host='0.0.0.0', port=flask_port, debug=True)
    #app.run(host='127.0.0.1', port=5005, debug=True)