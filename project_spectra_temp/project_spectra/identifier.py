# -*- coding: utf-8 -*-

"""This module is used to request protein database for protein identification.
   We'll get a dictionary which contains accession number as key, protein full name and organism in a list as value.
"""

from typing import List
import json
import logging
import requests
from sqlalchemy.orm import sessionmaker
from project_spectra.constants import PEPTIDE_MATCH_API
from project_spectra.models import Base, engine, Protein, Peptide


class Identifier:
	"""Retrieves information from a given protein sequence"""
	def __init__(self, sequence:List):
		self.sequence = sequence
		self.results = []
		self.stout = []
		self.all_results = []

		# Run methods automatically
		self.prot_identifier()

	def protein_match_api(self, url_API: str, query: str):
		"""Get protein information from bioinformatics.udel.edu via API"""
		requestURL = url_API + str(query)
		response_API = requests.get(requestURL, headers={"Accept": "application/json"})

		# Check the response status code
		if response_API.status_code == 200:
			logging.info('Successfully received data. API: {}'.format(requestURL))
		else:
			logging.error('Unsuccessful request. API: {} Status Code:{}'.format(requestURL, response_API.status_code))

		# Get the data from the API
		data = response_API.content
		text = response_API.text

		# Parse it
		json_data = json.loads(data)

		return json_data, text

	def extract(self, json_data):
		"""Extracts information from peptide match json to list of dictionaries"""
		if json_data.get('results'):
			for each in json_data['results']:
				if each.get('proteins'):
					for prot in each.get('proteins'):
						# Add everything to dic
						new_prot = {
								"ac": prot.get("ac"),
								"id": prot.get("id"),
								"name": prot.get("name"),
								"organism": prot.get("orgName"),
								"taxon": prot.get("orgTaxonId"),
								"sequence": prot.get("sequence")
							}
						self.results.append(new_prot)

	def prot_identifier(self):
		"""Identifies a protein given a peptide sequence"""
		# Start database session
		Base.metadata.create_all(engine, checkfirst=True)
		Session = sessionmaker(bind=engine)
		session = Session()

		# Check if peptide in already in db
		query = session.query(Peptide).filter(Peptide.seq == ','.join(self.sequence)).first()

		if query:
			self.all_results = [{'accession_number': each.accession_number,
								'name' : each.name,
								'full_name' : each.full_name,
								'organism': each.organism,
								'tax_id': each.tax_id ,
								'sequence' : each.sequence
								 } for each in query.proteins ]

			#self.stout = [each.get('name') for each in self.all_results]
			self.stout = {each.get('name'): [each.get('full_name'), each.get('organism')] for each in self.all_results}
			return

		# Run peptide match API
		all_peptides = '%2C'.join(self.sequence) # Join with a comma value

		query = all_peptides + "&swissprot=false&isoform=false&uniref100=false&leqi=false&offset=0"
		json_data, text = self.protein_match_api(url_API=PEPTIDE_MATCH_API, query=query)

		# Check the response is not empty
		if json_data.get('numberFound') == 0:
			logging.warning('No results found for this peptide: {}.'.format(self.sequence))
			return

		# Extract data into 'self.results'
		self.extract(json_data=json_data)

		# Add peptide to the database
		if self.results:
			peptide = Peptide(seq=','.join(self.sequence))
			session.add(peptide)
			for result in self.results:
				prot = Protein(accession_number=result.get("ac"),
							name=result.get("id"),
							full_name=result.get("name"),
							organism=result.get("organism"),
							tax_id=result.get("taxon"),
							sequence=result.get("sequence"))
				session.add(prot)
				peptide.proteins.append(prot)

				self.all_results.append({'accession_number': result.get("ac"),
										 'name': result.get("id"),
										 'full_name': result.get("name"),
										 'organism': result.get("organism"),
										 'tax_id': result.get("taxon"),
										 'sequence': result.get("sequence")
										 })
			session.commit()
			logging.info('Successfully added peptides to the database. Peptide: {}'.format(self.sequence))
			#self.stout = [each.get('name') for each in self.all_results]
		self.stout = {each.get('name'):[each.get('full_name'), each.get('organism')] for each in self.all_results}
		return

