# -*- coding: utf-8 -*-
from pathlib import Path

from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String, ForeignKey, Table
from sqlalchemy.orm import declarative_base, relationship

# Init project directory
HIDDEN_FOLDER = Path(Path.home(), '.projectSpectra')
HIDDEN_FOLDER.mkdir(parents=True, exist_ok=True)
DATABASE_PATH = Path(HIDDEN_FOLDER, 'spectra_db.db')

# Create engine
engine = create_engine('sqlite:///'+str(DATABASE_PATH), echo = False)
Base = declarative_base()

# Association tables
relationship_table = Table('tags_association', Base.metadata,
    Column('pept_id', ForeignKey('peptide.id'), primary_key=True),
    Column('prot_id', ForeignKey('protein.id'), primary_key=True),
)


class Peptide(Base):
	__tablename__ = 'peptide'
	id = Column(Integer, primary_key=True)
	seq = Column(String)
	proteins = relationship("Protein", secondary=relationship_table, back_populates="peptides")  # many to many rel

	def __init__(self, seq):
		self.seq = seq

class Protein(Base):
	__tablename__ = 'protein'
	id = Column(Integer, primary_key=True)
	accession_number = Column(String(255))
	name = Column(String(255))
	full_name = Column(String(255))
	organism = Column(String(225))
	tax_id = Column(Integer)
	sequence = Column(String)
	peptides = relationship("Peptide", secondary=relationship_table, back_populates="proteins")  # many to many rel

	def __init__(self, accession_number, name, full_name, organism, tax_id, sequence):
		self.accession_number = accession_number
		self.name = name
		self.full_name = full_name
		self.organism = organism
		self.tax_id = tax_id
		self.sequence = sequence