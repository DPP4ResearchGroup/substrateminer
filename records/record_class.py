#####################################################################################################################
# This module defines the Record, Feature, and Location classes for handling records.
# The Record class includes methods for writing records in various formats such as fasta, text, genbank, and swiss.
#####################################################################################################################

# Importing required modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

class filterOutput:
	def __init__(self, target_records: list, non_target_records: list) -> dict:
		self.target = target_records
		self.non_target = non_target_records
		# self.unknown = unknown_records
		pass	
 
class emblRecord:
	def __init__(self, seq: Seq, id: str, name: str, description, annotations, features) -> None:
		self.seq = seq
		self.id = id
		self.name = name
		self.description = description
		self.annotations = annotations
		self.features = features
		pass

# Definition of uniprot record class
class swissRecord:
	def __init__(self, accessions, entry_name, keywords, description, features, sequence) -> None:
		self.accessions = accessions
		self.entry_name = entry_name
		self.keywords = keywords
		self.description = description
		self.features = features
		self.sequence = sequence
		pass
	
class swissRecord2:
	def __init__(self, accessions, entry_name, keywords, description, features, sequence, organism) -> None:
		self.accessions = accessions
		self.entry_name = entry_name
		self.keywords = keywords
		self.description = description
		self.features = features
		self.sequence = sequence
		self.organism = organism
		pass
	
	def _sR2eR(self) -> list[emblRecord]:
		# Extract data from the Swiss-Prot record
		seq = Seq(self.sequence)
		description = self.description
		accession = self.accessions[0] if self.accessions else "UNKNOWN"
		name = self.entry_name
		organism = self.organism
		features = []

		# Map Swiss-Prot features to SeqFeature objects
		for feature in self.features:
			location = FeatureLocation(feature.location.start, feature.location.end)
			qualifiers = {
				"type": feature.type,
				"description": feature.type,
			}
			seq_feature = SeqFeature(location=location, type=feature.type, qualifiers=qualifiers)
			features.append(seq_feature)

		# Create a SeqRecord for EMBL conversion
		seq_records: emblRecord = [emblRecord(
		    seq=seq,
		    id=accession,
		    name=name,
		    description=description,
		    annotations={
          		"organism": organism,
            	"molecule_type": "protein",
            	"topology": "linear",
            },
		    features=features,
		)]

		return seq_records

	def _writeRecord(self, out_file, format) -> None:
		if format == "fasta":
			fasta_record = SeqRecord(Seq(self.sequence), id=self.accessions, name=self.entry_name, description=self.description)
			SeqIO.write(fasta_record, out_file, "fasta")

		elif format == "text" or format == "txt":
			with open(out_file, "w") as f:
				f.write(">{0}\n".format(self.accessions))
				f.write("\t{0}\n".format(self.entry_name))
				f.write("\t{0}\n".format(self.keywords))
				f.write("\t{0}\n".format(self.description))
				f.write("\t{0}\n".format(self.sequence))
				f.write("\n")
       
			print ("++++++++++++++++++++")
			print (f"Record: {self.entry_name} Registered ...")

		## TODO: genbank move to genomic class
		## TODO: this class not working under swiss conditions
		elif format == "genbank":
			gb_record = SeqRecord(Seq(self.sequence), id=self.accessions, name=self.entry_name, description=self.description)
			SeqIO.write(gb_record, out_file, "genbank")
   
			print ("++++++++++++++++++++")
			print (f"Record: {gb_record.entry_name} Registered ...")
   
		elif format == "swiss":
			embl_record_handle = self._sR2eR()
			SeqIO.write(embl_record_handle, out_file, "embl")
   
			print ("++++++++++++++++++++")
			print (f"Record: {record.entry_name} Registered ...")
   
		else:
			print (">>>>>>>>>>>>>>>>>>>>>>>>>")
			raise ValueError("Unknown format requested ...")
		

class genBank:
	def __init__(self, accessions, entry_name, keywords, description, features, sequence) -> None:
		self.accessions = accessions
		self.entry_name = entry_name
		self.keywords = keywords
		self.description = description
		self.features = features
		self.sequence = sequence
		pass

	def _writeRecord(self, record, out_file, format) -> None:
		if format == "fasta":
			fasta_record = SeqRecord(Seq(record.sequence), id=record.accessions, name=record.entry_name, description=record.description)
			SeqIO.write(fasta_record, out_file, "fasta")

		elif format == "text" or format == "txt":
			with open(out_file, "w") as f:
				f.write(">{0}\n".format(record.accessions))
				f.write("\t{0}\n".format(record.entry_name))
				f.write("\t{0}\n".format(record.keywords))
				f.write("\t{0}\n".format(record.description))
				f.write("\t{0}\n".format(record.sequence))
				f.write("\n")
	   
			print ("++++++++++++++++++++")
			print (f"Record: {record.entry_name} Registered ...")
	
		elif format == "genbank":
			SeqIO.write(record, out_file, "genbank")
   
			print ("++++++++++++++++++++")
			print (f"Record: {record.entry_name} Registered ...")
   
		elif format == "swiss":
			swiss_record = SeqRecord(Seq(record.sequence), id=record.accessions, name=record.entry_name, description=record.description)
			SeqIO.write(swiss_record, out_file, "swiss")
   
			print ("++++++++++++++++++++")
			print (f"Record: {swiss_record.entry_name} Registered ...")
   
		else:
			print (">>>>>>>>>>>>>>>>>>>>>>>>>")
			raise ValueError("Unknown format requested ...")

# Definition of feature class
class seqFeature:
	def __init__(self, type, location) -> None:
		self.type = type
		self.location = location
		pass

# Definition of location class
class seqLocation:
	def __init__(self, start, end) -> None:
		self.start = start
		self.end = end
		pass

# Definition of motif features class
class motifFeatures:
	def __init__(self, type, motif, terminal) -> None:
		self.type = type
		self.motif = motif
		self.terminal = terminal

		features = self._feature(type, motif)
		self.features = features
		pass

	def _feature(self, type, motif) -> dict:
		if type == "exo":
			motif_start = 0
			motif_end = len(motif)
			exo_feature = {"type": type, "start": motif_start, "end": motif_end}
			return exo_feature	
   
		elif type == "endo":
			motif_dict = {c: i for i, c in enumerate(motif) if c != 'X'}
			endo_feature = {"type": type, "residues": motif_dict}
			return endo_feature
		else:
			raise ValueError("Unknown motif type requested ...")