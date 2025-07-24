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
	"""
	A class to store filtered records from a database query.
	Attributes:
		target (list): A list of records that match the target criteria.
		non_target (list): A list of records that do not match the target criteria.
		# unknown (list): A list of records that are not classified as target or non-target.
	"""
 
	def __init__(self, target_records: list, non_target_records: list) -> dict:
		self.target = target_records
		self.non_target = non_target_records
		# self.unknown = unknown_records
		pass	
 
class emblRecord:
	"""
	A class to represent an EMBL record.
	Attributes:
		seq (Seq): The sequence of the record.
		id (str): The identifier of the record.
		name (str): The name of the record.
		description (str): The description of the record.
		annotations (dict): A dictionary of annotations for the record.
		features (list): A list of SeqFeature objects representing features of the record.
	"""
 
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
	"""
	A class to represent a Swiss-Prot record.
	Attributes:
		accessions (list): A list of accession numbers for the record.
		entry_name (str): The entry name of the record.
		keywords (list): A list of keywords associated with the record.
		description (str): The description of the record.
		features (list): A list of SeqFeature objects representing features of the record.
		sequence (str): The sequence of the record.
	"""
  
	def __init__(self, accessions, entry_name, keywords, description, features, sequence) -> None:
		self.accessions = accessions
		self.entry_name = entry_name
		self.keywords = keywords
		self.description = description
		self.features = features
		self.sequence = sequence
		pass
	
class swissRecord2:
	"""
	A class to represent a Swiss-Prot record with additional organism information.
	Attributes:
		accessions (list): A list of accession numbers for the record.
		entry_name (str): The entry name of the record.
		keywords (list): A list of keywords associated with the record.
		description (str): The description of the record.
		features (list): A list of SeqFeature objects representing features of the record.
		sequence (str): The sequence of the record.
		organism (str): The organism associated with the record.
	"""
 
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
		"""
		Converts a Swiss-Prot record to an EMBL record format.
		Returns:
			list[emblRecord]: A list containing an emblRecord object representing the Swiss-Prot record.
		"""
  
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
		"""
		Writes the Swiss-Prot record to a file in the specified format.
		Args:
			out_file (str): The output file path where the record will be written.
			format (str): The format in which to write the record. Options are "fasta", "text", "genbank", or "swiss".
		Raises:
			ValueError: If an unknown format is requested.
		"""
  
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
	"""
	A class to represent a GenBank record.
	Attributes:
		accessions (list): A list of accession numbers for the record.
		entry_name (str): The entry name of the record.
		keywords (list): A list of keywords associated with the record.
		description (str): The description of the record.
		features (list): A list of SeqFeature objects representing features of the record.
		sequence (str): The sequence of the record.
	"""
 
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
	"""
	A class to represent motif features in a sequence.
	Attributes:
		type (str): The type of motif, either "exo" or "endo".
		motif (str): The motif sequence.
		terminal (bool): Indicates if the motif is terminal.
		features (dict): A dictionary containing features of the motif.
	"""
 
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