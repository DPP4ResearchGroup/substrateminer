##################################################
# This module applies filter functions including
# - Mit or Cyto into designed categories
# - Extracellular Proteins
# - Non-Target/Target consensus Proteins
##################################################

import sys
import os

from substrate_miner.records.record_class import swissRecord, swissRecord2
from substrate_miner.records.record2seq import record2seq

from Bio import SeqIO


global GLOBAL_RECORD_TYPES
GLOBAL_RECORD_TYPES = ["text", "txt", "genbank", "swiss", "fasta"]
global GLOBAL_TARGET_PEPTIDES_TYPE
GLOBAL_TARGET_PEPTIDES_TYPE = ['SIGNAL', 'CHAIN', 'PEPTIDE']

# Define return class
class filteredRecords:
	def __init__(self) -> None:
		self.extracted_target_records = []
		self.extracted_non_target_records = []
		pass
	
	def _writeRecords(self, records: list, nontarget_records: list) -> dict[str, list]:
		return {"target": records, "nontarget": nontarget_records}

## Aux functions
def seq_error(seq_position) -> bool:
	'''
	Args:
		seq_position (int): The sequence position to check for errors.

	Returns:
		bool: True if an error is detected, False otherwise.
	Raises:
		ValueError: If an error is detected in the sequence position.
	'''

	if "UnknownPosition()" in str(seq_position):
		return True
	return False

# SeparationFilter Separation Functions in this Module
def seperate_target_size(records: list, target_size_min, target_size_max,\
	output_base: str, output_format: str, output_file_name, output_file_format,\
    stats) -> dict[str, list]:
	'''
	Args:
		records (list): A list of protein records in the form of UniProt records.
		target_size_min (int): The minimum target size to extract proteins from, 0 if not specified.
		target_size_max (int): The maximum target size to extract proteins from, 0 if not specified.
		output_format (str): The output format to save the extracted proteins.
		output_file_name (str): The output file name to save the extracted proteins.
		output_file_format (str): The output file format to save the extracted proteins.

	Returns:
		dict[str, list]: A dictionary containing the extracted target and non-target records.    
	Raises:
		None
	'''
 
	## params sanity checks
	# allowed ouput_formats: inline, file, all
	if output_format not in ["inline", "file", "all"]:
		raise ValueError("Unknown output format requested, allowed formats: inline, file, all.")

	if output_format == "file" and not output_file_name:
		raise ValueError("Output file name is required for file output format.")

	if output_format == "file" and (not output_file_format or output_file_format not in GLOBAL_RECORD_TYPES):
		raise ValueError(f"Unknown output file format requested, allowed formats: {GLOBAL_RECORD_TYPES}.")
 
	## params init
	# counter init
	i = 0
	j = 0
	
	# return variables init
	extracted_target_records = []
	extracted_non_target_records = []

	# output variables init
	unknown_records = []
 
	out_file = False
	if output_format != "inline":
		out_file = True
 
	# Extract records based on target size
	for record in records:
		# TODO: record conform
  		# record = swissRecord(record.accessions, record.entry_name, record.keywords, record.description, record.features, record.sequence)
     
		features_list = record.features
		for feature in features_list:
      
			## search sequence features in CHAIN, PEPTIDE, SIGNAL
			if feature.type in GLOBAL_TARGET_PEPTIDES_TYPE:
				if seq_error(feature.location.start) or seq_error(feature.location.end):
					print (">>>>>>>>>>>>>>>>>>>>>>>>>")
					print("Unknown Position detected ...")
					print(">>>>>>>>>>>>>>>>>>>>>>>>>")
     
					unknown_records.append(record)

					continue
				
    			# debug block
				#print("%%%%%%%%%%%%%%%%%%%%%%%%")
				#print(feature.location.start, feature.location.end)
				#print(target_size_max, target_size_min, abs(feature.location.end - feature.location.start))
				#print("%%%%%%%%%%%%%%%%%%%%%%%%")
    
				if (target_size_max != 0 and target_size_min == 0 and abs(feature.location.end - feature.location.start) <= target_size_max) or\
					(target_size_max == 0 and target_size_min != 0 and abs(feature.location.end - feature.location.start) >= target_size_min) or\
					(target_size_min != 0 and target_size_max != 0 and abs(feature.location.end - feature.location.start) >= target_size_min and\
         				abs(feature.location.end - feature.location.start) <= target_size_max):
         
					extracted_target_records.append(record)
					i += 1

				else:
					extracted_non_target_records.append(record)
					j += 1
    
    ## output block
	target_unique = f"{target_size_min}_{target_size_max}"
	if out_file:
		write_to_file_batch(target_unique,
							extracted_target_records, extracted_non_target_records, unknown_records,
							output_base, output_file_name, output_file_format)
		
	## status block
	if stats:
		print("++++++++++++++++++++++++++++++++++++++++++++")
		print(f"Total Target Records: {i}")
		print(f"Total Non-Target Records: {j}")
		print("++++++++++++++++++++++++++++++++++++++++++++")

	if out_file and len(unknown_records) > 0:
		print("++++++++++++++++++++++++++++++++++++++++++++")
		print(f"Total Unknown Records: {len(unknown_records)} ...")
		print("++++++++++++++++++++++++++++++++++++++++++++")

	## conform return variables before exporting
	filtered_records = filteredRecords()
	for record in extracted_target_records:
		print(f"Record: {record.entry_name} Registered ...")
  
	return filtered_records._writeRecords(records=extracted_target_records, nontarget_records=extracted_non_target_records)



def seperate_cellular_location(records: list, target_orgin: str,\
  output_base: str, output_format: str, output_file_name, output_file_format,\
  stats) -> dict[str, list]:
	'''
	Args:
		target_orgin (str): The target cellular location to extract proteins from.
		output_format (str): The format of the output, either "inline" or "file".
		output_file_name (str): The name of the output file if output_format is "file".
		extracted_target_records (list): A list to store protein records that match the target origin.
		extracted_non_target_records (list): A list to store protein records that do not match the target origin.
	
	Returns:
		dict[str, list]: A dictionary containing the extracted target and non-target records.
	Raises:
		None
	'''

	## params sanity checks
	out_file = False
	if output_format != "inline":
		out_file = True
	
	# allowed ouput_formats: inline, file, all
	if output_format not in ["inline", "file", "all"]:
		raise ValueError("Unknown output format requested, allowed formats: inline, file, all.")

	if output_format == "file" and not output_file_name:
		raise ValueError("Output file name is required for file output format.")

	if output_format == "file" and (not output_file_format or output_file_format not in GLOBAL_RECORD_TYPES):
		raise ValueError(f"Unknown output file format requested, allowed formats: {GLOBAL_RECORD_TYPES}.")
	
	# return variables init
	extracted_target_records = []
	extracted_non_target_records = []
 
	# output variables init
	unknown_records = []
 
	# counter init
	i = 0
	j = 0
	
	# allowed ouput_formats: inline, file, all
	if output_format != "inline" and not output_file_name:
		raise ValueError("Output file name is required for file output format ...")

	if output_format != "inline":
		out_file = True
  
	# Extract records based on target origin
	for record in records:
		if target_orgin in record.keywords or target_orgin.capitalize() in record.keywords or target_orgin.upper() in record.keywords:
			extracted_target_records.append(record)
			i += 1

		else:
			extracted_non_target_records.append(record) 
			j += 1
    
    ## output block
	if out_file:
		write_to_file_batch(target_orgin,\
      		extracted_target_records, extracted_non_target_records, unknown_records,\
            output_base, output_file_name, output_file_format)
		
	## status block
	if stats:
		print("++++++++++++++++++++++++++++++++++++++++++++")
		print(f"Total Target Records: {i}")
		print(f"Total Non-Target Records: {j}")
		print("++++++++++++++++++++++++++++++++++++++++++++")

	if out_file and len(unknown_records) > 0:
		print("++++++++++++++++++++++++++++++++++++++++++++")
		print(f"Total Unknown Records: {len(unknown_records)} ...")
		print("++++++++++++++++++++++++++++++++++++++++++++")
 
    ## conform return variables before exporting
	filtered_records = filteredRecords()
	for record in extracted_target_records:
		print(f"Record: {record.entry_name} Registered ...")

	return filtered_records._writeRecords(records=extracted_target_records, nontarget_records=extracted_non_target_records)


def seperate_target_motif(records: list,\
  target_motif_features: dict,\
  output_base: str, output_format: str, output_file_name, output_file_format,\
  stats) -> dict[str, list]:
	'''
    Separates records based on a target motif and outputs the results in the specified format.
    
	Args:
		records (list): List of swissRecord objects to be processed.
		target_motif (str): The motif sequence to search for within the records.
		target_start_position (int): The starting position of the target motif.
		target_motif_type (str): The type of motif to search for (e.g., 'SIGNAL', 'CHAIN', 'PEPTIDE').
		output_format (str): The format for output ('inline', 'file', 'all').
		output_file_name (str): The name of the output file if output_format is 'file' or 'all'.
		output_file_format (str): The format of the output file.
		stats: Additional statistics or metadata to be included in the output.
	
	Returns:
		dict[str, list]: A dictionary with keys 'extracted_target_records' and 'extracted_non_target_records',
		each containing a list of swissRecord objects that match or do not match the target motif, respectively.

	Raises:
		ValueError: If an unknown output format is requested.
		ValueError: If the output file name is required for file output format.
	'''

	## params sanity checks
	out_file = False
	if output_format != "inline":
		out_file = True

	# allowed ouput_formats: inline, file, all
	if output_format not in ["inline", "file", "all"]:
		raise ValueError("Unknown output format requested, allowed formats: inline, file, all.")

	if output_format == "file" and not output_file_name:
		raise ValueError("Output file name is required for file output format.")
	
	if output_format == "file" and (not output_file_format or output_file_format not in GLOBAL_RECORD_TYPES):
		raise ValueError(f"Unknown output file format requested, allowed formats: {GLOBAL_RECORD_TYPES}.")
	
	# counter init
	i = 0
	j = 0

	# return variables init
	extracted_target_records = []
	extracted_non_target_records = []
	
	# output variables init
	unknown_records = []
 	
	# debug block
	#print("%%%%%%%%%%%%%%%%%%%%%%%%")
	#print(target_motif_features)
	#print(target_motif_features.motif)
	#print("%%%%%%%%%%%%%%%%%%%%%%%%")
 	
	# Extract target motif
	target_motif = target_motif_features.motif
	motif_len = len(target_motif)  # Removed unused variable

	# Ensure the motif contains only valid amino acids
	valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWYX")
	if not all(aa in valid_amino_acids for aa in target_motif):
		raise ValueError("The target motif contains invalid amino acids. For unknown amino acids, use 'X'.")
 
	for record in records:
		# record conform
		# record = swissRecord2(record.accessions, record.entry_name, record.keywords, record.description, record.features, record.sequence, record.organism)
  
		features_list = record.features
		for feature in features_list:
			if feature.type in GLOBAL_TARGET_PEPTIDES_TYPE:
				# motif dependent search block
				if target_motif_features.type == "exo":
					if target_motif_features.terminal == "N":
						feature_start_position = feature.location.start
						feature_end_position = motif_len

					elif target_motif_features.terminal == "C":
						feature_start_position = feature.location.end - motif_len
						feature_end_position = feature.location.end
					else:
						raise ValueError("Unknown terminal requested, for exo-terminal search, use 'N' or 'C', checkin config file ...")

					# debug block
					#print("%%%%%%%%%%%%%%%%%%%%%%%%")
					#print(seq_error(feature_start_position), feature_end_position)
					#print(record.sequence[feature_start_position:feature_end_position])
					#print(target_motif)
					#print("%%%%%%%%%%%%%%%%%%%%%%%%")
	
					try:
						if seq_error(feature_start_position) or seq_error(feature_end_position):
							unknown_records.append(record)

							print ("++++++++++++++++++++++++")
							print (f"Found Incomplete Record: {record.entry_name}")
							print ("++++++++++++++++++++++++")
							continue
					except ValueError:
  							raise ValueError("Unable to handle Position detected ...")

					# debug block
					#print("%%%%%%%%%%%%%%%%%%%%%%%%")
					#print(feature_start_position, feature_end_position)
					#print(record.sequence[feature_start_position:feature_end_position])
					#print(target_motif)
					#print("%%%%%%%%%%%%%%%%%%%%%%%%")

					if record.sequence[feature_start_position:feature_end_position] == target_motif:
						extracted_target_records.append(record)
						i += 1

					else:
						extracted_non_target_records.append(record)
						j += 1
       
					# debug block
					#print("%%%%%%%%%%%%%%%%%%%%%%%%")
					#print(f"Target: {extracted_target_records}")
					#print(f"Non-Target: {extracted_non_target_records}")
					#print("%%%%%%%%%%%%%%%%%%%%%%%%")

				elif target_motif_features.type == "endo":
					motif_status = True
					for key in target_motif_features.residue:
						if motif_status and record.sequence[key] == target_motif_features.residue[key]:
							pass
						else:
							extracted_non_target_records.append(record)
							j += 1

					if motif_status:
						extracted_target_records.append(record)
						i += 1

				else:
					raise ValueError("Unknown motif type requested, for exo-terminal search, use 'exo' or 'endo', checkin config file ...")

	## output block
	if out_file:
		write_to_file_batch(target_motif,\
      		extracted_target_records, extracted_non_target_records, unknown_records,\
            output_base, output_file_name, output_file_format)
		
	## status block
	if stats:
		print("++++++++++++++++++++++++++++++++++++++++++++")
		print(f"Total Target Records: {i}")
		print(f"Total Non-Target Records: {j}")
		print("++++++++++++++++++++++++++++++++++++++++++++")

	if out_file and len(unknown_records) > 0:
		print("++++++++++++++++++++++++++++++++++++++++++++")
		print(f"Total Unknown Records: {len(unknown_records)} ...")
		print("++++++++++++++++++++++++++++++++++++++++++++")

	## conform return variables before exporting
	filtered_records = filteredRecords()
	for record in extracted_target_records:
		print(f"Record: {record.entry_name} Registered ...")

	return filtered_records._writeRecords(records=extracted_target_records, nontarget_records=extracted_non_target_records)


def write_to_file(records, output_file_name: str, output_file_format: str) -> None:
	'''
	Args:
		records (list): A list of protein records to be written to a file.
		output_file_name (str): The name of the output file to write the records to.
		output_file_format (str): The format of the output file to write the records to.
	
	Returns:
		None
	Raises:
		None
	'''

	if output_file_format == "text" or output_file_format == "txt":
		with open(output_file_name + ".txt", "w") as f:
			for record in records:
				f.write(f"{record.accessions}\n")
				f.write(f"{record.entry_name}\n")
				f.write(f"{record.keywords}\n")
				f.write(f"{record.description}\n")
				f.write(f"{record.sequence}\n")
				f.write("\n")

    ## TODO: to be tested "genbank" format
	elif output_file_format == "genbank":
		with open(output_file_name + ".gb", "w") as f:
			SeqIO.write(records, f, "genbank")
	elif output_file_format == "fasta":
		seqobjs = []
		with open(output_file_name + ".fasta", "w") as f:
			for record in records:
				seqobj = record2seq(record)
				seqobjs.append(seqobj)
			SeqIO.write(seqobjs, f, "fasta")

	else:
		raise ValueError("Unknown format requested ...")


def write_to_file_batch(target: str,\
    extracted_target_records, extracted_non_target_records, unknown_records,\
    output_base: str, output_file_name: str, output_file_format: str) -> None:
	'''
	Args:
		target_motif (str): The target motif to search for within the records.
		extracted_target_records (list): A list of protein records that match the target motif.
		extracted_non_target_records (list): A list of protein records that do not match the target motif.
		unknown_records (list): A list of protein records with unknown positions.
		output_base (str): The base path to save the output files.
		output_file_name (str): The name of the output file to write the records to.
		output_file_format (str): The format of the output file to write the records to.
	
	Returns:
		None
	Raises:
		None
	'''

	output_path_target = os.path.join(output_base, output_file_name + f"_target_{target}" + ".out")
	output_path_nontarget = os.path.join(output_base, output_file_name + f"_nontarget_{target}" + ".out")
	output_path_unknown = os.path.join(output_base, output_file_name + f"_unknown_{target}" + ".out")
 
	write_to_file(extracted_target_records,\
     	output_path_target, output_file_format)
	write_to_file(extracted_non_target_records,\
     	output_path_nontarget, output_file_format)
	if len(unknown_records) > 0:
		write_to_file(unknown_records,\
      		output_path_unknown, output_file_format)