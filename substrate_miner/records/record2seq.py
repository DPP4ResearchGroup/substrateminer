from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def record2seq(swissprot_record) -> SeqRecord:
    """
    Convert a SwissProt record to a SeqRecord object.
    Args:
        swissprot_record: A SwissProt record object containing protein sequence data and annotations.
    Returns:
        SeqRecord: A SeqRecord object containing the sequence, ID, name, description, annotations, and features from the SwissProt record.
    Annotations include:
        - organism: The organism from which the protein is derived.
        - taxonomy: The taxonomic classification of the organism.
        - created: The date the record was created.
        - sequence_length: The length of the protein sequence.
        - data_class: The data class of the record.
        - keywords: Keywords associated with the record.
    Features include:
        - location: The start and end positions of the feature.
        - type: The type of the feature.
        - description: The description of the feature.
    """
    
    # Extract necessary fields from the SwissProt record
    seq = Seq(swissprot_record.sequence)
    id = swissprot_record.accessions[0] if swissprot_record.accessions else "UNKNOWN"
    name = swissprot_record.entry_name
    description = swissprot_record.description
    annotations = {
        "organism": swissprot_record.organism,
        "taxonomy": swissprot_record.organism_classification,
        # "sequence_version": swissprot_record.sequence_version,
        # "entry_version": swissprot_record.entry_version,
        "created": swissprot_record.created,
        "sequence_length": swissprot_record.seqinfo[0],
        "data_class": swissprot_record.data_class,
        "keywords": swissprot_record.keywords
    }

    # Create a SeqRecord object
    seq_record = SeqRecord(
        seq=seq,
        id=id,
        name=name,
        description=swissprot_record.description,
        annotations=annotations
    )

    for feature in swissprot_record.features:
        location = FeatureLocation(start=feature.location.start, end=feature.location.end)
        qualifiers = {
            "type": feature.type,
            "description": swissprot_record.description
        }
        seq_feature = SeqFeature(location=location, type=feature.type, qualifiers=qualifiers)
        seq_record.features.append(seq_feature)
    
    return seq_record
