from src.base_modifications import BaseModification
from src.single_stranded_nucleic_acids import SingleStrandNucleicAcidSequence
from src.sequence_annotations import SequenceAnnotation


def test_make_ss_dna():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
    assert dna_seq.sequence == "ATGCGTAA"
    assert dna_seq.nucleic_acid_type == "DNA"
    assert dna_seq.circular == False
    assert dna_seq.strand_direction == "forward"


def test_make_ss_circular_rna():
    rna_seq = SingleStrandNucleicAcidSequence(sequence="AUGCGUAA", nucleic_acid_type="RNA", circular=True, strand_direction="forward")
    assert rna_seq.sequence == "AUGCGUAA"
    assert rna_seq.nucleic_acid_type == "RNA"
    assert rna_seq.circular == True
    assert rna_seq.strand_direction == "forward"


def test_make_ss_reverse_dna_with_ambiguous_base():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATNNNTAA", nucleic_acid_type="DNA", circular=False, strand_direction="reverse")
    assert dna_seq.sequence == "ATNNNTAA"
    assert dna_seq.nucleic_acid_type == "DNA"
    assert dna_seq.circular == False
    assert dna_seq.strand_direction == "reverse"


def test_make_ss_dna_with_annotations():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGNNNNNTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward",
                                              annotations=[{"name": "Gene1", "start": 0, "end": 3, "note": "This is Gene1"},
                                                            {"name": "random_barcode", "start": 5, "end": 9, "note": "This is a random barcode"}])
    assert dna_seq.annotations["Gene1"].name == "Gene1"
    assert dna_seq.annotations["Gene1"].start == 0
    assert dna_seq.annotations["Gene1"].end == 3
    assert dna_seq.annotations["Gene1"].note == "This is Gene1"
    assert dna_seq.annotations["random_barcode"].name == "random_barcode"
    assert dna_seq.annotations["random_barcode"].start == 5
    assert dna_seq.annotations["random_barcode"].end == 9
    assert dna_seq.annotations["random_barcode"].note == "This is a random barcode"


def test_make_ss_dna_with_base_modifications():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGNNNNNTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward",
                                              base_modifications=[{"position": 3, "modification_type": "5-mC"},
                                                                  {"position": 0, "modification_type": "6-mA"}])
    assert dna_seq.base_modifications[3].position == 3
    assert dna_seq.base_modifications[3].modification_type == "5-mC"
    assert dna_seq.base_modifications[0].position == 0
    assert dna_seq.base_modifications[0].modification_type == "6-mA"


def test_make_ss_dna_with_cutsites_and_base_mods():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGGAATTCTAGCATGCAAATT", nucleic_acid_type="DNA", circular=False, strand_direction="forward",
                                              base_modifications=[{"position": 10, "modification_type": "5-mC"},
                                                                  {"position": 0, "modification_type": "6-mA"}])
    assert dna_seq.cut_sites[13].start == 13
    assert dna_seq.cut_sites[13].end == 18
    assert dna_seq.cut_sites[13].cut_position == 18
    assert dna_seq.cut_sites[13].recognition_sequence == "GCATGC"
    assert dna_seq.cut_sites[13].restriction_enzyme == "SphI"
    assert dna_seq.base_modifications[10].position == 10
    assert dna_seq.base_modifications[10].modification_type == "5-mC"
    assert dna_seq.base_modifications[0].position == 0
    assert dna_seq.base_modifications[0].modification_type == "6-mA"


def test_make_ss_dna_with_cutsites_then_add_edit_and_remove_base_mods():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGGAATTCTAGCATGCAAATT", nucleic_acid_type="DNA", circular=False, strand_direction="forward")

    dna_seq.add_base_modifications(
        [{"position": 10, "modification_type": "5-mC"},
         {"position": 0, "modification_type": "6-mA"}])
    dna_seq.edit_base_modification(0, position=15)
    dna_seq.remove_base_modifications([10])

    assert dna_seq.cut_sites[5].start == 5
    assert dna_seq.cut_sites[5].end == 10
    assert dna_seq.cut_sites[5].cut_position == 6
    assert dna_seq.cut_sites[5].recognition_sequence == "GAATTC"
    assert dna_seq.cut_sites[5].restriction_enzyme == "EcoRI"



# # Sample test cases

# # Valid DNA sequence
# try:
#     dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
#     print("Valid DNA Sequence Created:", dna_seq)
# except Exception as e:
#     print("Error creating DNA sequence:", e)

# # Valid RNA sequence
# try:
#     rna_seq = SingleStrandNucleicAcidSequence(sequence="AUGCGUAA", nucleic_acid_type="RNA", circular=True, strand_direction="forward")
#     print("Valid RNA Sequence Created:", rna_seq)
# except Exception as e:
#     print("Error creating RNA sequence:", e)

# # Invalid sequence characters
# try:
#     invalid_seq = SingleStrandNucleicAcidSequence(sequence="ATBZ", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
#     print("Invalid Sequence Created:", invalid_seq)
# except Exception as e:
#     print("Error creating invalid sequence:", e)

# # DNA sequence with 'U'
# try:
#     dna_with_u = SingleStrandNucleicAcidSequence(sequence="ATGCUA", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
#     print("DNA Sequence with U Created:", dna_with_u)
# except Exception as e:
#     print("Error creating DNA sequence with U:", e)

# # RNA sequence with 'T'
# try:
#     rna_with_t = SingleStrandNucleicAcidSequence(sequence="AUGCTA", nucleic_acid_type="RNA", circular=False, strand_direction="forward")
#     print("RNA Sequence with T Created:", rna_with_t)
# except Exception as e:
#     print("Error creating RNA sequence with T:", e)

# try:
#     ssDNA = SingleStrandNucleicAcidSequence(sequence="ATGCGNNNNNTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
#     ssDNA.add_annotations([
#         {"name": "Gene1", "start": 0, "end": 3, "note": "This is Gene1"},
#         {"name": "random_barcode", "start": 5, "end": 9, "note": "This a random barcode"},
#     ])
#     print("Annotations added successfully:", ssDNA.annotations)
# except Exception as e:
#     print("Error adding annotations:", e)

# # Test editing an annotation
# try:
#     ssDNA.edit_annotation("Gene1", start=1, end=4)
#     print("Annotation edited successfully:", ssDNA.annotations['Gene1'])
# except Exception as e:
#     print("Error editing annotation:", e)

# # Test removing an annotation
# try:
#     ssDNA.remove_annotations(["Gene1"])
#     print("Annotation removed successfully, remaining annotations:", ssDNA.annotations)
# except Exception as e:
#     print("Error removing annotation:", e)

# # Test adding base modifications
# try:
#     ssDNA.add_base_modifications([
#         {"position": 3, "modification_type": "5-mC"}
#     ])
#     print("Base modifications added successfully:", ssDNA.base_modifications)
# except Exception as e:
#     print("Error adding base modifications:", e)

# # Test editing a base modification
# try:
#     ssDNA.edit_base_modification(3, modification_type="5-hmC")
#     print("Base modification edited successfully:", ssDNA.base_modifications[3])
# except Exception as e:
#     print("Error editing base modification:", e)

# # Test removing a base modification
# try:
#     ssDNA.remove_base_modifications([3])
#     print("Base modification removed successfully, remaining:", ssDNA.base_modifications)
# except Exception as e:
#     print("Error removing base modification:", e)

# # Test reverse complement
# try:
#     print("Original sequence:", ssDNA.sequence)
#     ssDNA.reverse_complement()
#     print("Reverse complement successful, new sequence:", ssDNA.sequence)
# except Exception as e:
#     print("Error computing reverse complement:", e)

# # Test setting circular
# try:
#     ssDNA.set_circular()
#     print("Set circular successful, circular status:", ssDNA.circular)
# except Exception as e:
#     print("Error setting circular:", e)

# # Test removing circular
# try:
#     ssDNA.remove_circular()
#     print("Circular removed successfully, circular status:", ssDNA.circular)
# except Exception as e:
#     print("Error removing circular:", e)

# print(ssDNA)

# ssDNA.reverse_complement()
# ssDNA.add_annotations([
#     {"name": "Annot1", "start": 0, "end": 3, "note": "This is Gene1"},
# ])
# print(ssDNA)

# ssDNA.annotations["wee"] = "TEST TEST TEST"

# print("test after", ssDNA)
