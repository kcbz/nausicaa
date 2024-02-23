from src.single_stranded_nucleic_acids import SingleStrandNucleicAcidSequence


# Sample test cases

# Valid DNA sequence
try:
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
    print("Valid DNA Sequence Created:", dna_seq)
except Exception as e:
    print("Error creating DNA sequence:", e)

# Valid RNA sequence
try:
    rna_seq = SingleStrandNucleicAcidSequence(sequence="AUGCGUAA", nucleic_acid_type="RNA", circular=True, strand_direction="forward")
    print("Valid RNA Sequence Created:", rna_seq)
except Exception as e:
    print("Error creating RNA sequence:", e)

# Invalid sequence characters
try:
    invalid_seq = SingleStrandNucleicAcidSequence(sequence="ATBZ", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
    print("Invalid Sequence Created:", invalid_seq)
except Exception as e:
    print("Error creating invalid sequence:", e)

# DNA sequence with 'U'
try:
    dna_with_u = SingleStrandNucleicAcidSequence(sequence="ATGCUA", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
    print("DNA Sequence with U Created:", dna_with_u)
except Exception as e:
    print("Error creating DNA sequence with U:", e)

# RNA sequence with 'T'
try:
    rna_with_t = SingleStrandNucleicAcidSequence(sequence="AUGCTA", nucleic_acid_type="RNA", circular=False, strand_direction="forward")
    print("RNA Sequence with T Created:", rna_with_t)
except Exception as e:
    print("Error creating RNA sequence with T:", e)

try:
    ssDNA = SingleStrandNucleicAcidSequence(sequence="ATGCGNNNNNTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
    ssDNA.add_annotations([
        {"name": "Gene1", "start": 0, "end": 3, "note": "This is Gene1"},
        {"name": "random_barcode", "start": 5, "end": 9, "note": "This a random barcode"},
    ])
    print("Annotations added successfully:", ssDNA.annotations)
except Exception as e:
    print("Error adding annotations:", e)

# Test editing an annotation
try:
    ssDNA.edit_annotation("Gene1", start=1, end=4)
    print("Annotation edited successfully:", ssDNA.annotations['Gene1'])
except Exception as e:
    print("Error editing annotation:", e)

# Test removing an annotation
try:
    ssDNA.remove_annotations(["Gene1"])
    print("Annotation removed successfully, remaining annotations:", ssDNA.annotations)
except Exception as e:
    print("Error removing annotation:", e)

# Test adding base modifications
try:
    ssDNA.add_base_modifications([
        {"position": 3, "modification_type": "5-mC"}
    ])
    print("Base modifications added successfully:", ssDNA.base_modifications)
except Exception as e:
    print("Error adding base modifications:", e)

# Test editing a base modification
try:
    ssDNA.edit_base_modification(3, modification_type="5-hmC")
    print("Base modification edited successfully:", ssDNA.base_modifications[3])
except Exception as e:
    print("Error editing base modification:", e)

# Test removing a base modification
try:
    ssDNA.remove_base_modifications([3])
    print("Base modification removed successfully, remaining:", ssDNA.base_modifications)
except Exception as e:
    print("Error removing base modification:", e)

# Test reverse complement
try:
    print("Original sequence:", ssDNA.sequence)
    ssDNA.reverse_complement()
    print("Reverse complement successful, new sequence:", ssDNA.sequence)
except Exception as e:
    print("Error computing reverse complement:", e)

# Test setting circular
try:
    ssDNA.set_circular()
    print("Set circular successful, circular status:", ssDNA.circular)
except Exception as e:
    print("Error setting circular:", e)

# Test removing circular
try:
    ssDNA.remove_circular()
    print("Circular removed successfully, circular status:", ssDNA.circular)
except Exception as e:
    print("Error removing circular:", e)

print(ssDNA)

ssDNA.reverse_complement()
ssDNA.add_annotations([
    {"name": "Annot1", "start": 0, "end": 3, "note": "This is Gene1"},
])
print(ssDNA)

ssDNA.annotations["wee"] = "TEST TEST TEST"

print("test after", ssDNA)
