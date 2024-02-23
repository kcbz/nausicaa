from src.double_stranded_nucleic_acids import DoubleStrandNucleicAcidSequence
from src.single_stranded_nucleic_acids import SingleStrandNucleicAcidSequence


# Create a forward single-stranded DNA sequence
forward_ssDNA = SingleStrandNucleicAcidSequence(
    sequence="ATGCGTAA", 
    nucleic_acid_type="DNA", 
    circular=False, 
    strand_direction="forward"
)

# Create a reverse single-stranded DNA sequence manually
reverse_ssDNA = SingleStrandNucleicAcidSequence(
    sequence="TACGCATTAAAAA",
    nucleic_acid_type="DNA",
    circular=False,
    strand_direction="reverse"
)

# Create a DoubleStrandNucleicAcidSequence using the forward sequence and a string for the reverse sequence
dsDNA = DoubleStrandNucleicAcidSequence(
    forward_sequence=forward_ssDNA, 
    nucleic_acid_type="DNA",
    circular=False,
    note="Example DS DNA",
)
print("Original Double-Stranded DNA:\n", dsDNA)

# Reverse complement the double-stranded DNA
dsDNA.reverse_complement()
print("\nReverse Complemented Double-Stranded DNA:\n", dsDNA)

# Set the double-stranded DNA to circular
dsDNA.set_circular()
print("\nCircular Double-Stranded DNA:\n", dsDNA)

# Remove circularity from the double-stranded DNA
dsDNA.remove_circular()
print("\nLinear Double-Stranded DNA:\n", dsDNA)

dsDNA2 = DoubleStrandNucleicAcidSequence(
    forward_sequence=forward_ssDNA,
    reverse_sequence=reverse_ssDNA,
    nucleic_acid_type="DNA",
    circular=False,
    note="Example DS DNA",
)
print("Original Double-Stranded DNA2:\n", dsDNA2)

# Reverse complement the double-stranded DNA
dsDNA2.reverse_complement()
print("\nReverse Complemented Double-Stranded DNA2:\n", dsDNA2)

dsDNA3 = DoubleStrandNucleicAcidSequence(
    forward_sequence="ATGCTAAA",
    # reverse_sequence="CCCCTACGATTTTTT",
    nucleic_acid_type="DNA",
    circular=False,
    reverse_sequence_start=-4,
    note="Example DS DNA",
)
print("Original Double-Stranded DNA3:\n", dsDNA3)
