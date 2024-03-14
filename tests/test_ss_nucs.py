from src.single_stranded_nucleic_acids import SingleStrandNucleicAcidSequence


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


def test_make_ss_dna_with_annotations_then_edit_annotation():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGNNNNNTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward",
                                              annotations=[{"name": "Gene1", "start": 0, "end": 3, "note": "This is Gene1"},
                                                            {"name": "random_barcode", "start": 5, "end": 9, "note": "This is a random barcode"}])
    dna_seq.edit_annotation("Gene1", name="Lac", end=4)

    assert dna_seq.annotations["Lac"].name == "Lac"
    assert dna_seq.annotations["Lac"].start == 0
    assert dna_seq.annotations["Lac"].end == 4
    assert dna_seq.annotations["Lac"].note == "This is Gene1"
    assert dna_seq.annotations["random_barcode"].name == "random_barcode"
    assert dna_seq.annotations["random_barcode"].start == 5
    assert dna_seq.annotations["random_barcode"].end == 9
    assert dna_seq.annotations["random_barcode"].note == "This is a random barcode"


def test_make_ss_dna_with_base_modifications_then_edit_base_modification():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGNNNNNTCCAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward",
                                              base_modifications=[{"position": 3, "modification_type": "5-mC"},
                                                                  {"position": 0, "modification_type": "6-mA"}])
    dna_seq.edit_base_modification(3, position=12, modification_type="5-hmC")

    assert dna_seq.base_modifications[12].position == 12
    assert dna_seq.base_modifications[12].modification_type == "5-hmC"
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


def test_reverse_sequence():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward",
                                              base_modifications=[{"position": 0, "modification_type": "6-mA"}])
    dna_seq.reverse()

    assert dna_seq.sequence == "AATGCGTA"
    assert dna_seq.base_modifications == {}


def test_compliment_sequence():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward",
                                              base_modifications=[{"position": 0, "modification_type": "6-mA"}])
    dna_seq.complement()

    assert dna_seq.sequence == "TACGCATT"
    assert dna_seq.base_modifications == {}


def test_reverse_compliment_sequence():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGTAA", nucleic_acid_type="DNA", circular=False, strand_direction="forward",
                                              base_modifications=[{"position": 0, "modification_type": "6-mA"}])
    dna_seq.reverse_complement()

    assert dna_seq.sequence == "TTACGCAT"
    assert dna_seq.base_modifications == {}


def test_change_circular():
    dna_seq = SingleStrandNucleicAcidSequence(sequence="ATGCGGAATTCTAGCATGCAAATT", nucleic_acid_type="DNA", circular=False, strand_direction="forward")
    dna_seq.set_circular()

    assert dna_seq.circular == True
