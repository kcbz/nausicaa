class RestrictionEnzymeSite:
    restriction_enzyme_types = {
        "EcoRI": {"recognition_sequence": "GAATTC", "cut_position": 1},  # Cuts between G and A
        "HindIII": {"recognition_sequence": "AAGCTT", "cut_position": 1},  # Cuts between A and A
        "BamHI": {"recognition_sequence": "GGATCC", "cut_position": 1},  # Cuts between G and G
        "NotI": {"recognition_sequence": "GCGGCCGC", "cut_position": 2},  # Cuts between G and C
        "XhoI": {"recognition_sequence": "CTCGAG", "cut_position": 1},  # Cuts between C and T
        "PstI": {"recognition_sequence": "CTGCAG", "cut_position": 5},  # Cuts after the G
        "SmaI": {"recognition_sequence": "CCCGGG", "cut_position": 3},  # Cuts in the middle
        "KpnI": {"recognition_sequence": "GGTACC", "cut_position": 5},  # Cuts after the first G
        "SacI": {"recognition_sequence": "GAGCTC", "cut_position": 5},  # Cuts after the G
        "SpeI": {"recognition_sequence": "ACTAGT", "cut_position": 1},  # Cuts between A and C
        "NcoI": {"recognition_sequence": "CCATGG", "cut_position": 1},  # Cuts between C and C
        "SalI": {"recognition_sequence": "GTCGAC", "cut_position": 1},  # Cuts between G and T
        "XbaI": {"recognition_sequence": "TCTAGA", "cut_position": 1},  # Cuts between T and C
        "SphI": {"recognition_sequence": "GCATGC", "cut_position": 5}  # Cuts after the G
    }

    def __init__(self, parent_nucleic_acid_sequence, restriction_enzyme):
        enzyme_data = self.validate_restriction_site(self, restriction_enzyme)


        self._parent_nucleic_acid_sequence_id = parent_nucleic_acid_sequence.id
        self.restriction_enzyme = restriction_enzyme
        self.start = ...
        self.end = ...
        # self.cut_position = cut_position

    def __repr__(self):
        return (f"RestrictionEnzyme(restriction_enzyme='{self.restriction_enzyme}', "
                f"recognition_sequence='{self.recognition_sequence}', "
                f"cut_position={self.cut_position})")

    @property
    def parent_nucleic_acid_sequence_id(self):
        return self._parent_nucleic_acid_sequence_id

    @property
    def restriction_enzyme(self):
        return self._restriction_enzyme

    @property
    def recognition_sequence(self):
        return self._recognition_sequence

    @property
    def cut_position(self):
        return self._cut_position

    def validate_restriction_site(self, restriction_enzyme):
        enzyme_data = self.restriction_enzyme_types.get(restriction_enzyme, None)
        if enzyme_data is None:
            raise ValueError(f"{restriction_enzyme} is not a recognized restriction enzyme.")
        return enzyme_data

    def set_restriction_site(self, enzyme_data, parent_nucleic_acid_sequence):
        rec_seq = enzyme_data["recognition_sequence"]
        cut_pos = enzyme_data["cut_position"]
        start = template_seq.index(trimmed_query_seq)
