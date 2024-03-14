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

class RestrictionEnzymeCutSite:

    def __init__(self, parent_ss_nucleic_acid, start, restriction_enzyme_type):
        cut_site_data = self.validate_restriction_site(parent_ss_nucleic_acid, start, restriction_enzyme_type)

        self._parent_nucleic_acid_sequence_id = parent_ss_nucleic_acid.id
        self._start = start
        self._end = start + len(cut_site_data["recognition_sequence"]) - 1
        self._cut_position = start + cut_site_data["cut_position"]
        self._recognition_sequence = cut_site_data["recognition_sequence"]
        self._restriction_enzyme = restriction_enzyme_type

    def __repr__(self):
        return (f"RestrictionEnzymeCutSite(start='{self.start}, end='{self.end}', "
                f"cut_position='{self.cut_position}', "
                f"recognition_sequence='{self.recognition_sequence}', "
                f"restriction_enzyme='{self.restriction_enzyme}')")

    @property
    def parent_nucleic_acid_sequence_id(self):
        return self._parent_nucleic_acid_sequence_id

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def cut_position(self):
        return self._cut_position

    @property
    def recognition_sequence(self):
        return self._recognition_sequence

    @property
    def restriction_enzyme(self):
        return self._restriction_enzyme

    def validate_restriction_site(self, parent_ss_nucleic_acid, start, restriction_enzyme_type):
        if not isinstance(start, int):
            raise TypeError("Cannot add restriction site, start must be an integer.")
        if not isinstance(restriction_enzyme_type, str):
            raise TypeError("Cannot add restriction site, restriction_enzyme_type must be a string.")
        if restriction_enzyme_type not in restriction_enzyme_types.keys():
            raise ValueError("Cannot add restriction site, invalid restriction site type. "
                             f"Options are {restriction_enzyme_types.keys()}")
        cut_site_data = restriction_enzyme_types[restriction_enzyme_type]
        cut_site_sequence = cut_site_data["recognition_sequence"]
        found_cut_site_sequence = parent_ss_nucleic_acid.sequence[start:start + len(cut_site_sequence)]

        print(start)
        print(start + len(cut_site_sequence))
        print(found_cut_site_sequence)

        if found_cut_site_sequence != cut_site_sequence:
            raise TypeError(f"Cannot add restriction site, sequence for {restriction_enzyme_type}"
                            f"is {cut_site_sequence}, found {found_cut_site_sequence}.")
        return cut_site_data
