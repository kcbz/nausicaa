class RestrictionEnzymes:

    EcoRI = {"recognition_sequence": "GAATTC", "cut_position": 1}  # Cuts between G and A
    HindIII = {"recognition_sequence": "AAGCTT", "cut_position": 1}  # Cuts between A and A
    BamHI = {"recognition_sequence": "GGATCC", "cut_position": 1}  # Cuts between G and G
    NotI = {"recognition_sequence": "GCGGCCGC", "cut_position": 2}  # Cuts between G and C
    XhoI = {"recognition_sequence": "CTCGAG", "cut_position": 1}  # Cuts between C and T
    PstI = {"recognition_sequence": "CTGCAG", "cut_position": 5}  # Cuts after the G
    SmaI = {"recognition_sequence": "CCCGGG", "cut_position": 3}  # Cuts in the middle
    KpnI = {"recognition_sequence": "GGTACC", "cut_position": 5}  # Cuts after the first G
    SacI = {"recognition_sequence": "GAGCTC", "cut_position": 5}  # Cuts after the G
    SpeI = {"recognition_sequence": "ACTAGT", "cut_position": 1}  # Cuts between A and C
    NcoI = {"recognition_sequence": "CCATGG", "cut_position": 1}  # Cuts between C and C
    SalI = {"recognition_sequence": "GTCGAC", "cut_position": 1}  # Cuts between G and T
    XbaI = {"recognition_sequence": "TCTAGA", "cut_position": 1}  # Cuts between T and C
    SphI = {"recognition_sequence": "GCATGC", "cut_position": 5}  # Cuts after the G

    def __init__(self, restriction_enzyme, recognition_sequence, cut_position):
        self.restriction_enzyme = restriction_enzyme
        self.recognition_sequence = recognition_sequence
        self.cut_position = cut_position
    
    def __repr__(self):
        return (f"RestrictionEnzyme(restriction_enzyme='{self.restriction_enzyme}', "
                f"recognition_sequence='{self.recognition_sequence}', "
                f"cut_position={self.cut_position})")

    @classmethod
    def list_names(cls):
        return [member.name for member in cls]

    @classmethod
    def list_info(cls):
        return [{member.name: {"recognition_sequence": member.value["recognition_sequence"], "cut_position": member.value["cut_position"]}} for member in cls]

    @classmethod
    def get_info(cls, enzyme_name):
        # Get the recognition sequence and cut position of the given enzyme.
        # If the enzyme name is not found, return None.
        return cls[enzyme_name].value if enzyme_name in cls._member_names_ else None
