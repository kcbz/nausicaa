from src.util_classes import NucleicAcidTypes


class NucleicAcidSequence:

    def __init__(self, nucleic_acid_type, circular):
        self.validate_nucleic_acid_type(nucleic_acid_type)
        self.validate_circular(circular)
        self._nucleic_acid_type = nucleic_acid_type
        self._circular = circular

    @property
    def nucleic_acid_type(self):
        return self._nucleic_acid_type

    @property
    def circular(self):
        return self._circular

    @staticmethod
    def validate_nucleic_acid_type(nucleic_acid_type):
        if nucleic_acid_type not in NucleicAcidTypes.list_values():
            raise ValueError(f"Invalid nucleic acid type: {nucleic_acid_type}. "
                             f"Options are {NucleicAcidTypes.list_values()}.")

    @staticmethod
    def validate_circular(circular):
        if not isinstance(circular, bool):
            raise ValueError("Circular value must be boolean.")


def reverse_sequence(input_sequence: str) -> str:
    return "".join(reversed(input_sequence))

def complement_sequence(input_sequence: str, nuc_type: str = "DNA") -> str:
    if nuc_type not in (NucleicAcidTypes.DNA.value, NucleicAcidTypes.RNA.value):
        raise ValueError(f"Incorrect nuc type, must be DNA or RNA, not {nuc_type}.")
    comp_map = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "U": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
    if nuc_type == "RNA":
        comp_map["A"] = "U"
    return "".join(comp_map[nuc] for nuc in input_sequence)

def reverse_complement_sequence(input_sequence: str, nuc_type: str = "DNA") -> str:
    rev_seq = reverse_sequence(input_sequence)
    return complement_sequence(rev_seq, nuc_type=nuc_type)
