import random
from enum import Enum


class IUPACCodes(Enum):
    A = ["A"]  # Adenine
    C = ["C"]  # Cytosine
    G = ["G"]  # Guanine
    T = ["T"]  # Thymine
    U = ["U"]  # Uracil
    R = ["A", "G"]  # Purine (A or G)
    Y = ["C", "T"]  # Pyrimidine (C or T)
    S = ["G", "C"]  # Strong interaction (G or C)
    W = ["A", "T"]  # Weak interaction (A or T)
    K = ["G", "T"]  # Keto (G or T)
    M = ["A", "C"]  # Amino (A or C)
    B = ["C", "G", "T"]  # Not A (C or G or T)
    D = ["A", "G", "T"]  # Not C (A or G or T)
    H = ["A", "C", "T"]  # Not G (A or C or T)
    V = ["A", "C", "G"]  # Not T (A or C or G)
    N = ["A", "C", "G", "T"]  # Any (A or C or G or T)

    @classmethod
    def list_names(cls):
        return [member.name for member in cls]

    @classmethod
    def list_values(cls):
        return [member.value for member in cls]


class NucleicAcidTypes(Enum):
    DNA = "DNA"
    RNA = "RNA"

    @classmethod
    def list_names(cls):
        return [member.name for member in cls]

    @classmethod
    def list_values(cls):
        return [member.value for member in cls]


class StrandDirections(Enum):
    FWD_STRAND = "forward"
    REV_STRAND = "reverse"

    @classmethod
    def list_names(cls):
        return [member.name for member in cls]

    @classmethod
    def list_values(cls):
        return [member.value for member in cls]


class Colors(Enum):
    light_blue = "#7fdbff"
    blue = "#0074d9"
    green = "#2ecc40"
    yellow = "#ffdc00"
    orange = "#ff851b"
    neon_green = "#01ff70"
    purple = "#f012be"
    dark_purple = "#b10dc9"
    red = "#ff4136"
    maroon = "#85144b"
    teal = "#39cccc"
    olive = "#3d9970"
    navy = "#001f3f"
    grey = "#aaaaaa"
    white = "#ffffff"
    light_grey = "#dddddd"

    @classmethod
    def list_names(cls):
        return [member.name for member in cls]

    @classmethod
    def list_values(cls):
        return [member.value for member in cls]

    @classmethod
    def random_color(cls):
        return random.choice(list(cls.__members__.values())).value

    @classmethod
    def get_color(cls, code):
    # Get the color associated with the given IUPAC code.
    # If the code is not found, return a default color (e.g., black).
        return cls[code].value if code in cls._member_names_ else "#111111"
