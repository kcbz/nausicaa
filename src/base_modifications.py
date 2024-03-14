base_modification_types = {
    "5-mC": "C",
    "5-hmC": "C",
    "5-fC": "C",
    "5-caC": "C",
    "6-mA": "A",
}


class BaseModification:

    def __init__(self, parent_ss_nucleic_acid, position, modification_type):
        self.validate_base_modification(parent_ss_nucleic_acid, position, modification_type)

        self._parent_ss_nucleic_acid_id = parent_ss_nucleic_acid.id
        self._position = position
        self._base = parent_ss_nucleic_acid.sequence[position]
        self._modification_type = modification_type

    def __repr__(self):
        return (f"BaseModification(position={self.position}, base='{self.base}', "
                f"modification_type='{self.modification_type}')")

    @property
    def parent_ss_nucleic_acid_id(self):
        return self._parent_ss_nucleic_acid_id

    @property
    def position(self):
        return self._position

    @property
    def modification_type(self):
        return self._modification_type

    @property
    def base(self):
        return self._base

    def validate_base_modification(self, parent_ss_nucleic_acid, position, modification_type):
        if not isinstance(position, int):
            raise TypeError("Cannot add base modification, position must be an integer.")
        if not isinstance(modification_type, str):
            raise TypeError("Cannot add base modification, modification_type must be a string.")
        if modification_type not in base_modification_types.keys():
            raise ValueError("Cannot add base modification, invalid modification type. "
                             f"Options are {base_modification_types.keys()}")
        base = parent_ss_nucleic_acid.sequence[position]
        if base != base_modification_types[modification_type]:
            raise ValueError(
                "Cannot add modification, invalid base for specified modification type."
                f"Base for {modification_type} should be "
                f"{base_modification_types[modification_type]}")
        if position < 0 or position >= len(parent_ss_nucleic_acid.sequence):
            raise ValueError(
                "Cannot add modification, must be within index bounds of the parent "
                "nucleic acid sequence.")
