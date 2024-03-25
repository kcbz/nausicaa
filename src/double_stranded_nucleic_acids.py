import uuid
from dna_features_viewer import CircularGraphicRecord, GraphicFeature, GraphicRecord

from src.nucleic_acids import NucleicAcidSequence, complement_sequence
from src.single_stranded_nucleic_acids import SingleStrandNucleicAcidSequence
from src.util_classes import NucleicAcidTypes, StrandDirections, Colors


class DoubleStrandNucleicAcidSequence(NucleicAcidSequence):

    def __init__(
        self,
        forward_sequence,
        nucleic_acid_type=NucleicAcidTypes.DNA.value,
        circular=False,
        note="",
        reverse_sequence=None,
        reverse_sequence_start=None,
    ):
        self._validate_nucleic_acid_type_with_sequences(forward_sequence, reverse_sequence, nucleic_acid_type)
        self._validate_circular_with_sequences(forward_sequence, reverse_sequence, reverse_sequence_start, circular)

        super().__init__(nucleic_acid_type, circular)
        self._id = str(uuid.uuid4())
        self.note = note
        (self._forward_sequence,
         self._reverse_sequence,
         self._reverse_sequence_start) = self._validate_and_set_sequences(forward_sequence,
                                                                          reverse_sequence,
                                                                          reverse_sequence_start,
                                                                          nucleic_acid_type,
                                                                          circular)

    def __repr__(self):
        return (f"DoubleStrandNucleicAcidSequence(id='{self.id}', "
                f"forward_sequence='{self.forward_sequence}', "
                f"reverse_sequence='{self.reverse_sequence}', "
                f"reverse_sequence_start='{self.reverse_sequence_start}', "
                f"nucleic_acid_type='{self.nucleic_acid_type}', "
                f"circular={self.circular}, "
                f"note='{self.note}')")

    @property
    def id(self):
        return self._id

    @property
    def forward_sequence(self):
        return self._forward_sequence

    @property
    def reverse_sequence(self):
        return self._reverse_sequence

    @property
    def reverse_sequence_start(self):
        return self._reverse_sequence_start

    @property
    def nucleic_acid_type(self):
        return self._nucleic_acid_type

    @property
    def circular(self):
        return self._circular

    @property
    def note(self):
        return self._note

    @note.setter
    def note(self, value):
        self.validate_note(value)
        self._note = value

    def _validate_nucleic_acid_type_with_sequences(self, fwd_seq, rev_seq, nuc_acid_type):
        error_types = []
        if isinstance(fwd_seq, SingleStrandNucleicAcidSequence) and fwd_seq.nucleic_acid_type != nuc_acid_type:
            error_types.append("forward")
        if rev_seq:
            if isinstance(rev_seq, SingleStrandNucleicAcidSequence) and rev_seq.nucleic_acid_type != nuc_acid_type:
                error_types.append("reverse")
        if error_types:
            sequence_str = " and ".join(error_types)
            raise ValueError(
                f"{sequence_str} sequence nucleic acid type(s) don't equal the specified "
                "nucleic acid type for the DoubleStrandNucleicAcidSequence object.")

    def _validate_circular_with_sequences(self, fwd_seq, rev_seq, rev_seq_start, circ):
        error_types = []
        if isinstance(fwd_seq, SingleStrandNucleicAcidSequence):
            fwd_seq_len = len(fwd_seq.sequence)
            if fwd_seq.circular != circ:
                error_types.append("forward")
        else:
            fwd_seq_len = len(fwd_seq)
        if rev_seq:
            if isinstance(rev_seq, SingleStrandNucleicAcidSequence):
                rev_seq_len = len(rev_seq.sequence)
                if rev_seq.circular != circ:
                    error_types.append("reverse")
            else:
                rev_seq_len = len(rev_seq)
        else:
            rev_seq_len = fwd_seq_len
        if error_types:
            sequence_str = " and ".join(error_types)
            raise ValueError(
                f"{sequence_str} sequence circular value(s) don't equal the specified "
                "circular value for the DoubleStrandNucleicAcidSequence object.")
        if circ:
            if fwd_seq_len != rev_seq_len:
                raise ValueError(
                    "Cannot set circular to True when the reverse sequence length is not "
                    "equal to the forward sequence length.")
            if rev_seq_start != 0:
                raise ValueError(
                    "Cannot set circular to True when the reverse sequence start is not "
                    "0.")

    def _validate_and_set_sequences(self, fwd_seq, rev_seq, rev_seq_start, nuc_acid_type, circ):
        fwd_seq_obj = self._create_or_validate_strand(fwd_seq, StrandDirections.FWD_STRAND.value, nuc_acid_type, circ)
        rev_seq_obj, rev_seq_start_value = self._handle_reverse_sequence(rev_seq, rev_seq_start, nuc_acid_type, circ, fwd_seq_obj)
        fwd_seq_obj._set_is_part_of_dsDNA_true
        rev_seq_obj._set_is_part_of_dsDNA_true
        return fwd_seq_obj, rev_seq_obj, rev_seq_start_value

    def _handle_reverse_sequence(self, rev_seq, rev_seq_start, nuc_acid_type, circ, fwd_seq_obj):
        if not rev_seq and rev_seq_start:
            raise ValueError("Cannot set reverse sequence start if no reverse sequence provided.")
        if rev_seq:
            rev_seq_obj = self._create_or_validate_strand(rev_seq, StrandDirections.REV_STRAND.value, nuc_acid_type, circ)
            rev_seq_start_value = self._validate_rev_seq_start(rev_seq_start, fwd_seq_obj, rev_seq_obj)
        else:
            rev_strand_seq = complement_sequence(fwd_seq_obj.sequence)
            rev_seq_obj = SingleStrandNucleicAcidSequence(
                sequence=rev_strand_seq,
                nucleic_acid_type=nuc_acid_type,
                circular=circ,
                strand_direction=StrandDirections.REV_STRAND.value
            )
            rev_seq_start_value = 0
        return rev_seq_obj, rev_seq_start_value

    def _create_or_validate_strand(self, seq, expected_strand_direction, nuc_acid_type, circ):
        if isinstance(seq, SingleStrandNucleicAcidSequence):
            self._validate_ssDNA_attributes(seq, expected_strand_direction)
            return seq
        else:
            return SingleStrandNucleicAcidSequence(
                sequence=seq,
                nucleic_acid_type=nuc_acid_type,
                circular=circ,
                strand_direction=expected_strand_direction
            )

    def _validate_ssDNA_attributes(self, ssDNA, strand_direction):
        if ssDNA.strand_direction != strand_direction:
            expected_direction = "forward" if strand_direction == StrandDirections.FWD_STRAND.value else "reverse"
            raise ValueError(
                "Cannot set a SingleStrandNucleicAcidSequence with incorrect strand direction. "
                f"Expected {expected_direction}.")

    def _validate_rev_seq_start(self, rev_seq_start, fwd_seq_obj, rev_seq_obj):
        if rev_seq_start is not None:
            if not isinstance(rev_seq_start, int):
                raise ValueError("Reverse sequence start must be an integer.")
            if rev_seq_start >= len(fwd_seq_obj.sequence):
                raise ValueError(
                    "Cannot set reverse sequence start to an integer value greater than "
                    "the length of the forward sequence.")
            if rev_seq_start <= -1 * len(rev_seq_obj.sequence):
                raise ValueError(
                    "Cannot set reverse sequence start to an integer value less than "
                    "the negative length of the reverse sequence.")
        else:
            rev_seq_start = 0
        return rev_seq_start

    def validate_note(self, note):
        if not isinstance(note, str):
            raise TypeError("Cannot make sequence, note must be a string.")
    
    def copy(self):
        new_ds_seq = DoubleStrandNucleicAcidSequence(
            forward_sequence=self.forward_sequence,
            reverse_sequence=self.reverse_sequence,
            rev_sequence_start=self.reverse_sequence_start,
            nucleic_acid_type=self.nucleic_acid_type,
            circular=self.circular,
            note=self.note,
        )
        return new_ds_seq

    def reverse(self):
        with self._forward_sequence._unlocked(), self._reverse_sequence._unlocked():
            self._forward_sequence.reverse()
            self._reverse_sequence.reverse()
            self._reverse_sequence_start = -1 * (self.reverse_sequence_start - (
                len(self.forward_sequence.sequence) - len(self.reverse_sequence.sequence)))

    def complement(self):
        with self._forward_sequence._unlocked(), self._reverse_sequence._unlocked():
            self._forward_sequence.complement()
            self._reverse_sequence.complement()

    def reverse_complement(self):
        with self._forward_sequence._unlocked(), self._reverse_sequence._unlocked():
            self._forward_sequence.reverse_complement()
            self._reverse_sequence.reverse_complement()
            self._reverse_sequence_start = -1 * (self.reverse_sequence_start - (
                len(self.forward_sequence.sequence) - len(self.reverse_sequence.sequence)))

    def set_circular(self):
        with self._forward_sequence._unlocked(), self._reverse_sequence._unlocked():
            if self.circular is True:
                raise Exception("circular is already True.")
            else:
                self._validate_circular_with_sequences(
                    self.forward_sequence,
                    self.reverse_sequence,
                    self.reverse_sequence_start,
                    self.circular)
                self._forward_sequence.set_circular()
                self._reverse_sequence.set_circular()
                self._circular = True

    def remove_circular(self):
        with self._forward_sequence._unlocked(), self._reverse_sequence._unlocked():
            if self.circular is False:
                raise Exception("circular is already False.")
            else:
                self._forward_sequence.remove_circular()
                self._reverse_sequence.remove_circular()
                self._circular = False

    def view(self):
        features = []
        for annotation in self.forward_sequence.annotations.values():
            label=annotation.name
            if annotation.note:
                label=f"{annotation.name}, {annotation.note}"
            features.append(GraphicFeature(
                label=label,
                start=annotation.start,
                end=annotation.end,
                strand=+1,
                color=Colors.random_color(),
            ))
        for annotation in self.reverse_sequence.annotations.values():
            label=annotation.name
            if annotation.note:
                label=f"{annotation.name}, {annotation.note}"
            features.append(GraphicFeature(
                label=label,
                start=annotation.start,
                end=annotation.end,
                strand=-1,
                color=Colors.random_color(),
            ))
        seq_length = len(self.forward_sequence.sequence)
        if seq_length < len(self.reverse_sequence.sequence):
            seq_length = len(self.reverse_sequence.sequence)
        if self.circular:
            record = CircularGraphicRecord(sequence_length=seq_length, features=features)
        else:
            record = GraphicRecord(sequence_length=seq_length, features=features)
        record.plot(figure_width=5)