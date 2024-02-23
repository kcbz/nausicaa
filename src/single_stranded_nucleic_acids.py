import logging
import uuid
from contextlib import contextmanager
from dna_features_viewer import CircularGraphicRecord, GraphicFeature, GraphicRecord
from typing import Optional

from src.base_modifications import BaseModification
from src.nucleic_acids import NucleicAcidSequence, reverse_sequence, complement_sequence
from src.sequence_annotation import SequenceAnnotation
from src.util_classes import IUPACCodes, NucleicAcidTypes, StrandDirections, Colors

logging.basicConfig(level=logging.INFO)


class SingleStrandNucleicAcidSequence(NucleicAcidSequence):

    def __repr__(self):
        return (f"SingleStrandNucleicAcidSequence(id='{self.id}', "
                f"is_part_of_dsDNA='{self._is_part_of_dsDNA}', "
                f"sequence='{self.sequence}', "
                f"nucleic_acid_type='{self.nucleic_acid_type}', "
                f"circular={self.circular}, "
                f"strand_direction='{self.strand_direction}', "
                f"annotations={self.annotations}, "
                f"base_modifications={self.base_modifications}, "
                f"note='{self.note}')")

    def __init__(
        self,
        sequence,
        nucleic_acid_type=NucleicAcidTypes.DNA.value,
        circular=False,
        strand_direction=StrandDirections.FWD_STRAND.value,
        note="",
        annotations=None,
        base_modifications=None,
    ):
        self._validate_sequence(sequence.upper())
        self._validate_nucleic_acid_type_with_sequence(nucleic_acid_type, sequence)
        self._validate_strand_direction(strand_direction)
        self._validate_note(note)

        super().__init__(nucleic_acid_type, circular)
        self._id = uuid.uuid4().hex
        self._is_part_of_dsDNA = False
        self._sequence = sequence.upper()
        self._strand_direction = strand_direction
        self.note = note
        self._annotations = {}
        self._base_modifications = {}

        if annotations:
            self.add_annotations(annotations)
        if base_modifications:
            self.add_base_modifications(base_modifications)

    @property
    def id(self):
        return self._id

    @property
    def is_part_of_dsDNA(self):
        return self._is_part_of_dsDNA

    @property
    def sequence(self):
        return self._sequence

    @property
    def nucleic_acid_type(self):
        return self._nucleic_acid_type

    @property
    def circular(self):
        return self._circular

    @property
    def strand_direction(self):
        return self._strand_direction

    @property
    def note(self):
        return self._note

    @note.setter
    def note(self, value):
        self._validate_note(value)
        self._note = value

    @property
    def annotations(self):
        return self._annotations.copy()

    @property
    def base_modifications(self):
        return self._base_modifications.copy()

    @contextmanager
    def _unlocked(self):
        original_state = self._is_part_of_dsDNA
        self._is_part_of_dsDNA = False
        try:
            yield
        finally:
            self._is_part_of_dsDNA = original_state

    def _validate_sequence(self, sequence):
        if len(sequence) == 0:
            raise ValueError("Cannot make sequence from an empty string.")
        if not all(nuc in IUPACCodes.list_names() for nuc in sequence.upper()):
            raise ValueError(
                f"Cannot make sequence, invalid nucleotide found in sequence: {sequence}"
            )

    def _validate_nucleic_acid_type_with_sequence(self, nucleic_acid_type, sequence):
        if nucleic_acid_type == NucleicAcidTypes.DNA.value and "U" in sequence:
            raise ValueError(
                f"Cannot be nucleic acid type {NucleicAcidTypes.DNA.value} "
                "when U in sequence.")
        elif nucleic_acid_type == NucleicAcidTypes.RNA.value and "T" in sequence:
            raise ValueError(
                f"Cannot be nucleic acid type {NucleicAcidTypes.RNA.value} "
                "when T in sequence.")

    def _validate_strand_direction(self, strand_direction):
        if strand_direction not in StrandDirections.list_values():
            raise ValueError(
                f"Invalid strand direction: {strand_direction}. "
                f"Options are {StrandDirections.list_values()}")

    def _validate_note(self, note):
        if not isinstance(note, str):
            raise TypeError("Note must be a string.")

    def _validate_annotations(self, annot_list):
        annot_names_set = set()
        for annot in annot_list:
            name = annot["name"]
            if name in annot_names_set:
                raise ValueError(
                    "Cannot add annotation, multiple annotations share the same name "
                    f"{name} in the provided annotations list.")
            annot_names_set.add(name)
            if name in self.annotations:
                raise ValueError(
                    f"Cannot add annotation {name}, annotation with the same name "
                    " already exists.")

    def _validate_base_modifications(self, base_mod_list):
        base_mod_pos_set = set()
        for base_mod in base_mod_list:
            pos = base_mod["position"]
            if pos in base_mod_pos_set:
                raise ValueError(
                    "Cannot add base modification, multiple base modifications share the "
                    f"same position {pos} in the provided base modification list.")
            base_mod_pos_set.add(pos)
            if pos in self.base_modifications:
                raise ValueError(
                    "Cannot add base modification, base modification already exists at "
                    f"position {pos}.")

    def add_annotations(self, annot_list):
        self._validate_annotations(annot_list)
        for annot in annot_list:
            name = annot["name"]
            start = annot["start"]
            end = annot["end"]
            note = annot.get("note", "")
            new_annot = SequenceAnnotation(self, name, start, end, note)
            self._annotations[name] = new_annot

    def edit_annotation(self, annot_name, **kwargs):
        if annot_name not in self.annotations:
            raise KeyError(f"Annotation with name {annot_name} does not exist.")
        current_annot = self.annotations[annot_name]
        new_name = kwargs.get("name", current_annot.name)
        new_start = kwargs.get("start", current_annot.start)
        new_end = kwargs.get("end", current_annot.end)
        new_note = kwargs.get("note", current_annot.note)
        if ("name" in kwargs) and (new_name != annot_name) and (new_name in self.annotations):
            raise ValueError(f"Cannot change annotation name to {new_name} as it already exists.")
        valid_keys = ["name", "start", "end", "note"]
        for key in kwargs:
            if key not in valid_keys:
                raise AttributeError(f"Annotation {annot_name} has no attribute {key}.")
        updated_annot = SequenceAnnotation(self, new_name, new_start, new_end, new_note)
        if new_name != annot_name:
           del self._annotations[annot_name]
        self._annotations[new_name] = updated_annot

    def remove_annotations(self, annot_names):
        for annot_name in annot_names:
            if annot_name in self.annotations:
                del self._annotations[annot_name]
            else:
                raise KeyError(f"Annotation with name {annot_name} does not exist.")

    def remove_all_annotations(self):
        annot_names = list(self.annotations.keys())
        self.remove_annotations(annot_names)

    def add_base_modifications(self, base_mod_list):
        self._validate_base_modifications(base_mod_list)
        for base_mod in base_mod_list:
            pos = base_mod["position"]
            mod_type = base_mod["modification_type"]
            new_mod = BaseModification(self, pos, mod_type)
            self._base_modifications[pos] = new_mod

    def edit_base_modification(self, base_mod_pos, **kwargs):
        if base_mod_pos not in self.base_modifications:
            raise KeyError(f"Base modification with position {base_mod_pos} does not exist.")
        current_base_mod = self.base_modifications[base_mod_pos]
        new_pos = kwargs.get("position", current_base_mod.position)
        new_mod_type = kwargs.get("modification_type", current_base_mod._modification_type)
        if ("position" in kwargs) and (new_pos != base_mod_pos) and (new_pos in self.base_modifications):
            raise ValueError(
                f"Cannot change base modification position to {new_pos} as it already exists.")
        valid_keys = ["position", "modification_type"]
        for key in kwargs:
            if key == "base":
                raise AttributeError("Unable to change the base attribute directly. "
                                     f"Must use {valid_keys}.")
            elif key not in valid_keys:
                raise AttributeError(f"Base modification at {base_mod_pos} has no attribute {key}.")
        updated_base_mod = BaseModification(self, new_pos, new_mod_type)
        if new_pos != base_mod_pos:
            del self._base_modifications[base_mod_pos]
        self._base_modifications[new_pos] = updated_base_mod

    def remove_base_modifications(self, base_mod_pos_list):
        for base_mod_pos in base_mod_pos_list:
            if base_mod_pos in self.base_modifications:
                del self._base_modifications[base_mod_pos]
            else:
                raise KeyError(f"Base modification with position {base_mod_pos} does not exist.")

    def remove_all_base_modifications(self):
        base_mod_pos_list = list(self.base_modifications.keys())
        self.remove_base_modifications(base_mod_pos_list)

    def copy(self):
        new_ss_seq = SingleStrandNucleicAcidSequence(
            sequence=self.sequence,
            nucleic_acid_type=self.nucleic_acid_type,
            strand_direction=self.strand_direction,
            circular=self.circular,
            note=self.note,
            annotations=self.annotations,
            base_modifications=self.base_modifications,
        )
        return new_ss_seq

    def change_strand_direction(self):
        if self._is_part_of_dsDNA:
            raise Exception("Operation not allowed directly on ssDNA part of dsDNA. "
                            "Use dsDNA methods.")
        if self.strand_direction == StrandDirections.FWD_STRAND.value:
                self._strand_direction = StrandDirections.REV_STRAND.value
        else:
            self._strand_direction = StrandDirections.FWD_STRAND.value

    def reverse(self, change_strand_dir=False):
        if self._is_part_of_dsDNA:
            raise Exception("Operation not allowed directly on ssDNA part of dsDNA. "
                            "Use dsDNA methods.")
        if change_strand_dir:
            self.change_strand_direction()
        self.remove_all_annotations()
        self.remove_all_base_modifications()
        self._sequence = reverse_sequence(self.sequence)

    def complement(self, change_strand_dir=False):
        if self._is_part_of_dsDNA:
            raise Exception("Operation not allowed directly on ssDNA part of dsDNA. "
                            "Use dsDNA methods.")
        if change_strand_dir:
            self.change_strand_direction()
        self.remove_all_annotations()
        self.remove_all_base_modifications()
        self._sequence = complement_sequence(self._sequence, nuc_type=self.nucleic_acid_type)

    def reverse_complement(self, change_strand_dir=False):
        self.reverse(change_strand_dir=change_strand_dir)
        self.complement(change_strand_dir=change_strand_dir)

    def set_circular(self):
        if self._is_part_of_dsDNA:
            raise Exception("Operation not allowed directly on ssDNA part of dsDNA. "
                            "Use dsDNA methods.")
        if self.circular is True:
            logging.info("circular is already True")
        else:
            self._circular = True

    def remove_circular(self):
        if self._is_part_of_dsDNA:
            raise Exception("Operation not allowed directly on ssDNA part of dsDNA. "
                            "Use dsDNA methods.")
        if self.circular is False:
            logging.info("circular is already False")
        else:
            self._circular = False
        annot_names_to_remove = [
            name for name, annot in self.annotations.items()
            if annot.end < annot.start
        ]
        self.remove_annotations(annot_names_to_remove)

    def view(self):
        strand = +1
        if self.strand_direction == StrandDirections.REV_STRAND.value:
            strand = -1
        features = []
        for annotation in self.annotations.values():
            label=annotation.name
            if annotation.note:
                label=f"{annotation.name}, {annotation.note}"
            features.append(GraphicFeature(
                label=label,
                start=annotation.start,
                end=annotation.end,
                strand=strand,
                color=Colors.random_color(),
            ))
        if self.circular:
            record = CircularGraphicRecord(sequence_length=len(self.sequence), features=features)
        else:
            record = GraphicRecord(sequence_length=len(self.sequence), features=features)
        record.plot(figure_width=5)


def transfer_annotations(
    source_annotation_dict: dict,
    target_strand: SingleStrandNucleicAcidSequence,
    higher_bound_index: Optional[int] = None,
    lower_bound_index: int = 0,
    index_adjustment: int = 0) -> SingleStrandNucleicAcidSequence:
    for annotation in source_annotation_dict.values():
        if not higher_bound_index:
            higher_bound_index = len(target_strand.sequence)
        if lower_bound_index <= annotation.start and annotation.end <= higher_bound_index:
            new_annot_start = annotation.start + index_adjustment
            new_annot_end = annotation.end + index_adjustment
            target_strand.add_annotations([{"name": annotation.name,
                                            "start": new_annot_start,
                                            "end": new_annot_end,
                                            "note": annotation.note}])
    return target_strand
