import logging
import uuid

from src.nucleic_acids import reverse_sequence, complement_sequence
from src.double_stranded_nucleic_acids import DoubleStrandNucleicAcidSequence
from src.single_stranded_nucleic_acids import SingleStrandNucleicAcidSequence, transfer_annotations
from src.util_classes import NucleicAcidTypes, StrandDirections
from src.util_functions import kmer_trimming_search

logging.basicConfig(level=logging.INFO)


class NucleicAcidReaction:
    def __init__(self, reaction_type, inputs, **kwargs):
        self.id = uuid.uuid4()
        self.reaction_type = reaction_type
        self.inputs = inputs
        self.kwargs = kwargs
        self.outputs = {}
        self.execute_reaction()

    def execute_reaction(self):
        """
        Dispatches the reaction execution to the appropriate method
        based on the reaction type.
        """
        method_name = f"_perform_{self.reaction_type}"
        method = getattr(self, method_name, None)
        if callable(method):
            self.outputs = method(**self.kwargs)
        else:
            raise NotImplementedError(f"Reaction type '{self.reaction_type}' is not implemented.")

    # Reaction methods
    # TODO: some code in here can probably go into functions to be shared between methods
    def _perform_pcr(self, keep_primer_annotations=False):
        self._validate_pcr()
        template_molecule = self.inputs["template"]
        fwd_primer_seq = self.inputs["forward_primer"].sequence
        rev_primer_seq = self.inputs["reverse_primer"].sequence

        if isinstance(template_molecule, SingleStrandNucleicAcidSequence):
            template_fwd_seq = template_molecule.sequence
            template_rev_seq = complement_sequence(template_fwd_seq)
        else:
            template_fwd_seq = template_molecule.forward_sequence.sequence
            template_rev_seq = template_molecule.reverse_sequence.sequence

        # Find forward primer last base binding index
        fwd_primer_start, fwd_primer_end  = kmer_trimming_search(template_fwd_seq, fwd_primer_seq)
        # Find reverse primer first base binding index
        rev_primer_start, rev_primer_end = kmer_trimming_search(template_rev_seq, rev_primer_seq, trim_front=False)

        pcr_fwd_seq = fwd_primer_seq + template_fwd_seq[fwd_primer_end + 1:rev_primer_start] + complement_sequence(rev_primer_seq)
        pcr_rev_seq = complement_sequence(fwd_primer_seq) + template_rev_seq[fwd_primer_end + 1:rev_primer_start] + rev_primer_seq

        pcr_fwd_strand = SingleStrandNucleicAcidSequence(
            sequence=pcr_fwd_seq,
            strand_direction=StrandDirections.FWD_STRAND.value,
        )
        pcr_rev_strand = SingleStrandNucleicAcidSequence(
            sequence=pcr_rev_seq,
            strand_direction=StrandDirections.REV_STRAND.value,
        )

        # Collect annotations
        fwd_overhang_len = len(fwd_primer_seq) - len(range(fwd_primer_start, fwd_primer_end + 1))
        new_rev_primer_start = rev_primer_start + fwd_overhang_len
        if isinstance(template_molecule, SingleStrandNucleicAcidSequence):
            pcr_fwd_strand = transfer_annotations(
                template_molecule.annotations,
                pcr_fwd_strand,
                lower_bound_index=fwd_primer_start,
                higher_bound_index=rev_primer_end,
                index_adjustment=fwd_overhang_len,
            )
        else:
            pcr_fwd_strand = transfer_annotations(
                template_molecule.forward_sequence.annotations,
                pcr_fwd_strand,
                lower_bound_index=fwd_primer_start,
                higher_bound_index=rev_primer_end,
                index_adjustment=fwd_overhang_len,
            )
            pcr_rev_strand = transfer_annotations(
                template_molecule.reverse_sequence.annotations,
                pcr_rev_strand,
                lower_bound_index=fwd_primer_start,
                higher_bound_index=rev_primer_end,
                index_adjustment=fwd_overhang_len,
            )
        if keep_primer_annotations:
            pcr_fwd_strand = transfer_annotations(
                self.inputs["forward_primer"].annotations,
                pcr_fwd_strand
            )
            pcr_rev_strand = transfer_annotations(
                self.inputs["reverse_primer"].annotations,
                pcr_rev_strand,
                index_adjustment=new_rev_primer_start,
            )

        pcr_construct = DoubleStrandNucleicAcidSequence(
            forward_sequence=pcr_fwd_strand,
            reverse_sequence=pcr_rev_strand,
        )
        return {"pcr_construct": pcr_construct}

    # TODO: Fill this out
    def _perform_ligation(self):
        # Implement ligation logic here
        # Use self.inputs for ligation inputs
        return {"ligated_sequence": "Simulated ligation output"}

    # TODO: confirm structure of the final smrtbell_library_construct 
    def _perform_pacbio_smrtbell_library_prep(self):
        self._validate_pacbio_smrtbell_library_prep()
        template_fwd_seq = self.inputs["template"].forward_sequence.sequence
        template_rev_seq = self.inputs["template"].reverse_sequence.sequence
        front_adapt_seq = self.inputs["front_adapter"].sequence
        back_adapt_seq = self.inputs["back_adapter"].sequence

        smrtbell_library_construct = SingleStrandNucleicAcidSequence(
            sequence=front_adapt_seq + template_fwd_seq + reverse_sequence(back_adapt_seq) + reverse_sequence(template_rev_seq),
            nucleic_acid_type=NucleicAcidTypes.DNA.value,
            circular=True,
            note="",
        )
        return {"smrtbell_library_construct": smrtbell_library_construct}

    def _validate_pcr(self):
        self._validate_input_keys(["template", "forward_primer", "reverse_primer"])
        if not isinstance(self.inputs["template"], (SingleStrandNucleicAcidSequence, DoubleStrandNucleicAcidSequence)):
            raise ValueError("The template must be a SingleStrandNucleicAcidSequence or "
                             "DoubleStrandNucleicAcidSequence object.")
        if self.inputs["template"].nucleic_acid_type != NucleicAcidTypes.DNA.value:
            raise ValueError("The template cannot be RNA, it must be DNA.")
        if isinstance(self.inputs["template"], SingleStrandNucleicAcidSequence) and self.inputs["template"].strand_direction != StrandDirections.FWD_STRAND.value:
            raise ValueError("Currently, if the template is ssDNA, it must have forward (5->3) strand direction.")
        if not isinstance(self.inputs["forward_primer"], SingleStrandNucleicAcidSequence):
            raise ValueError("The forward primer must be a SingleStrandNucleicAcidSequence object.")
        if self.inputs["forward_primer"].circular:
            raise ValueError("The forward primer cannot be circular, it must be linear.")
        if self.inputs["forward_primer"].nucleic_acid_type != NucleicAcidTypes.DNA.value:
            raise ValueError("The forward primer cannot be RNA, it must be DNA.")
        if not isinstance(self.inputs["reverse_primer"], SingleStrandNucleicAcidSequence):
            raise ValueError("The reverse primer must be a SingleStrandNucleicAcidSequence object.")
        if self.inputs["reverse_primer"].circular:
            raise ValueError("The reverse primer cannot be circular, it must be linear.")
        if self.inputs["reverse_primer"].nucleic_acid_type != NucleicAcidTypes.DNA.value:
            raise ValueError("The reverse primer cannot be RNA, it must be DNA.")

    def _validate_pacbio_smrtbell_library_prep(self):
        self._validate_input_keys(["template", "front_adapter", "back_adapter"])
        if not isinstance(self.inputs["template"], DoubleStrandNucleicAcidSequence):
            raise ValueError("The template must be a DoubleStrandNucleicAcidSequence object.")
        if self.inputs["template"].nucleic_acid_type != NucleicAcidTypes.DNA.value:
            raise ValueError("The template cannot be RNA, it must be DNA.")
        if not isinstance(self.inputs["front_adapter"], SingleStrandNucleicAcidSequence):
            raise ValueError("The front adapter must be a SingleStrandNucleicAcidSequence object.")
        if self.inputs["front_adapter"].circular:
            raise ValueError("The front adapter cannot be circular, it must be linear.")
        if self.inputs["front_adapter"].nucleic_acid_type != NucleicAcidTypes.DNA.value:
            raise ValueError("The front adapter cannot be RNA, it must be DNA.")
        if not isinstance(self.inputs["back_adapter"], SingleStrandNucleicAcidSequence):
            raise ValueError("The back adapter must be a SingleStrandNucleicAcidSequence object.")
        if self.inputs["back_adapter"].circular:
            raise ValueError("The back adapter cannot be circular, it must be linear.")
        if self.inputs["back_adapter"].nucleic_acid_type != NucleicAcidTypes.DNA.value:
            raise ValueError("The back adapter cannot be RNA, it must be DNA.")

    def _validate_input_keys(self, required_keys):
        missing_keys = [key for key in required_keys if key not in self.inputs]
        if missing_keys:
            raise KeyError(f"Missing required keys for a PCR reaction: {', '.join(missing_keys)}.")

    def __repr__(self):
        return f"NucleicAcidReaction(id={self.id}, type={self.reaction_type}, outputs={self.outputs})"

