from src.nucleic_acid_reactions import NucleicAcidReaction
from src.single_stranded_nucleic_acids import SingleStrandNucleicAcidSequence
from src.double_stranded_nucleic_acids import DoubleStrandNucleicAcidSequence
from src.util_classes import NucleicAcidTypes, StrandDirections


def test_pcr_reaction_1():
    template_seq = DoubleStrandNucleicAcidSequence(
        forward_sequence="ATGCGTAATAAGC",
    )
    fwd_primer = SingleStrandNucleicAcidSequence(
        sequence="ATGC",
    )
    rev_primer = SingleStrandNucleicAcidSequence(
        sequence="TTCG",
    )

    pcr_reaction = NucleicAcidReaction(reaction_type="pcr", inputs={
        "template": template_seq,
        "forward_primer": fwd_primer,
        "reverse_primer": rev_primer
    })

    assert pcr_reaction.outputs["pcr_construct"].forward_sequence.sequence == "ATGCGTAATAAGC"
    assert pcr_reaction.outputs["pcr_construct"].reverse_sequence.sequence == "TACGCATTATTCG"

def test_pcr_reaction_2():
    template_seq = DoubleStrandNucleicAcidSequence(
        forward_sequence="ATGCGTAATAAGC",
    )
    fwd_primer = SingleStrandNucleicAcidSequence(
        sequence="GGGGGGGATGC",
    )
    rev_primer = SingleStrandNucleicAcidSequence(
        sequence="TTCGTTT",
    )

    pcr_reaction = NucleicAcidReaction(reaction_type="pcr", inputs={
        "template": template_seq,
        "forward_primer": fwd_primer,
        "reverse_primer": rev_primer
    })

    assert pcr_reaction.outputs["pcr_construct"].forward_sequence.sequence == "GGGGGGGATGCGTAATAAGCAAA"
    assert pcr_reaction.outputs["pcr_construct"].reverse_sequence.sequence == "CCCCCCCTACGCATTATTCGTTT"


def test_pcr_reaction_3():
    template_seq = DoubleStrandNucleicAcidSequence(
        forward_sequence="TAGTAGCCCCCGGGGGAAAAATTTTTAAAAAAAAAA",
    )
    fwd_primer = SingleStrandNucleicAcidSequence(
        sequence="CGCGCGTAGTAG",
    )
    rev_primer = SingleStrandNucleicAcidSequence(
        sequence="TTTTTTTTTTTTTTTTTTTT",
    )

    pcr_reaction = NucleicAcidReaction(reaction_type="pcr", inputs={
        "template": template_seq,
        "forward_primer": fwd_primer,
        "reverse_primer": rev_primer
    })

    # print("template", template_seq)
    # print("construct", pcr_reaction.outputs["pcr_construct"])

    assert pcr_reaction.outputs["pcr_construct"].forward_sequence.sequence == "CGCGCGTAGTAGCCCCCGGGGGAAAAATTTTTAAAAAAAAAAAAAAAAAAAA"
    assert pcr_reaction.outputs["pcr_construct"].reverse_sequence.sequence == "GCGCGCATCATCGGGGGCCCCCTTTTTAAAAATTTTTTTTTTTTTTTTTTTT"


def test_pcr_reaction_annotations():
    template_seq = DoubleStrandNucleicAcidSequence(
        forward_sequence="ATGCGTAATAAGC",
    )
    template_seq.forward_sequence.add_annotations([
        {"name": "Annot1", "start": 0, "end": 5, "note": "This is a test"}])
    fwd_primer = SingleStrandNucleicAcidSequence(
        sequence="ATGC",
    )
    rev_primer = SingleStrandNucleicAcidSequence(
        sequence="TTCG",
        annotations=[
            {"name": "Annot2", "start": 0, "end": 3, "note": "This is the rev binding portion"}]
    )

    pcr_reaction = NucleicAcidReaction(
        reaction_type="pcr",
        keep_primer_annotations=True,
        inputs={
        "template": template_seq,
        "forward_primer": fwd_primer,
        "reverse_primer": rev_primer
        })

    print(pcr_reaction)

# TODO: Add more test functions for different scenarios
