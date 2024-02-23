from src.nucleic_acids import (
    reverse_sequence,
    complement_sequence,
    reverse_complement_sequence,
)


def test_reverse_sequence():
    seq = "ATGC"
    expected_rev_seq = "CGTA"
    rev_seq = reverse_sequence(seq)
    assert rev_seq == expected_rev_seq

def test_compliment_sequence():
    seq = "ATGC"
    expected_comp_seq = "TACG"
    comp_seq = complement_sequence(seq)
    assert comp_seq == expected_comp_seq

def test_compliment_sequence_RNA():
    seq = "ATGC"
    expected_comp_seq = "UACG"
    comp_seq = complement_sequence(seq, nuc_type="RNA")
    assert comp_seq == expected_comp_seq

def test_reverse_complement_sequence():
    seq = "ATGC"
    expected_revcomp_seq = "GCAT"
    revcomp_seq = reverse_complement_sequence(seq)
    assert revcomp_seq == expected_revcomp_seq
