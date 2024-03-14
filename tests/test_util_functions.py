from src.util_functions import kmer_trimming_search


def test_kmer_trimming_search():
    seq = "ATGCTTTGGCGATCGA"
    fwd_query_seq = "ATGCTT"
    expected_fwd_search_index_start = 0
    expected_fwd_search_index_end = 5
    rev_query_seq = "CGATC"
    expected_rev_search_index_start = 9
    expected_rev_search_index_end = 13

    fwd_search_index_start, fwd_search_index_end = kmer_trimming_search(seq, fwd_query_seq)
    rev_search_index_start, rev_search_index_end = kmer_trimming_search(seq, rev_query_seq, trim_front=False)
    assert fwd_search_index_start == expected_fwd_search_index_start
    assert fwd_search_index_end == expected_fwd_search_index_end
    assert rev_search_index_start == expected_rev_search_index_start
    assert rev_search_index_end == expected_rev_search_index_end


def test_kmer_trimming_search_overhangs():
    seq = "ATGCTTTGGCGATCGA"
    fwd_query_seq = "GGTGCCAAGTGCAGTATGCTT"
    expected_fwd_search_index_start = 0
    expected_fwd_search_index_end = 5
    rev_query_seq = "CGATCAAAAAAAAAAAAAAACG"
    expected_rev_search_index_start = 9
    expected_rev_search_index_end = 13

    fwd_search_index_start, fwd_search_index_end = kmer_trimming_search(seq, fwd_query_seq)
    rev_search_index_start, rev_search_index_end = kmer_trimming_search(seq, rev_query_seq, trim_front=False)
    assert fwd_search_index_start == expected_fwd_search_index_start
    assert fwd_search_index_end == expected_fwd_search_index_end
    assert rev_search_index_start == expected_rev_search_index_start
    assert rev_search_index_end == expected_rev_search_index_end
