from typing import Union


def kmer_trimming_search(template_seq: str, query_seq: str, trim_front=True) -> Union[int, int]:
    """
    Searches for a primer binding site on a template sequence by progressively trimming
    the primer sequence from the start or end until a perfect match is found.

    Parameters:
    - template_seq (str): The DNA sequence of the template.
    - primer_seq (str): The DNA sequence of the primer.
    - trim_direction (str): Direction to trim the primer sequence. 'start' or 'end'.

    Returns:
    - tuple (int, int): THe index positions of the first and last base that the query sequence matches
        to the template sequence.Returns -1 if no binding site is found.
    """

    template_seq = template_seq.upper()
    query_seq = query_seq.upper()
    if trim_front:
        for i in range(len(query_seq)):
            trimmed_query_seq = query_seq[i:]
            if trimmed_query_seq in template_seq:
                start = template_seq.index(trimmed_query_seq)
                end = template_seq.index(trimmed_query_seq) + len(trimmed_query_seq) - 1
                return (start, end)
    # return first base index of match 
    else:
        for i in range(len(query_seq), 0, -1):
            trimmed_query_seq = query_seq[:i]
            if trimmed_query_seq in template_seq:
                start = template_seq.index(trimmed_query_seq)
                end = template_seq.index(trimmed_query_seq) + len(trimmed_query_seq) - 1
                return (start, end)
    return -1
