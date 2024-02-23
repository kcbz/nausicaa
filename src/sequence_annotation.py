class SequenceAnnotation:

  def __repr__(self):
    return (f"SequenceAnnotation(name='{self.name}', start={self.start}, "
            f"end={self.end}, note='{self.note}')")

  def __init__(self, parent_nucleic_acid_sequence, name, start, end, note=""):
    self.validate_annotation(parent_nucleic_acid_sequence, name, start, end,
                             note)

    self._parent_nucleic_acid_sequence_id = parent_nucleic_acid_sequence.id
    self._name = name
    self._start = start
    self._end = end
    self._note = note

  @property
  def parent_nucleic_acid_sequence_id(self):
    return self._parent_nucleic_acid_sequence_id

  @property
  def name(self):
    return self._name

  @property
  def start(self):
    return self._start

  @property
  def end(self):
    return self._end

  @property
  def note(self):
    return self._note

  def validate_annotation(self, parent_nucleic_acid_sequence, name, start, end,
                          note):
    if not isinstance(name, str):
      raise TypeError("Cannot add annotation, name must be a string.")
    if not isinstance(start, int):
      raise TypeError("Cannot add annotation, start must be an integer.")
    if not isinstance(end, int):
      raise TypeError("Cannot add annotation, end must be an integer.")
    if not isinstance(note, str):
      raise TypeError("Cannot add annotation, note must be a string.")
    if (start == end) or (start < 0) or (end <= 0) or (
        start >= len(parent_nucleic_acid_sequence.sequence) - 1) or (
            end > len(parent_nucleic_acid_sequence.sequence) - 1):
      raise ValueError(
          "Cannot add annotation, start and end must be within index bounds of the "
          "parent nucleic acid sequence.")
    if (parent_nucleic_acid_sequence.circular is False and start > end):
      raise ValueError(
          "Cannot add annotation, start must be an integer less than end for linear "
          "sequences. Only possible with circular sequences.")
