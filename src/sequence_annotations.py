class SequenceAnnotation:

  def __init__(self, parent_ss_nucleic_acid, name, start, end, note=""):
    self.validate_annotation(parent_ss_nucleic_acid, name, start, end,
                             note)

    self._parent_ss_nucleic_acid_id = parent_ss_nucleic_acid.id
    self._name = name
    self._start = start
    self._end = end
    self._note = note

  def __repr__(self):
    return (f"SequenceAnnotation(name='{self.name}', start={self.start}, "
            f"end={self.end}, note='{self.note}')")

  @property
  def parent_ss_nucleic_acid_id(self):
    return self._parent_ss_nucleic_acid_id

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

  def validate_annotation(self, parent_ss_nucleic_acid, name, start, end,
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
        start >= len(parent_ss_nucleic_acid.sequence) - 1) or (
            end > len(parent_ss_nucleic_acid.sequence) - 1):
      raise ValueError(
          "Cannot add annotation, start and end must be within index bounds of the "
          "parent nucleic acid sequence.")
    if (parent_ss_nucleic_acid.circular is False and start > end):
      raise ValueError(
          "Cannot add annotation, start must be an integer less than end for linear "
          "sequences. Only possible with circular sequences.")
