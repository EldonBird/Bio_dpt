# pcr_lib.pyi
from typing import List

class primer:
    rsid: str
    snp_id: str
    allele: str
    sequence: str
    direction: str
    score: float
    gc: float
    tm: float
    hairpin: float
    homodimer: float
    position: int
    length: int



def reverse_complement(s: str) -> str: ...
def introduce_missmatch(s: str) -> primer: ... #add two more mismatch options 4/10
def generate_allele_specific_primers(data: list[primer], min_length: int, max_length: int) -> list[primer]: ...
def filter_primers(data: list[primer], tm_min: int, tm_max: int, hairpin_max: float, homodimer_max: float) -> list[primer]: ... #redoing 8/10
def rank_primers(data: list[primer], target_tm: float, target_gc: float) -> list[primer]: ... #done pretty sure
def generate_matching_primers(data: list[primer], allele_specific_primers: list[primer], min_distance: int, max_distance: int) -> list[primer]: ... #needs testing 5/10
def check_multiplex_compatibility(data: list[primer], heterodimer_max: float) -> list[primer]: ... # big one 10/10


# introduce mismatch added to 4/10
# filter primers redone 8/10
# generate matching primers tested 3/10
# check multiplex compatibility tested 4/10
# multiplexing 10/10
# look over code to see if it's right - Sophia
# unit tests 6/10 (thorough understanding)
# research how to do what we're doing?


