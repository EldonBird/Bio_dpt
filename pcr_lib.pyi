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
def introduce_missmatch(s: str) -> primer: ...
def generate_allele_specific_primers(data: list[primer], min_length: int, max_length: int) -> list[primer]: ...
def filter_primers(data: list[primer], tm_min: int, tm_max: int, hairpin_max: float, homodimer_max: float) -> list[primer]: ...
def rank_primers(data: list[primer], target_tm: float, target_gc: float) -> list[primer]: ...
def generate_matching_primers(data: list[primer], allele_specific_primers: list[primer], min_distance: int, max_distance: int) -> list[primer]: ...
def check_multiplex_compatibility(data: list[primer], heterodimer_max: float) -> list[primer]: ...
