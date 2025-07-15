import pandas as pd
import primer3
from Bio import Entrez, Seq
import requests
import json
from typing import List, Dict, Tuple, Optional

# Set your email for NCBI Entrez
Entrez.email = "your.email@example.com"  # Replace with your actual email

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return str(Seq.Seq(sequence).reverse_complement())

def fetch_snp_data(rsids: List[str], flank_length: int = 800) -> pd.DataFrame:
    """
    Fetch SNP data from dbSNP (mocked for demo).
    TODO: Replace with real Entrez/API call to fetch rsIDs, alleles, and flanking sequences from dbSNP.
    - Use Bio.Entrez.esearch and efetch to query SNP database.
    - Parse XML/JSON output for alleles and flanking sequences.
    - Handle errors for invalid rsIDs or missing data.
    """
    try:
        snp_data = []
        for rsid in rsids:
            # Mock data: assumes two alleles and 800 bp flanks
            alleles = ["A", "G"]
            flank_seq = "A" * flank_length + "[" + "/".join(alleles) + "]" + "T" * flank_length
            for allele in alleles:
                snp_data.append({
                    "snpID": rsid,
                    "allele": allele,
                    "sequence": flank_seq.replace(f'[]', allele),
                    "position": flank_length  # SNP position
                })
        return pd.DataFrame(snp_data)
    except Exception as e:
        print(f"Error fetching SNP data: {e}")
        return pd.DataFrame()

def introduce_mismatch(sequence: str, position: int) -> List[str]:
    """
    Introduce a mismatch at the antepenultimate position (3rd from last).
    TODO: Validate mismatch rule for biological accuracy.
    - Consider additional mismatch types if needed for specificity.
    - Add checks for invalid sequences (e.g., non-ACGT bases).
    """
    mismatch_rules = {"A": "G", "G": "A", "C": "T", "T": "C"}
    primers = []
    if position < 3 or position > len(sequence):
        return [sequence]
    target_base = sequence[position - 3]
    new_base = mismatch_rules.get(target_base, "A")
    primer = sequence[:position - 3] + new_base + sequence[position - 2:]
    primers.append(primer)
    return primers

def generate_allele_specific_primers(snp_data: pd.DataFrame, min_len: int = 18, max_len: int = 28) -> pd.DataFrame:
    """
    Generate allele-specific primers (forward/reverse) ending at the SNP.
    TODO: Optimize for large SNP sets.
    - Use parallel processing (e.g., multiprocessing) for many SNPs.
    - Add validation for sequence length and SNP position.
    """
    primers = []
    for _, row in snp_data.iterrows():
        snp_id = row["snpID"]
        allele = row["allele"]
        sequence = row["sequence"]
        center = row["position"]
        for length in range(min_len, max_len + 1):
            # Forward primer: upstream sequence ending at SNP
            forward = sequence[center - length + 1:center + 1]
            forward_mismatches = introduce_mismatch(forward, length)
            for fm in forward_mismatches:
                primers.append({
                    "snpID": snp_id,
                    "allele": allele,
                    "primer_sequence": fm,
                    "direction": "forward",
                    "length": length
                })
            # Reverse primer: downstream sequence, reverse complemented
            reverse = reverse_complement(sequence[center:center + length])
            reverse_mismatches = introduce_mismatch(reverse, length)
            for rm in reverse_mismatches:
                primers.append({
                    "snpID": snp_id,
                    "allele": allele,
                    "primer_sequence": rm,
                    "direction": "reverse",
                    "length": length
                })
    return pd.DataFrame(primers)

def evaluate_primer(primer_seq: str) -> Dict:
    """
    Evaluate primer quality using primer3-py.
    TODO: Enhance error handling.
    - Handle primer3-py failures gracefully.
    - Add logging for failed evaluations.
    """
    try:
        result = primer3.bindings.designPrimers(
            {
                "SEQUENCE_TEMPLATE": primer_seq,
                "SEQUENCE_PRIMER": primer_seq
            },
            {
                "PRIMER_OPT_SIZE": 20,
                "PRIMER_MIN_SIZE": 18,
                "PRIMER_MAX_SIZE": 28,
                "PRIMER_OPT_TM": 62.5,
                "PRIMER_MIN_TM": 60.0,
                "PRIMER_MAX_TM": 65.0,
                "PRIMER_MIN_GC": 40.0,
                "PRIMER_MAX_GC": 60.0,
                "PRIMER_MAX_HAIRPIN_TH": 45.0,
                "PRIMER_MAX_SELF_ANY_TH": 45.0
            }
        )
        return {
            "tm": result.get("PRIMER_LEFT_0_TM", 0),
            "gc_content": result.get("PRIMER_LEFT_0_GC_PERCENT", 0),
            "hairpin": result.get("PRIMER_LEFT_0_HAIRPIN_TH", 0),
            "homodimer": result.get("PRIMER_LEFT_0_SELF_ANY_TH", 0)
        }
    except Exception:
        return {"tm": 0, "gc_content": 0, "hairpin": 999, "homodimer": 999}

def filter_primers(primers: pd.DataFrame, tm_min: float = 60.0, tm_max: float = 65.0, hairpin_max: float = 45.0, homodimer_max: float = 45.0) -> pd.DataFrame:
    """
    Filter primers based on quality metrics.
    TODO: Add fallback for strict filtering.
    - If no primers pass, consider relaxing criteria or selecting best available.
    """
    quality_metrics = []
    for _, row in primers.iterrows():
        metrics = evaluate_primer(row["primer_sequence"])
        quality_metrics.append({
            "snpID": row["snpID"],
            "allele": row["allele"],
            "primer_sequence": row["primer_sequence"],
            "direction": row["direction"],
            "length": row["length"],
            "tm": metrics["tm"],
            "gc_content": metrics["gc_content"],
            "hairpin": metrics["hairpin"],
            "homodimer": metrics["homodimer"]
        })
    df = pd.DataFrame(quality_metrics)
    return df[
        (df["tm"].between(tm_min, tm_max)) &
        (df["gc_content"].between(40, 60)) &
        (df["hairpin"] < hairpin_max) &
        (df["homodimer"] < homodimer_max)
        ]

def rank_primers(primers: pd.DataFrame) -> pd.DataFrame:
    """
    Rank primers based on Tm proximity to 62.5Â°C and GC content.
    TODO: Refine ranking criteria.
    - Consider weighting Tm vs. GC scores.
    - Add user-configurable ranking metrics.
    """
    primers["tm_score"] = abs(primers["tm"] - 62.5)
    primers["gc_score"] = abs(primers["gc_content"] - 50.0)
    primers["score"] = primers["tm_score"] + primers["gc_score"] + primers["hairpin"] + primers["homodimer"]
    return primers.sort_values("score").groupby(["snpID", "allele", "direction"]).head(5)

def generate_matching_primers(snp_data: pd.DataFrame, allele_specific_primers: pd.DataFrame, min_dist: int = 100, max_dist: int = 500) -> pd.DataFrame:
    """
    Generate matching primers for top 5 allele-specific primers.
    TODO: Optimize primer pairing.
    - Use primer3-py's designPrimers for more efficient pairing.
    - Add checks for primer pair compatibility (e.g., Tm difference < 5Â°C).
    """
    matching_primers = []
    for _, asp in allele_specific_primers.iterrows():
        snp_id = asp["snpID"]
        allele = asp["allele"]
        direction = asp["direction"]
        sequence = snp_data[snp_data["snpID"] == snp_id]["sequence"].iloc[0]
        center = snp_data[snp_data["snpID"] == snp_id]["position"].iloc[0]
        for dist in range(min_dist, max_dist + 1, 10):
            for length in range(18, 29):
                if direction == "forward":
                    # Matching reverse primer: downstream, reverse complemented
                    start = center + dist
                    primer_seq = reverse_complement(sequence[start:start + length])
                else:
                    # Matching forward primer: upstream
                    start = center - dist - length
                    primer_seq = sequence[start:start + length]
                metrics = evaluate_primer(primer_seq)
                if (60.0 <= metrics["tm"] <= 65.0 and
                        40.0 <= metrics["gc_content"] <= 60.0 and
                        metrics["hairpin"] < 45.0 and
                        metrics["homodimer"] < 45.0):
                    matching_primers.append({
                        "snpID": snp_id,
                        "allele": allele,
                        "asp_sequence": asp["primer_sequence"],
                        "matching_sequence": primer_seq,
                        "matching_direction": "reverse" if direction == "forward" else "forward",
                        "tm": metrics["tm"],
                        "gc_content": metrics["gc_content"],
                        "hairpin": metrics["hairpin"],
                        "homodimer": metrics["homodimer"],
                        "distance": dist
                    })
    return pd.DataFrame(matching_primers).groupby(["snpID", "allele", "asp_sequence"]).head(1)

def check_multiplex_compatibility(primer_pairs: pd.DataFrame, heterodimer_max: float = 50.0) -> pd.DataFrame:
    """
    Check primer sets for multiplex compatibility.
    TODO: Enhance for multiple SNPs.
    - Extend to check cross-SNP interactions (current checks within SNP).
    - Optimize for large primer sets using batch dimer calculations.
    """
    compatible_pairs = []
    for snp_id in primer_pairs["snpID"].unique():
        snp_pairs = primer_pairs[primer_pairs["snpID"] == snp_id]
        for i, pair1 in snp_pairs.iterrows():
            is_compatible = True
            for j, pair2 in snp_pairs.iterrows():
                if i != j:
                    result = primer3.bindings.calcHeterodimer(pair1["asp_sequence"], pair2["asp_sequence"])
                    if result.tm > heterodimer_max:
                        is_compatible = False
                        break
                    result = primer3.bindings.calcHeterodimer(pair1["matching_sequence"], pair2["matching_sequence"])
                    if result.tm > heterodimer_max:
                        is_compatible = False
                        break
            if is_compatible:
                compatible_pairs.append(pair1)
    return pd.DataFrame(compatible_pairs)

def export_results(primer_pairs: pd.DataFrame, output_format: str = "csv", output_file: str = "primer_sets") -> None:
    """
    Export primer sets to CSV or JSON.
    TODO: Enhance output flexibility.
    - Add support for custom metadata fields.
    - Implement JSON output validation.
    """
    export_data = primer_pairs[["snpID", "allele", "asp_sequence", "matching_sequence", "direction", "matching_direction", "tm", "gc_content", "hairpin", "homodimer", "distance"]]
    if output_format == "csv":
        export_data.to_csv(f"{output_file}.csv", index=False)
    elif output_format == "json":
        export_data.to_json(f"{output_file}.json", orient="records", indent=2)
    # TODO: Add logging for successful export or errors

def main(rsids: List[str], output_format: str = "csv"):
    """
    Main function to generate, filter, rank, pair, and export primers.
    TODO: Add comprehensive error handling and logging.
    - Log progress and errors to a file.
    - Add input validation for rsIDs and output format.
    TODO: Test with real SNP data.
    - Validate output with biological experts.
    - Benchmark performance for large SNP sets.
    """
    # Fetch SNP data
    snp_data = fetch_snp_data(rsids)
    if snp_data.empty:
        print("No SNP data retrieved.")
        return

    # Generate allele-specific primers
    allele_specific_primers = generate_allele_specific_primers(snp_data)
    if allele_specific_primers.empty:
        print("No allele-specific primers generated.")
        return

    # Filter primers
    filtered_primers = filter_primers(allele_specific_primers)
    if filtered_primers.empty:
        print("No primers passed quality filters.")
        return

    # Rank primers
    ranked_primers = rank_primers(filtered_primers)
    if ranked_primers.empty:
        print("No primers available after ranking.")
        return

    # Generate matching primers
    primer_pairs = generate_matching_primers(snp_data, ranked_primers)
    if primer_pairs.empty:
        print("No matching primers generated.")
        return

    # Check multiplex compatibility
    compatible_pairs = check_multiplex_compatibility(primer_pairs)
    if compatible_pairs.empty:
        print("No primer sets are multiplex compatible.")
        return

    # Export results
    export_results(compatible_pairs, output_format)
    print(f"Results exported to primer_sets.{output_format}")

if __name__ == "__main__":
    rsids = ["rs9462492", "rs58318008", "rs1421085", "rs9939609", "rs1121980"]
    main(rsids)