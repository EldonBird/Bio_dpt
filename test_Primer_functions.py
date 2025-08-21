import pytest
import pandas as pd
import sys
import os
from unittest.mock import patch, MagicMock
import logging

# Add the directory containing Primer_functions.py to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import the functions from Primer_functions.py
from Primer_functions import (
    Introduce_Mismatch,
    Calc_GC_Content,
    Evaluate_Primers,
    rank_primers,
    Filter_Primers,
    Generate_Allele_Specific_Primers,
    Make_Primers
)

class TestIntroduceMismatch:
    """Test suite for Introduce_Mismatch function"""
    
    def test_introduce_mismatch_valid_sequence(self):
        """Test mismatch introduction with valid sequence"""
        result = Introduce_Mismatch("ATCGATCG")
        # 3rd from last is 'T', should become 'C'
        assert result == "ATCGACCG"
    
    def test_introduce_mismatch_purine_to_purine(self):
        """Test A->G and G->A transformations"""
        assert Introduce_Mismatch("ATCGA") == "ATCGG"  # A->G
        assert Introduce_Mismatch("ATCGG") == "ATCGA"  # G->A
    
    def test_introduce_mismatch_pyrimidine_to_pyrimidine(self):
        """Test C->T and T->C transformations"""
        assert Introduce_Mismatch("ATCGC") == "ATCGT"  # C->T
        assert Introduce_Mismatch("ATCGT") == "ATCGC"  # T->C
    
    def test_introduce_mismatch_short_sequence(self, capsys):
        """Test with sequence too short for mismatch"""
        result = Introduce_Mismatch("AT")
        captured = capsys.readouterr()
        assert "Warning: Primer too short for mismatch" in captured.out
        assert result == "AT"
    
    def test_introduce_mismatch_empty_sequence(self, capsys):
        """Test with empty sequence"""
        result = Introduce_Mismatch("")
        captured = capsys.readouterr()
        assert "Warning: Invalid primer input" in captured.out
        assert result == ""
    
    def test_introduce_mismatch_invalid_characters(self, capsys):
        """Test with invalid DNA characters"""
        result = Introduce_Mismatch("ATCGN")
        captured = capsys.readouterr()
        assert "Warning: Invalid characters in primer" in captured.out
        assert result == "ATCGN"
    
    def test_introduce_mismatch_lowercase_input(self):
        """Test with lowercase input"""
        result = Introduce_Mismatch("atcgatcg")
        assert result == "ATCGACCG"  # Should convert to uppercase and process
    
    def test_introduce_mismatch_none_input(self, capsys):
        """Test with None input"""
        result = Introduce_Mismatch(None)
        captured = capsys.readouterr()
        assert "Warning: Invalid primer input" in captured.out
        assert result is None


class TestCalcGCContent:
    """Test suite for Calc_GC_Content function"""
    
    def test_calc_gc_content_all_gc(self):
        """Test with sequence containing only G and C"""
        result = Calc_GC_Content("GCGCGC")
        assert result == 1.0
    
    def test_calc_gc_content_no_gc(self):
        """Test with sequence containing no G or C"""
        result = Calc_GC_Content("ATATAT")
        assert result == 0.0
    
    def test_calc_gc_content_half_gc(self):
        """Test with sequence containing 50% GC"""
        result = Calc_GC_Content("ATGC")
        assert result == 0.5
    
    def test_calc_gc_content_mixed(self):
        """Test with mixed sequence"""
        result = Calc_GC_Content("ATCGATCG")
        # 4 GC out of 8 total = 0.5
        assert result == 0.5
    
    def test_calc_gc_content_single_nucleotide(self):
        """Test with single nucleotides"""
        assert Calc_GC_Content("G") == 1.0
        assert Calc_GC_Content("A") == 0.0


class TestEvaluatePrimers:
    """Test suite for Evaluate_Primers function"""
    
    @pytest.fixture
    def sample_primer_dict(self):
        return {
            "snpID": "rs123",
            "allele": "A",
            "primer_sequence": "ATCGATCGATCGATCG",
            "direction": "forward",
            "length": 16
        }
    
    @patch('primer3.bindings.calc_tm')
    @patch('primer3.bindings.calc_hairpin')
    @patch('primer3.bindings.calc_homodimer')
    def test_evaluate_primers_success(self, mock_homodimer, mock_hairpin, mock_tm, sample_primer_dict):
        """Test successful primer evaluation"""
        # Mock primer3 results
        mock_tm.return_value = 60.5
        mock_hairpin.return_value = MagicMock(tm=2.0)
        mock_homodimer.return_value = MagicMock(tm=1.5)
        
        result = Evaluate_Primers(sample_primer_dict)
        
        assert result["tm"] == 60.5
        assert result["gc_content"] == 0.5  # 8 GC out of 16 total
        assert result["hairpin_dg"] == 2.0
        assert result["homodimer_dg"] == 1.5
        assert result["success"] is True
    
    @patch('primer3.bindings.calc_tm')
    def test_evaluate_primers_exception(self, mock_tm, sample_primer_dict):
        """Test primer evaluation with exception"""
        mock_tm.side_effect = Exception("Primer3 error")
        
        result = Evaluate_Primers(sample_primer_dict)
        
        assert result["tm"] == 0
        assert result["gc_content"] == 0
        assert result["hairpin_dg"] == 999
        assert result["homodimer_dg"] == 999
        assert result["success"] is False
    
    @patch('primer3.bindings.calc_tm')
    @patch('primer3.bindings.calc_hairpin')
    @patch('primer3.bindings.calc_homodimer')
    def test_evaluate_primers_no_tm_attribute(self, mock_homodimer, mock_hairpin, mock_tm, sample_primer_dict):
        """Test primer evaluation when hairpin/homodimer objects lack tm attribute"""
        mock_tm.return_value = 60.5
        mock_hairpin.return_value = MagicMock()
        del mock_hairpin.return_value.tm  # Remove tm attribute
        mock_homodimer.return_value = MagicMock()
        del mock_homodimer.return_value.tm  # Remove tm attribute
        
        result = Evaluate_Primers(sample_primer_dict)
        
        assert result["hairpin_dg"] == 0
        assert result["homodimer_dg"] == 0
        assert result["success"] is True


class TestRankPrimers:
    """Test suite for rank_primers function"""
    
    @pytest.fixture
    def sample_primers(self):
        return [
            {"primer_sequence": "PRIMER1", "tm": 60.0, "gc_content": 50.0, "hairpin_dg": -5.0, "homodimer_dg": -3.0},
            {"primer_sequence": "PRIMER2", "tm": 65.0, "gc_content": 45.0, "hairpin_dg": -8.0, "homodimer_dg": -6.0},
            {"primer_sequence": "PRIMER3", "tm": 62.5, "gc_content": 50.0, "hairpin_dg": -2.0, "homodimer_dg": -1.0},
        ]
    
    def test_rank_primers_default_targets(self, sample_primers):
        """Test primer ranking with default targets"""
        result = rank_primers(sample_primers.copy())
        
        # PRIMER3 should rank highest (closest to targets)
        assert result[0]["primer_sequence"] == "PRIMER3"
        assert len(result) <= 5  # Default optimism value
    
    def test_rank_primers_custom_targets(self, sample_primers):
        """Test primer ranking with custom targets"""
        result = rank_primers(sample_primers.copy(), target_tm=60.0, target_gc=45.0, optimism=2)
        
        assert len(result) <= 2
        # Verify scores are calculated
        for primer in result:
            assert "tm_score" in primer
            assert "gc_score" in primer
            assert "score" in primer
    
    def test_rank_primers_empty_list(self):
        """Test ranking with empty primer list"""
        result = rank_primers([])
        assert result == []
    
    def test_rank_primers_single_primer(self):
        """Test ranking with single primer"""
        primers = [{"primer_sequence": "PRIMER1", "tm": 60.0, "gc_content": 50.0, "hairpin_dg": -5.0, "homodimer_dg": -3.0}]
        result = rank_primers(primers)
        assert len(result) == 1


class TestFilterPrimers:
    """Test suite for Filter_Primers function"""
    
    @pytest.fixture
    def sample_primers_df(self):
        return pd.DataFrame([
            {"primer_sequence": "PRIMER1", "tm": 62.0, "gc_content": 50.0, "hairpin_dg": -5.0, "homodimer_dg": -3.0},
            {"primer_sequence": "PRIMER2", "tm": 70.0, "gc_content": 30.0, "hairpin_dg": -15.0, "homodimer_dg": -20.0},
            {"primer_sequence": "PRIMER3", "tm": 63.5, "gc_content": 55.0, "hairpin_dg": -2.0, "homodimer_dg": -1.0},
        ])
    
    def test_filter_primers_strict_pass(self, sample_primers_df):
        """Test filtering with primers that pass strict criteria"""
        result = Filter_Primers(sample_primers_df)
        
        # PRIMER1 and PRIMER3 should pass, PRIMER2 should fail
        assert len(result) == 2
        assert "filter_level" in result.columns
        assert (result["filter_level"] == "strict").all()
    
    def test_filter_primers_empty_dataframe(self):
        """Test filtering with empty DataFrame"""
        empty_df = pd.DataFrame()
        result = Filter_Primers(empty_df)
        assert result.empty
    
    def test_filter_primers_no_strict_pass_with_fallback(self, capsys):
        """Test filtering when no primers pass strict but some pass relaxed"""
        # Create primers that fail strict but might pass relaxed
        bad_primers_df = pd.DataFrame([
            {"primer_sequence": "PRIMER1", "tm": 58.0, "gc_content": 35.0, "hairpin_dg": -12.0, "homodimer_dg": -15.0},
        ])
        
        result = Filter_Primers(bad_primers_df, use_fallback=True)
        captured = capsys.readouterr()
        
        assert "WARNING: No primers passed strict filtering" in captured.out
        if not result.empty:
            assert "filter_level" in result.columns
            assert (result["filter_level"] == "relaxed").all()
    
    def test_filter_primers_no_fallback(self, capsys):
        """Test filtering with fallback disabled"""
        bad_primers_df = pd.DataFrame([
            {"primer_sequence": "PRIMER1", "tm": 50.0, "gc_content": 20.0, "hairpin_dg": -20.0, "homodimer_dg": -25.0},
        ])
        
        result = Filter_Primers(bad_primers_df, use_fallback=False)
        assert result.empty
    
    @patch('Primer_functions.Evaluate_Primers')
    def test_filter_primers_missing_metrics(self, mock_evaluate, sample_primers_df):
        """Test filtering when metrics need to be calculated"""
        # Remove metric columns to trigger evaluation
        df_no_metrics = sample_primers_df[["primer_sequence"]].copy()
        mock_evaluate.return_value = {"tm": 62.0, "gc_content": 50.0, "hairpin_dg": -5.0, "homodimer_dg": -3.0}
        
        result = Filter_Primers(df_no_metrics)
        
        # Should have called Evaluate_Primers for each row
        assert mock_evaluate.call_count == len(df_no_metrics)


class TestMakePrimers:
    """Test suite for Make_Primers function"""
    
    @patch('Primer_functions.Evaluate_Primers')
    def test_make_primers_valid_sequence(self, mock_evaluate):
        """Test making primers from valid sequence"""
        mock_evaluate.return_value = {
            "snpID": "rs123", "allele": "A", "primer_sequence": "ATCGATCG",
            "direction": "forward", "length": 8, "tm": 60.0, "success": True
        }
        
        result = Make_Primers("ATCGATCGATCG", 6, 10, "rs123", "A")
        
        # Should create primers of different lengths
        assert len(result) > 0
        assert all(primer["snpID"] == "rs123" for primer in result)
        assert all(primer["allele"] == "A" for primer in result)
        assert all(primer["direction"] == "forward" for primer in result)
    
    @patch('Primer_functions.Evaluate_Primers')
    def test_make_primers_reverse_direction(self, mock_evaluate):
        """Test making primers in reverse direction"""
        mock_evaluate.return_value = {
            "snpID": "rs123", "allele": "A", "primer_sequence": "ATCGATCG",
            "direction": "reverse", "length": 8, "tm": 60.0, "success": True
        }
        
        result = Make_Primers("ATCGATCGATCG", 6, 10, "rs123", "A", "reverse")
        
        assert all(primer["direction"] == "reverse" for primer in result)
    
    def test_make_primers_sequence_too_short(self, capsys):
        """Test with sequence shorter than minimum length"""
        result = Make_Primers("ATCG", 10, 15, "rs123", "A")
        captured = capsys.readouterr()
        
        assert "The length of your forward primer wasn't long enough" in captured.out
        assert result == []


class TestGenerateAlleleSpecificPrimers:
    """Test suite for Generate_Allele_Specific_Primers function"""
    
    @pytest.fixture
    def sample_snp_data(self):
        return [
            {
                'snpID': 'rs1799971',
                'allele': 'A',
                'sequence': 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG',
                'position': 25
            }
        ]
    
    @patch('Primer_functions.Make_Primers')
    @patch('Primer_functions.Introduce_Mismatch')
    def test_generate_allele_specific_primers(self, mock_mismatch, mock_make_primers, sample_snp_data):
        """Test generation of allele-specific primers"""
        mock_mismatch.return_value = "MISMATCH_PRIMER"
        mock_make_primers.return_value = [{"primer": "test"}]
        
        result = Generate_Allele_Specific_Primers(sample_snp_data)
        
        # Should return list of lists (one per SNP, containing forward/reverse lists)
        assert len(result) == 1  # One SNP
        assert len(result[0]) == 2  # Forward and reverse directions
        
        # Should have called Make_Primers twice per SNP (forward and reverse)
        assert mock_make_primers.call_count == 2
    
    def test_generate_allele_specific_primers_empty_list(self):
        """Test with empty SNP list"""
        result = Generate_Allele_Specific_Primers([])
        assert result == []
    
    @patch('Primer_functions.Make_Primers')
    @patch('Primer_functions.Introduce_Mismatch')
    def test_generate_allele_specific_primers_custom_lengths(self, mock_mismatch, mock_make_primers, sample_snp_data):
        """Test with custom min/max lengths"""
        mock_mismatch.return_value = "MISMATCH_PRIMER"
        mock_make_primers.return_value = [{"primer": "test"}]
        
        result = Generate_Allele_Specific_Primers(sample_snp_data, min_len=20, max_len=30)
        
        # Check that Make_Primers was called with adjusted lengths
        # (min_len-2, max_len-1 according to the function logic)
        calls = mock_make_primers.call_args_list
        for call in calls:
            args, kwargs = call
            assert args[1] == 18  # min_len - 2
            assert args[2] == 29  # max_len - 1


class TestDataValidation:
    """Test data validation for SNP input"""
    
    def test_validate_snp_structure(self):
        """Test that SNP data has required fields"""
        snp_data = [
            {'snpID': 'rs1799971', 'allele': 'A', 'sequence': 'ATCG', 'position': 2}
        ]
        
        for snp in snp_data:
            assert 'snpID' in snp
            assert 'allele' in snp
            assert 'sequence' in snp
            assert 'position' in snp
    
    def test_sequence_is_valid_dna(self):
        """Test that sequences contain only valid DNA bases"""
        valid_sequence = "ATCGATCGATCG"
        valid_bases = set('ATCG')
        sequence_bases = set(valid_sequence.upper())
        assert sequence_bases.issubset(valid_bases)
    
    def test_position_is_valid(self):
        """Test that positions are positive integers"""
        position = 500
        assert isinstance(position, int)
        assert position > 0
    
    def test_allele_is_valid_dna_base(self):
        """Test that alleles are valid DNA bases"""
        valid_alleles = {'A', 'T', 'C', 'G'}
        test_allele = 'A'
        assert test_allele.upper() in valid_alleles


# Test runner configuration
if __name__ == '__main__':
    # Run tests with verbose output
    pytest.main(['-v', __file__])