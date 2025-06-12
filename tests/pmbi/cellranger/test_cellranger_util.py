import os
import tempfile
import pandas as pd
import pytest
from pathlib import Path

from pmbi.cellranger.util import read_10x_index_sheets

def test_read_10x_index_sheets():
    # Create temporary test data
    test_data = pd.DataFrame({
        'index_name': ['idx1', 'idx2', 'idx3'],
        'index(i7)': ['ATCGATCGAT', 'GCTAGCTAGT', 'CGATCGATCG'],
        'index2_workflow_a(i5)': ['TCGATCGATC', 'ATCGATCGAT', 'TAGCTAGCTA'],
        'index2_workflow_b(i5)': ['CGATCGATCG', 'TCGATCGATC', 'ATCGATCGAT']
    })
    
    # Create temporary files
    with tempfile.TemporaryDirectory() as tmpdir:
        file1_path = Path(tmpdir) / "test_indices1.csv"
        file2_path = Path(tmpdir) / "test_indices2.csv"
        
        # Write same data to both files
        test_data.to_csv(file1_path, index=False)
        test_data.to_csv(file2_path, index=False)
        
        # Test workflow A
        result_a = read_10x_index_sheets(file1_path, file2_path, workflow='a')
        assert len(result_a) == 6  # Should have 3 rows from each file
        assert list(result_a.columns) == ['index_name', 'i7', 'i5_workflow_a']
        
        # Test workflow B
        result_b = read_10x_index_sheets(file1_path, file2_path, workflow='b')
        assert len(result_b) == 6
        assert list(result_b.columns) == ['index_name', 'i7', 'i5_workflow_b']
        
        # Test invalid workflow
        with pytest.raises(ValueError, match="Invalid workflow:.*"):
            read_10x_index_sheets(file1_path, workflow='c')
        
        # Test missing columns
        bad_data = pd.DataFrame({
            'index_name': ['idx1'],
            'wrong_column': ['ATCGATCGAT']
        })
        bad_file = Path(tmpdir) / "bad_indices.csv"
        bad_data.to_csv(bad_file, index=False)
        
        with pytest.raises(ValueError, match="Input sheet missing expected column names:.*"):
            read_10x_index_sheets(bad_file, workflow='a')
