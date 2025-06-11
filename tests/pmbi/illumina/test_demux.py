import pytest
import pandas as pd
from pathlib import Path
import tempfile
import os
from pmbi.illumina.demux import read_samples_from_samplesheet

@pytest.fixture
def sample_sheet_file():
    """Create a temporary sample sheet file for testing."""
    content = """[Header]
IEMFileVersion,5
Date,2023-06-11

[Reads]
151
151

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,index,index2,Sample_Project
Sample1,Sample1_Name,ATCACG,GTCGTC,Project1
Sample2,Sample2_Name,CGATGT,ACGTAC,Project1
Sample3,Sample3_Name,TTAGGC,TGCATG,Project1
"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        f.write(content)
        temp_file = f.name
    
    yield temp_file
    
    # Clean up the temporary file
    os.unlink(temp_file)

@pytest.fixture
def empty_data_section_sample_sheet_file():
    """Create a temporary sample sheet file with an empty [Data] section."""
    content = """[Header]
IEMFileVersion,5
Date,2023-06-11

[Reads]
151
151

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,index,index2,Sample_Project
"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        f.write(content)
        temp_file = f.name
    
    yield temp_file
    
    # Clean up the temporary file
    os.unlink(temp_file)

@pytest.fixture
def no_data_section_sample_sheet_file():
    """Create a temporary sample sheet file without a [Data] section."""
    content = """[Header]
IEMFileVersion,5
Date,2023-06-11

[Reads]
151
151

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        f.write(content)
        temp_file = f.name
    
    yield temp_file
    
    # Clean up the temporary file
    os.unlink(temp_file)

def test_read_samples_from_samplesheet_valid(sample_sheet_file):
    """Test reading a valid sample sheet."""
    result = read_samples_from_samplesheet(sample_sheet_file)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 3
    assert 'Sample_ID' in result.columns
    assert 'index' in result.columns
    assert 'index2' in result.columns
    
    # Check specific values
    assert result['Sample_ID'].tolist() == ['Sample1', 'Sample2', 'Sample3']
    assert result['index'].tolist() == ['ATCACG', 'CGATGT', 'TTAGGC']
    assert result['index2'].tolist() == ['GTCGTC', 'ACGTAC', 'TGCATG']

def test_read_samples_from_samplesheet_empty_data(empty_data_section_sample_sheet_file):
    """Test reading a sample sheet with an empty [Data] section."""
    result = read_samples_from_samplesheet(empty_data_section_sample_sheet_file)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0  # Empty DataFrame but with columns

def test_read_samples_from_samplesheet_no_data_section(no_data_section_sample_sheet_file):
    """Test reading a sample sheet without a [Data] section."""
    with pytest.raises(ValueError, match="Data section missing"):
        read_samples_from_samplesheet(no_data_section_sample_sheet_file)

def test_read_samples_from_samplesheet_nonexistent_file():
    """Test reading a nonexistent sample sheet file."""
    with pytest.raises(FileNotFoundError):
        read_samples_from_samplesheet("nonexistent_file.csv")
