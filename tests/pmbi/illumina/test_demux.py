import pytest
import pandas as pd
from pathlib import Path
from io import StringIO, TextIOWrapper
import tempfile
import os
from pmbi.illumina.demux import read_samplesheet_data_section

@pytest.fixture
def sample_sheet():
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
    return StringIO(content)

@pytest.fixture
def empty_data_section_sample_sheet():
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
    return StringIO(content)

@pytest.fixture
def no_data_section_sample_sheet():
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
    return StringIO(content)

def test_read_samplesheet_data_section_valid(sample_sheet):
    """Test reading a valid sample sheet."""
    result = read_samplesheet_data_section(sample_sheet)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 3
    assert 'Sample_ID' in result.columns
    assert 'index' in result.columns
    assert 'index2' in result.columns
    
    # Check specific values
    assert result['Sample_ID'].tolist() == ['Sample1', 'Sample2', 'Sample3']
    assert result['index'].tolist() == ['ATCACG', 'CGATGT', 'TTAGGC']
    assert result['index2'].tolist() == ['GTCGTC', 'ACGTAC', 'TGCATG']

def test_read_samplesheet_data_section_empty_data(empty_data_section_sample_sheet):
    """Test reading a sample sheet with an empty [Data] section."""
    result = read_samplesheet_data_section(empty_data_section_sample_sheet)

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0  # Empty DataFrame but with columns

def test_read_samplesheet_data_section_no_data_section(no_data_section_sample_sheet):
    """Test reading a sample sheet without a [Data] section."""
    with pytest.raises(ValueError, match="Data section missing"):
        read_samplesheet_data_section(no_data_section_sample_sheet)

