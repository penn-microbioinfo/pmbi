import pytest
import pandas as pd
from pathlib import Path
from io import StringIO
import tempfile
from pmbi.illumina.demux import read_samplesheet_data_section, fastq_missing_for, get_run_id

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

@pytest.fixture
def temp_rundir_with_parameters():
    """Create a temporary run directory with RunParameters.xml file."""
    with tempfile.TemporaryDirectory() as temp_dir:
        run_dir = Path(temp_dir)
        # Create RunParameters.xml file with test data
        params_xml = run_dir / "RunParameters.xml"
        xml_content = """<?xml version="1.0"?>
<RunParameters xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <Setup>
    <ExperimentName>TestExperiment</ExperimentName>
  </Setup>
  <RunId>250817_B02389_0872_BJKT7LMXYZ</RunId>
  <BaseSpaceRunId>315628943</BaseSpaceRunId>
  <InstrumentName>B02389</InstrumentName>
</RunParameters>
"""
        with open(params_xml, 'w') as f:
            f.write(xml_content)
        
        yield run_dir

def test_get_run_id(temp_rundir_with_parameters):
    """Test get_run_id function with default tag."""
    # Test with default tag (BaseSpaceRunId)
    run_id = get_run_id(temp_rundir_with_parameters)
    assert run_id == "315628943"
    
    # Test with custom tag
    custom_run_id = get_run_id(temp_rundir_with_parameters, tag="RunId")
    assert custom_run_id == "250817_B02389_0872_BJKT7LMXYZ"

def test_get_run_id_file_not_found():
    """Test get_run_id function when RunParameters.xml doesn't exist."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with pytest.raises(FileNotFoundError):
            get_run_id(Path(temp_dir))

def test_get_run_id_missing_tag(temp_rundir_with_parameters):
    """Test get_run_id function when the specified tag doesn't exist."""
    with pytest.raises(ValueError, match="Unable to find NonExistentTag in RunParameters.xml"):
        get_run_id(temp_rundir_with_parameters, tag="NonExistentTag")


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

# %%

@pytest.fixture
def temp_fastq_dir_with_all_samples():
    """Create a temporary directory with FASTQ files for all samples."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create fastq files for all samples in the sample sheet
        sample_ids = ['Sample1', 'Sample2', 'Sample3']
        for sample_id in sample_ids:
            # Create R1 and R2 fastq files
            for read in ['R1', 'R2']:
                fastq_path = Path(temp_dir) / f"{sample_id}_S1_L001_{read}_001.fastq.gz"
                fastq_path.touch()
        
        yield Path(temp_dir)

@pytest.fixture
def temp_fastq_dir_with_missing_samples():
    """Create a temporary directory with missing FASTQ files for some samples."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create fastq files for only the first sample
        sample_id = 'Sample1'
        for read in ['R1', 'R2']:
            fastq_path = Path(temp_dir) / f"{sample_id}_S1_L001_{read}_001.fastq.gz"
            fastq_path.touch()
        yield Path(temp_dir)

def test_fastq_missing_for_all_present(sample_sheet, temp_fastq_dir_with_all_samples):
    """Test when all expected FASTQ files are present."""
    # Reset the sample sheet to the beginning
    sample_sheet.seek(0)
    
    # Check for missing files
    missing = fastq_missing_for(temp_fastq_dir_with_all_samples, sample_sheet)
    
    # Assert no samples are missing
    assert isinstance(missing, list)
    assert len(missing) == 0

def test_fastq_missing_for_some_missing(sample_sheet, temp_fastq_dir_with_missing_samples):
    """Test when some expected FASTQ files are missing."""
    # Reset the sample sheet to the beginning
    sample_sheet.seek(0)
    
    # Check for missing files
    missing = fastq_missing_for(temp_fastq_dir_with_missing_samples, sample_sheet)
    
    # Assert the correct samples are missing
    assert isinstance(missing, list)
    assert len(missing) == 2
    assert 'Sample2' in missing
    assert 'Sample3' in missing

