# NuXL to DIA-NN Spectral Library Converter

Python tool to convert NuXL DDA search results (crosslinks and linear peptides) to DIA-NN compatible spectral libraries.

## Features

- Converts NuXL TextExporter output (`.XLs.unknown` and `.peptides.unknown` files) to DIA-NN library format
- Handles RNA-protein crosslinks with proper modification notation (`[nucleotide]`)
- Converts CCS to ion mobility (1/K₀) using Mason-Schamp equation
- Supports piecewise linear iRT calibration
- Automatic handling of Oxidation + crosslink position conflicts
- 99.9% match rate with original R implementation

## Installation

```bash
pip install -r requirements.txt
```

### Dependencies

**Required:**
- pandas
- numpy

**Optional:**
- scipy (for linear iRT regression)
- pwlf (for piecewise iRT regression)
- pyopenms (for modification validation)

## Usage

### Command Line

```bash
# Basic conversion (no iRT)
python nuxl2dia.py -i data/*_XLs.unknown data/*_peptides.unknown -o library.tsv

# With piecewise iRT calibration
python nuxl2dia.py -i data/*.unknown -o library.tsv --irt piecewise --irt-ref hela_reference.tsv

# With linear iRT calibration
python nuxl2dia.py -i data/*.unknown -o library.tsv --irt linear --irt-ref hela_reference.tsv
```

### As Python Library

```python
from nuxl2dia import NuXLLibraryConverter

converter = NuXLLibraryConverter()
(converter
    .load_files(['exp1_XLs.unknown', 'exp1_peptides.unknown'])
    .filter_decoys()
    .filter_localization(min_score=0.0)
    .deduplicate()
    .convert_ccs_to_im()
    .build_modified_sequences()
    .parse_peak_annotations()
    .filter_fragment_ions(types=['b', 'y'])
    .convert_rt_to_irt(mode='piecewise', reference='hela_reference.tsv')
    .format_output()
    .save('output_library.tsv'))
```

## Output Format

The converter produces a 15-column TSV file compatible with DIA-NN:

| Column | Description |
|--------|-------------|
| ModifiedPeptideSequence | Sequence with modifications in `[nucleotide]` and `(UniMod:XX)` format |
| PrecursorCharge | Precursor charge state |
| AverageExperimentalRetentionTime | Retention time in seconds |
| PrecursorIonMobility | Ion mobility (1/K₀) in Vs/cm² |
| PeptideSequence | Stripped sequence (no modifications) |
| PrecursorMz | Precursor m/z |
| ProteinId | Protein accession(s) |
| Annotation | Fragment annotation |
| FragmentSeriesNumber | Fragment ion number |
| ProductMz | Fragment m/z |
| LibraryIntensity | Fragment intensity |
| FragmentCharge | Fragment charge |
| FragmentType | b or y |
| NormalizedRetentionTime | iRT value (if calibrated) |
| FragmenLossType | Loss type annotation |

## Comparison with R Implementation

This Python implementation replaces three R scripts:
- `UVECO_XL_lib_shot130925.R` - Crosslink processing
- `UVECO_PE_lib_shot130925.R` - Linear peptide processing
- `UV_iRT_lm.R` - iRT calibration

### Key Improvements

| Issue in R | Python Solution |
|------------|-----------------|
| Column access by index (fragile) | Named column access |
| Manual corrections needed (lines 175-178) | Automatic using position column |
| Complex offset calculations with bugs | Direct position-based insertion |
| Hardcoded Windows paths | Parameterized paths |
| No error handling | Comprehensive logging |
| Three separate scripts | Single module with fluent API |

## Reference R Scripts

The original R scripts are included for reference:
- `UVECO_XL_lib_shot130925.R`
- `UVECO_PE_lib_shot130925.R`
- `UV_iRT_lm.R`

## DIA-NN Modification Declaration

For DIA-NN to recognize the nucleotide modifications, add the contents of `declare_mod_UV_UCGA_ext.txt` to your DIA-NN modification settings.

## License

MIT License
