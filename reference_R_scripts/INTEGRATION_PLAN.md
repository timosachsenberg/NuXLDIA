# NuXL DIA Library Generation - Integration Plan

## Overview

Three R scripts convert NuXL DDA search results into a spectral library format compatible with DIA-NN v2.1.0.

**Pipeline Flow:**
```
NuXL TextExporter Output
        │
        ├──> Script 1: UVECO_XL_lib_shot130925.R (crosslinks)
        │                    │
        ├──> Script 2: UVECO_PE_lib_shot130925.R (linear peptides)
        │                    │
        │                    ▼
        │            Combined Library (UVECO_SHOT_SPLIBRARY_130126.tsv)
        │                    │
        │                    ▼
        └──> Script 3: UV_iRT_lm.R
                             │
                             ▼
              Final Library with iRT (UVECO_SHOT_LIB_iRT_130126.tsv)
              + DIA-NN/declare_mod_UV_UCGA_ext.txt (DIA-NN mod definitions)
```

---

## Script 1: UVECO_XL_lib_shot130925.R (Crosslinks)

### Input
- `*.XLs.unknown` files from NuXL TextExporter (TSV format)
- Key columns used: RT, mz, sequence, charge, accessions, peak_annotations, CCS, best_localization_position, localization_score, NuXL:score, NuXL:NA (nucleotide)

### Processing Steps
1. **Combine** all XL files from multiple runs
2. **Rename columns** to meaningful names (hardcoded column indices)
3. **Filter** by `locscore > 0`
4. **Deduplicate** by ID (sequence + nucleotide + charge), keeping best NuXL score then best loc score
5. **Convert CCS to Ion Mobility** using Mason-Schamp equation:
   ```r
   IM = CCS * 2.6867811e25 / (charge * 3/16 * 1e24 * sqrt(2*pi / (reduced_mass * kB * T)) * e)
   ```
   Where reduced_mass uses N2 (28.0134 Da), T=305K
6. **Insert nucleotide modification** into peptide sequence at localization position
7. **Reformat peak annotations** (swap delimiters: `|` → `;`, `,` → `|`)
8. **Expand to fragment rows** with columns: ModifiedPeptide, StrippedPeptide, PrecursorCharge, RT, IonMobility, ProteinID, PrecursorMz, FragmentMz, RelativeIntensity, FragmentCharge, FragmentAnnotation, FragmentType
9. **Filter** to b/y ions only

### MANUAL STEP (Lines 175-178)
```r
#MANUALLY CORRECT 3 LINES
uv_modo_nok [1,12] <- "ATLGEVGNAEH(U)M(Oxidation)LR"
uv_modo_nok [4,12] <- "ATLGEVGNAEH(C-H3N1)M(Oxidation)LR"
uv_modo_nok [5,12] <- "ATLGEVGNAEH(CG-H3N1-H2O1)M(Oxidation)LR"
```

**Root Cause:** The insertion logic fails when:
- Peptide has existing modification (e.g., `M(Oxidation)`)
- Crosslink position is BEFORE the existing modification
- The offset calculation adds `+11` (length of "(Oxidation)") but this is incorrect when inserting before it

**Fix Pattern:** These are cases where `local + 11` is used but the crosslink site (position 11) is actually before the M(Oxidation) at position 12. The nucleotide should be inserted at position 11, not at position 22.

### Output
- `UVECO_XL_lib_shot130925.tsv`

---

## Script 2: UVECO_PE_lib_shot130925.R (Linear Peptides)

### Input
- `*.peptides.unknown` files from NuXL TextExporter

### Processing Steps
Same as Script 1 but simpler:
1. Combine files
2. Rename columns
3. Convert CCS → IM
4. Deduplicate by best score
5. Expand to fragment rows
6. Filter to b/y ions

### Manual Steps
**NONE** - linear peptides don't have nucleotide insertions

### Output
- `UVECO_PE_lib_shot130925.tsv`
- Combined: `UVECO_SHOT_SPLIBRARY_130126.tsv`

---

## Script 3: UV_iRT_lm.R (iRT Conversion)

### Input
- HeLa QC runs from FragPipe: `psm.tsv`, `library.tsv`
- Combined library: `UVECO_SHOT_SPLIBRARY_130126.tsv`

### Processing Steps
1. **Load HeLa reference data** from FragPipe output
2. **Join** PSM retention times with library iRT values
3. **Fit piecewise linear regression** with 3 breakpoints using `segmented` package
4. **Apply model** to predict iRT for all library entries
5. **Format conversions** for DIA-NN compatibility:
   - `(C...)` → `[C...]` (nucleotide mods to square brackets)
   - `(U...)` → `[U...]`
   - `(A...)` → `[A...]`
   - `(G...)` → `[G...]`
   - `(Oxidation)` → `(UniMod:35)`
6. **Extract fragment series number** from annotation
7. **Rename columns** to DIA-NN expected names

### Manual Steps
**NONE** - but requires external HeLa data

### Output
- `UVECO_SHOT_LIB_iRT_130126.tsv`

---

## External Dependencies

### Required Files
1. **template_lib.tsv** - Template for library format (not in repo, needs to be provided or created)
2. **HeLa FragPipe output** - `psm.tsv` and `library.tsv` for iRT calibration
3. **DIA-NN/declare_mod_UV_UCGA_ext.txt** - DIA-NN modification declarations (provided)

### R Packages
```r
tidyverse, tidyr, dplyr, readr, ggplot2, ggpubr, arsenal, stringr, data.table, segmented
```

---

## Integration Plan for NuXL Webapp

### Option A: Keep R Backend (Recommended for quick integration)

1. **Package as single R script** with functions:
   ```r
   process_crosslinks(input_dir, output_file)
   process_peptides(input_dir, output_file)
   combine_libraries(xl_file, pep_file, output_file)
   add_irt(library_file, hela_dir, output_file)
   ```

2. **Automate the manual step** - Replace lines 175-178 with logic:
   ```r
   # Fix for peptides with existing modifications where XL position < mod position
   fix_mod_insertion <- function(sequence, nucleotide, local, existing_mods) {
     # Count modification characters before insertion point
     # Insert at correct position accounting for existing mods
   }
   ```

3. **Generate template_lib.tsv** programmatically (just needs column headers)

4. **Handle iRT options:**
   - Option 1: User provides HeLa reference runs
   - Option 2: Use AlphaPeptDeep for iRT prediction (as Timo suggested)
   - Option 3: Skip iRT (some DIA-NN workflows don't need it)

### Option B: Rewrite in Python (for better webapp integration)

If the webapp is Python-based, port the R logic to pandas/numpy. The core operations are straightforward data manipulation.

---

## Automating the Manual Step

The manual corrections at lines 175-178 follow a predictable pattern. Here's the automated fix:

```r
# Current problematic code:
for(i in 1:nrow(uv_modo_nok)) {
  uv_modo_nok[i,12] <- fun_insert(x = uv_modo_nok[i,4],
                              pos = uv_modo_nok[i,11] + 11,  # <-- BUG: +11 is wrong
                              insert = uv_modo_nok[i,6])
}

# Fixed version:
for(i in 1:nrow(uv_modo_nok)) {
  seq <- uv_modo_nok[i,4]  # e.g., "ATLGEVGNAEHM(Oxidation)LR"
  local <- uv_modo_nok[i,11]  # insertion position (1-indexed)
  nucleotide <- uv_modo_nok[i,6]  # e.g., "(U)"

  # Find position of existing modification
  mod_start <- regexpr("\\(", seq)

  # If insertion point is before the modification, don't add offset
  if (local < mod_start) {
    uv_modo_nok[i,12] <- fun_insert(x = seq, pos = local, insert = nucleotide)
  } else {
    # Calculate actual offset based on modification length
    mod_match <- regmatches(seq, regexpr("\\([^)]+\\)", seq))
    offset <- nchar(mod_match)
    uv_modo_nok[i,12] <- fun_insert(x = seq, pos = local + offset, insert = nucleotide)
  }
}
```

---

## Key Column Mappings

### NuXL TextExporter → Library Format

| NuXL Column (index) | Renamed To | Final Library Column |
|---------------------|------------|---------------------|
| #rt (1) | RT | AverageExperimentalRetentionTime |
| IM (29) | CCS → IM | PrecursorIonMobility |
| mz (10) | mz | PrecursorMz |
| sequence (5) | sequence | PeptideSequence |
| charge (6) | charge | PrecursorCharge |
| accessions (11) | accessions | ProteinId |
| peak_annotations (16) | peak_ann | → FragmentMz, LibraryIntensity, FragmentCharge, Annotation |
| NuXL:NA (36) | nucleic | → ModifiedPeptideSequence |
| NuXL:best_localization_position (47) | loc | (used for insertion) |
| NuXL:best_localization_score (48) | locscore | (used for filtering) |
| NuXL:score (68) | nuxl_score | (used for ranking) |

---

## Summary of Manual Steps to Automate

| Script | Line(s) | Manual Step | Automation Strategy |
|--------|---------|-------------|---------------------|
| UVECO_XL_lib_shot130925.R | 175-178 | Correct 3 peptide sequences | Fix insertion logic to handle pre-existing modifications |
| UV_iRT_lm.R | N/A | Requires HeLa QC runs | Provide option for AlphaPeptDeep or user-supplied reference |

---

## Recommended Next Steps

1. **Create template_lib.tsv** - Single header row with required columns
2. **Fix the insertion bug** in Script 1 (see automated fix above)
3. **Parameterize paths** - Replace hardcoded Windows paths with arguments
4. **Add iRT prediction option** - Integrate AlphaPeptDeep or make HeLa reference optional
5. **Wrap in single entry point** - One function/script that runs the full pipeline
6. **Add input validation** - Check for required columns in NuXL output
