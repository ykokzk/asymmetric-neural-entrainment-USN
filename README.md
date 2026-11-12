# Asymmetric neural entrainment at resonance frequencies underlies unilateral spatial neglect

This repository contains the source data and analysis code for:

Yuka O. Okazaki, Noriaki Hattori, Teiji Kawano, Megumi Hatakenaka, Ichiro Miyai, and Keiichi Kitajo. *Asymmetric neural entrainment at resonance frequencies underlies unilateral spatial neglect.* eLife. <https://doi.org/10.7554/eLife.107993>

The study recorded steady-state visual evoked potentials (SSVEPs; 3–30 Hz flicker) in stroke patients with unilateral spatial neglect (USN), stroke patients without USN, and healthy controls.

## Contents

```
.
├── .gitattributes
├── README.md
├── CITATION.cff
├── LICENSE
├── MANIFEST_public_release.txt
├── source_data/
│   ├── LICENSE
│   ├── Table1_source_data.xlsx
│   ├── Figure2_source_data.xlsx
│   ├── Figure3_source_data.xlsx
│   ├── Figure4_source_data.xlsx
│   ├── Figure5_source_data.xlsx
│   ├── Figure6_source_data.xlsx
│   ├── Figure7_source_data.xlsx
│   ├── Figure8_source_data.xlsx
│   └── Figure9_source_data.xlsx
├── transfer_entropy/
│   ├── run_transfer_entropy_analysis.m
│   └── TE.m
└── model/
    └── run_kuramoto_model.m
```

### Source data

`source_data/` contains processed source data underlying Table 1 and main Figures 2–9, with one `.xlsx` file per table or figure.

`Table1_source_data.xlsx` contains three worksheets:

| Worksheet | Contents |
| --- | --- |
| `Participant_values` | De-identified processed participant-level values underlying Table 1, one row per participant, identified by `subject_id` |
| `Table1_values` | Group-level summary statistics and p-values reported in Table 1 |
| `Notes` | Variable definitions and notes on the data |

Missing and non-applicable values are given as `NA`.

For figures involving outlier removal, both the plotted values and the values before outlier removal are provided in separate worksheets (`Figure_values` / `All_values`, or `A_Figure_values` / `A_All_values`). Excluded values are marked `Outlier` in the `Figure_values` worksheets.

### Transfer entropy analysis

```
transfer_entropy/
    run_transfer_entropy_analysis.m
    TE.m
```

`run_transfer_entropy_analysis.m` computes transfer entropy between the amplitude envelopes of band-pass filtered EEG signals for every channel pair and delay. The script implements the transfer entropy computation used in the study. 

`TE.m` is the transfer entropy class used by the script and was provided by H. Yokoyama. It estimates the required probability distributions using histograms.

`TE.m` contains the function `histndim`, taken from VOICEBOX, a MATLAB toolbox by Mike Brookes. The original copyright and GPL licence notices are preserved in `TE.m`. Accordingly, the MATLAB code in this repository is distributed under GPL-2.0-or-later.

### Computational model

```
model/
    run_kuramoto_model.m
```

`run_kuramoto_model.m` implements the coupled-oscillator (Kuramoto) model reported in Figure 9. It runs the asymmetric interhemispheric coupling, no-stimulation, asymmetric intrahemispheric coupling, and different-natural-frequency conditions and computes the order parameter of the left and right hemispheres.

The script uses MATLAB base functions only and requires no external input file.

The worksheets of `Figure9_source_data.xlsx` correspond to the simulation conditions as follows:

| Worksheet | Panel | Condition |
| --- | --- | --- |
| `A` | Figure 9A | `interhemispheric` and no-stimulation control `nostim` |
| `B` | Figure 9B | `intrahemispheric` |
| `C` | Figure 9C | `natural_frequency` |

## Requirements

The versions below are those used to run and check the code in this repository. Other versions may also work but have not been tested.

| Requirement | Version used | Needed by |
| --- | --- | --- |
| MATLAB | 25.2 (R2025b) Update 2 | all code |
| Signal Processing Toolbox | 25.2 (R2025b) | `run_transfer_entropy_analysis.m` |
| FieldTrip | git revision `da1b29fa4bbd3b48a43478e8f139df8725cbb209` (2025-02-11), 11 commits after release tag `20250123` | `run_transfer_entropy_analysis.m` |

FieldTrip is not included in this repository. Download it from <https://www.fieldtriptoolbox.org/> and add it to the MATLAB path before running the transfer entropy script.

`run_transfer_entropy_analysis.m` uses FieldTrip functions including `ft_preproc_bandpassfilter` and `nearest`.

`TE.m` and `run_kuramoto_model.m` require MATLAB only.

## Usage

### Transfer entropy analysis

`run_transfer_entropy_analysis.m` is a script, not a function. Load preprocessed EEG data into the workspace and then run:

```matlab
addpath('transfer_entropy');
addpath('<path to fieldtrip>'); ft_defaults;

data.segdata.<condition> = ...   % [channels x time x trials]
data.timeinfo            = ...   % [1 x time], seconds, 0 s at stimulus onset

run_transfer_entropy_analysis;
```

`data.segdata` contains one field per experimental condition. The first 19 channels are analysed, in the channel order of the recording montage. The original study used a 19-channel 10–20 montage sampled at 500 Hz.

The input should contain broadband EEG; band-pass filtering and Hilbert-envelope calculation are performed by the script. `data.timeinfo` is used to select the analysed segment (−10 to 20 s).

The result is stored in `output_data_te`. `TEy2x` and `TEx2y` contain transfer entropy values as `[epoch x channel pair x delay x condition]`, averaged over trials.

The participant-level EEG recordings from the original study are not included in this repository. The script is provided so that the computation can be inspected and applied to other data.

### Computational model

```matlab
run('model/run_kuramoto_model.m');
```

No input data are required. Results are stored in the struct `results`, with `rho_L` and `rho_R` for each simulation condition.

The model draws natural frequencies and initial phases at random and does not use a fixed random seed, as in the original analysis. Values therefore vary slightly between runs unless a seed is set explicitly. The values in `Figure9_source_data.xlsx` are those of the original run reported in the paper.

## Data availability

This repository contains de-identified processed source data underlying Table 1 and main Figures 2–9, together with the code used for the transfer entropy analysis and computational simulations.

Raw EEG recordings, original clinical records, and direct participant identifiers are not included.

## Licences

This repository uses separate licences for code and source data:

| Path | Licence |
| --- | --- |
| `transfer_entropy/`, `model/` | GPL-2.0-or-later — see [`LICENSE`](LICENSE) |
| `source_data/` | CC BY 4.0 — see [`source_data/LICENSE`](source_data/LICENSE) |

The MATLAB code is distributed under GPL-2.0-or-later. `TE.m` incorporates GPL-licensed code from VOICEBOX, with the original copyright and licence notices preserved.

The source-data files are distributed under the Creative Commons Attribution 4.0 International licence.

FieldTrip is not distributed in this repository and retains its own licence.

## Citation

If you use this material, please cite the paper:

> Okazaki, Y. O., Hattori, N., Kawano, T., Hatakenaka, M., Miyai, I., & Kitajo, K.  
> Asymmetric neural entrainment at resonance frequencies underlies unilateral spatial neglect. *eLife*. <https://doi.org/10.7554/eLife.107993>

See [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.
