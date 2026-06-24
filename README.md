This project contains Jupyter Notebooks for comparing peak coordinates and signals.

## Instructions

#### a) Create virtual environment from requirements
```bash
python -m venv .venv
source venv/bin/activate
pip install -r requirements.txt 
```

#### b) Load data in convenient folder structure
I recommend loading .bed and .wig files using the following structure

ChEC_Data/
├── homer_peaks_merged/
│   └── sample1.bed
├── homer_peaks_processed/
│   ├── sample1_A.bed
│   ├── sample1_B.bed
│   └── sample1_C.bed
├── wig/
│   ├── cpm/
│   │   ├── sample1_cpm_A.bed
│   │   ├── sample1_cpm_B.bed
│   │   └── sample1_cpm_C.bed
│   ├── cpm_mean/
│   │   └── sample1_cpm_mean.bed
│   ├── spikein/
│   │   ├── sample1_spikein_A.bed
│   │   ├── sample1_spikein_B.bed
│   │   └── sample1_spikein_C.bed
│   └── spikein_mean/
│       └── sample1_spikein_mean.bed

#### c) Run notebooks/compare_coordinates/1_compare_peaks.ipynb

#### d) Run notebooks/compare_signals/2_calculate_signals.ipynb

#### e) Plot signals notebooks/compare_signals/3_heatmaps.ipynb notebooks/compare_signals/3_scatterplot.ipynb