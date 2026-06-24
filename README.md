This project contains Jupyter Notebooks for comparing peak coordinates and signals.

## Instructions
#### a) Load data in convenient folder structure
I recommend loading .bed and .wig files using the following structure
```bash
ChEC_Data/
├── homer_peaks_merged/
│   └── sample1.bed
├── wig/
│   ├── cpm_mean/
│   │   └── sample1_cpm_mean.wig
│   └── spikein_mean/
│       └── sample1_spikein_mean.wig
```
#### b) Download version controlled repository
```bash
git clone https://github.com/plantazja/ChEC_Analysis
cd ChEC_Analysis
```

#### c) Create virtual environment from requirements.txt
```bash
python3 -m venv .venv
source .venv/bin/activate
pip3 install -r requirements.txt 
```

#### d) Run notebooks/compare_coordinates/1_compare_peaks.ipynb

#### e) Run notebooks/compare_signals/2_calculate_signals.ipynb

#### f) Plot signals notebooks/compare_signals/3_heatmaps.ipynb notebooks/compare_signals/3_scatterplot.ipynb