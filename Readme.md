# CCCN: Ambient Noise Cross-Correlation

CCCN is a Python toolkit for ambient noise cross-correlation, preprocessing, and stacking of correlation functions.

## Features

- Read waveform data (e.g., SAC) using ObsPy patterns.
- Preprocess traces for noise correlation:
  - trim by reference window (`day`/`hour`/`minute`)
  - detrend
  - time-domain normalization (running absolute mean)
  - spectral whitening
- Compute station-pair cross-correlation with MPI parallelism.
- Save outputs as per-pair HDF5 correlation files.
- Post-process and stack correlation functions:
  - linear stack
  - phase-weighted stack (PWS)
  - selective stack

## Installation

### Option 1: `pip` in a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
```

### Option 2: Conda

```bash
conda create -n cccn python=3.11 obspy mpi4py h5py scipy pyyaml -c conda-forge
conda activate cccn
pip install -e .
```

## Quick Start

```python
from cccn.noise import CrossCorrelation

cc = CrossCorrelation(level="INFO")

# Basic parameters
cc.para.datapath = "example/dataSAC/2008.010"
cc.para.outpath = "./CC_OUT"
cc.para.suffix = "sac"
cc.para.timeduration = 86400
cc.para.target_dt = 0.1
cc.para.freqmin = 0.02
cc.para.freqmax = 0.2
cc.para.maxlag = 500
cc.para.reftime = "day"   # one of: day, hour, minute

cc.read_data(matchstr="*")
cc.perwhiten()
cc.run_cc()
cc.finalize()
```

## Run with MPI

For parallel processing on multiple ranks:

```bash
mpiexec -n 4 python your_script.py
```

## Main Parameters (`Para`)

- `datapath`: input path or wildcard-compatible path
- `outpath`: output directory for correlation files
- `freqmin`, `freqmax`: whitening band (Hz)
- `nsmooth`: normalization smoothing half window (`0` means one-bit normalization)
- `timeduration`: expected segment duration (seconds)
- `cut_precentatge`: fraction trimmed at both sides before processing
- `suffix`: file suffix (e.g., `sac`, `SAC`)
- `target_dt`: target sampling interval (seconds), `None` keeps original
- `reftime`: reference alignment (`day`, `hour`, `minute`)
- `maxlag`: max lag time for correlation output (seconds)
- `src_mask`: optional source station list for partial pairing

## Output

Cross-correlation results are written to HDF5 files under `outpath`, named like:

```text
COR_<NET1_STA1>_<NET2_STA2>_<CHPAIR>.h5
```

Each file contains:

- file attributes: station/channel metadata, `delta`, `lag`
- datasets: one dataset per time tag (`YYYYMMDDHHMMSS`)

## Stacking and Post-processing

```python
from cccn.post_func import PostProcForNoise

pp = PostProcForNoise("your_noise_file.h5")
tr = pp.stack_all(method="linear", normalize=True, add_sac_header=True)
```

Supported stack methods in `stack_all(method=...)`:

- `linear`
- `pws`
- `selective`

## Notes

- Ensure all traces to be stacked have consistent sample count and sampling interval.
- `mpi4py` requires a working MPI runtime (e.g., OpenMPI or MPICH).
- If your data contains duplicates with identical `network_station_channel`, duplicates are removed automatically during reading.
