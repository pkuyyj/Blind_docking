# Blind_docking

## Installation

```bash
conda create -n blind_docking python==3.10 numpy 
conda activate blind_docking
conda install vina==1.2.3 -c conda-forge
pip install biopandas
```

## Example

```bash
python equi_bind_diffdock_box_vina.py
```
Then results including original output poses and RMSDS will be in the `results` folder.