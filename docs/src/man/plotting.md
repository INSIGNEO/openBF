# Plotting

## Python (matplotlib)

The `models/boileau2015/` directory contains `plot.py` and `plot_style.txt` for reproducing the benchmark figures. A model-specific script is also available for adan56 (`models/boileau2015/adan56/plot.py`).

To use them, run the simulation first, then call the plotting script from the model directory:

```bash
cd models/boileau2015/adan56
python plot.py
```

## Julia

Results can be loaded with `DelimitedFiles` and plotted with any Julia plotting library.

```julia
using DelimitedFiles
using Plots  # or CairoMakie, GLMakie, etc.

# load a .last file: rows = time-points, cols = [time, inlet, 1/4, mid, 3/4, outlet]
P = readdlm("adan56_results/aorta_P.last")

t = P[:, 1]           # time [s]
p_mid = P[:, 4] / 133.322  # mid-point pressure converted to mmHg

plot(t, p_mid, xlabel="t (s)", ylabel="P (mmHg)", label="aorta mid-point")
```

## Unit conversions

| Quantity | SI unit | Common display unit | Factor |
|----------|---------|---------------------|--------|
| Pressure | Pa | mmHg | ÷ 133.322 |
| Flow | m³/s | ml/s | × 10⁶ |
| Velocity | m/s | cm/s | × 100 |
| Area | m² | mm² | × 10⁶ |
| Radius | m | mm | × 1000 |
