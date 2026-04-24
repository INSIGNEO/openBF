# Configuration

The configuration file is a YAML document describing the vascular network, blood properties, and solver parameters. You can also use the [openBF webapp](https://openBF.streamlit.app) to build config files interactively.

---

## Top-level keys

```yaml
project_name: my_simulation    # used to name the results folder
write_results: ["P", "Q"]      # quantities to save: P, Q, u, A
output_directory: path/to/out  # optional; default is ./<project_name>_results
inlet_file: my_inlet.dat       # optional; default is <project_name>_inlet.dat
```

The inlet file is a two-column ASCII file with time (s) in column 1 and volumetric flow rate (m³/s) in column 2. The inlet vessel is the one whose source node (`sn`) is 1.

---

## `solver`

```yaml
solver:
  Ccfl: 0.9          # Courant number; keep ≤ 1.0
  cycles: 100        # maximum cardiac cycles before stopping
  convergence_tolerance: 5.0  # RMSE threshold in mmHg between consecutive cycles
  jump: 100          # number of time-points saved per cardiac cycle
```

The simulation stops when the pressure RMSE between two consecutive cycles drops below `convergence_tolerance`, or when `cycles` is reached.

---

## `blood`

```yaml
blood:
  rho: 1060.0   # density [kg/m³]
  mu:  0.004    # dynamic viscosity [Pa·s]
```

---

## `network`

A list of vessel entries. Each vessel connects a source node (`sn`) to a target node (`tn`). Node 1 is the network inlet.

### Mandatory parameters

```yaml
- label: aorta    # used to name output files
  sn: 1           # source node (integer)
  tn: 2           # target node (integer)
  L: 0.4          # length [m]
  E: 400000.0     # wall Young's modulus [Pa]
  R0: 0.015       # constant lumen radius [m]
  # — or, for a linearly tapered vessel:
  Rp: 0.016       # proximal radius at sn [m]
  Rd: 0.014       # distal radius at tn [m]
```

### Optional parameters

```yaml
  M: 40              # spatial divisions (default: L/M = 1 mm, minimum 5)
  h0: 0.0015         # wall thickness [m] (default: computed from R0)
  Pext: 0.0          # external pressure [Pa] (default: 0)
  gamma_profile: 9   # velocity profile parameter (default: 2, parabolic)
  to_save: true      # include vessel in output files (default: true)
  initial_pressure: 0.0  # initial pressure [Pa] (default: 0)
  initial_flow: 0.0      # initial flow [m³/s] (default: 0)
  visco-elastic: false   # enable visco-elastic wall model (default: false)
```

---

## Outlet boundary conditions

A vessel with no downstream connections requires an outlet BC. Three options:

### Reflection coefficient

```yaml
  Rt: 0.0    # -1.0 ≤ Rt ≤ 1.0; 0 = fully absorbing, 1 = fully reflecting
```

### Two-element Windkessel (WK2)

```yaml
  R1: 1.0e8   # peripheral resistance [Pa·s/m³]
  Cc: 1.0e-9  # peripheral compliance [m³/Pa]
```

### Three-element Windkessel (WK3)

WK3 is activated when `R2` is provided. `R1` here is the proximal (characteristic) impedance.

```yaml
  R1: 1.0e7              # proximal impedance [Pa·s/m³]
  R2: 9.0e8              # peripheral resistance [Pa·s/m³]
  Cc: 1.0e-9             # peripheral compliance [m³/Pa]
  inlet_impedance_matching: false  # auto-set R1 to match wave impedance (default: false)
  Pout: 0.0              # Windkessel outlet pressure [Pa] (default: 0)
```

---

## Full template

```yaml
project_name: <name>
write_results: ["P"]           # any subset of ["P", "Q", "u", "A"]
output_directory: <path>       # optional
inlet_file: <name>_inlet.dat   # optional

solver:
  Ccfl: 0.9
  cycles: 100
  convergence_tolerance: 5.0
  jump: 100

blood:
  rho: 1060.0
  mu: 0.004

network:
  - label: <name>
    sn: <int>
    tn: <int>
    L: <float>          # [m]
    E: <float>          # [Pa]
    R0: <float>         # [m]  — or use Rp/Rd for tapering
    M: <int>            # optional
    h0: <float>         # [m], optional
    Pext: <float>       # [Pa], optional
    gamma_profile: <float>   # optional, default 2
    to_save: true            # optional
    # outlet (choose one):
    Rt: <float>
    # or wk2:
    R1: <float>
    Cc: <float>
    # or wk3:
    R1: <float>
    R2: <float>
    Cc: <float>
    Pout: <float>            # optional
    inlet_impedance_matching: false  # optional
```
