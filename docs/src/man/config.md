Here is a description of the config file and its entries. You can also use openBF [webapp](https://openBF.streamlit.app) to quickly fill in config files.

---

The first line in the configuration file contains the `project_name` variable. 


```yaml
project_name: "my_simulation_name"
```

This is a string that will be used to name the results folder. At the end of the simulation, all the results files will be found in `project_name_results/` directory.

```yaml
inlet_file: "my_inlet.dat"
```

The inlet BC is given trough an ASCII file containing a list of values (pressure or volumetric flow rate) in time. For now, the inlet vessel is the one whose `sn` node is `1`.

```yaml
output_directory; "put/results/in/this/directory/please
```

You can specify where to save results, if not the default directory will be `./<project_name>_results`.

Following the `write_results` list is defined.

```yaml
write_results: ["P"] # ["P", "Q", "u", "A"]
```

Here you can specify which quantities to write in the output files (`P`ressure, `Q` flow, `u` velocity and `A`rea). These are all optional but at least one is recommended.

Then there are three main sections: `solver`, `blood`, and `network`.

## solver

```yaml
solver:
  Ccfl: 0.9
  cycles: 100
  convergence_tolerance: 5.0
  jump: 100
```

It contains the values for the numerical scheme parameters:

- `Ccfl` is the Courant's number used to compute the $\Delta t$ and it is usually taken equal to $0.9$.

- `cycles` is the maximum number of cardiac cycles to be simulated. This is used to stop openBF execution in case of non-converging simulations. An openBF simulation usually takes less than 20 cardiac cycles to converge.

- `convergence_tolerance` is the maximum error in $mmHg$ allowed between two consecutive cardiac cycle to claim convergence.

- `jump` is the number of time-points to be saved in the result files.

## blood

```yaml
blood:
  rho: 1060.0
  mu: 0.004
```

Blood rheological properties:

- `mu` dynamic viscosity in $Pa \cdot s$;
- `rho` density in $kg \cdot m^{-3}$.

## network

This contains the list of vessels in the network. Each vessel has the following __mandatory__ properties:

- `label` the vessel name that will be used to name result files;
- `sn` segment source node;
- `tn` segment target node;
- `L` vessel lenght in $m$;
- `E` wall Young's modulus in $Pa$;
- `R0` or `Rp` and `Rd` describe the lumen radius. If `R0` is defined, the vessel is assumed to have a constant reference lumen radius; it `Rp` and `Rd` are specified, the vessel is set to taper linearly from the proximal lumen radius (i.e. `Rp`, the lumen radius at `sn`) to the distal lumen radius (i.e. `Rd`, the lumen radius at `tn`).

__Optional__ parameters are:

- `to_save` a boolean flag (default `true`) to tell openBF to save results for the current vessel;

- `M` is the number of division along the vessels used to compute the artery $\Delta x$. When not specified, openBF automatically meshes the arteries so that $\Delta x$ is at least $1 mm$;

- `Pext` vessel external pressure in $Pa$, default $0.0 Pa$;

- `gamma_profile` is the radial velocity profile parameter used in the calculation of the viscous losses term, default $2$ (parabolic profile).


## Boundary conditions

The system boundary conditions (BCs) are applied to inlet vessel(s) and outlet vessel(s).

In the case a vessel outlet is not connected to any other vessel, an outlet BC must be assigned by imposing a reflection coefficient `Rt` or by coupling a _windkessel_ model. In case of `Rt`:

- `Rt` 1.0 ≤ Rt ≥ -1.0

In case of two-element windkessel:
- `R1` peripheral resistance
- `Cc` peripheral compliance

In case of three-element windkessel:
- `R1` first peripheral resistance
- `R2` second peripheral resistance
- `Cc` peripheral compliance

You can also enable `inlet_impedance_matching` to let openBF optimise `R1` at runtime. This will match the windkessel inlet impedance and minimise artificial reflections (recommended).

## Template

The configuration file template is reported below. This contains all the mandatory and optional paramenter that can be specified in openBF. The parameter type is reported for each line after the #.

```yml
project_name: <project name> # String
inlet_file: <inlet name>.dat
write_results: ["P"] # ["P", "Q", "u", "A"]

solver:
  Ccfl: <Courant's number> # 0.0 < Ccfl ≤ 1.0; Float
  cycles: <Max number of cardiac cycles> # Int
  convergence_tolerance: <Max RMSE between two cycles> # Pressure, Float
  jump: <Number of timepoints in result files; default 100> # Int

blood:
  mu:  <dynamic viscosity> # [Pa⋅s]; Float
  rho: <density> # [kg/m^3]; Float

network:
  - label: <vessel name> # String
    to_save: true # Bool (default true)
    sn: <source node> # Int
    tn: <target node> # Int

    L: <length> # [m]; Float
    M: <number of divisions along the vessel; default so that Δx = L/M = 1.0mm; minimum M=5> # Int (optional)
    E: <wall Young's modulus> # [Pa]; Float
    h0: <wall thickness; default computed as h0(x) = f(R0)> # [m]; Float (optional)

    R0: <constant lumen radius> # [m]; Float
    #------ OR ------ to assign a linear tapering
    Rp: <proximal lumen radius, i.e. lumen radius at sn>
    Rd: <distal lumen radius, i.e. lumen radius at tn>

    Pext: <external pressure; default 0.0 Pa> # [Pa]; Float (optional)

    gamma_profile: <radial velocity profile parameter; default 9> # ≥ 2.0; Float

    # outlet wk2
    R1: <windkessel peripheral resistance> # [Pa⋅s⋅m^-3]; Float
    Cc: <compliance> # [m^3⋅Pa^-1]; Float
    #------ OR ------ outlet wk3
    R1: <windkessel inlet impedance>
    R2: <peripheral resistance>
    Cc: <compliance>
    inlet_impedance_matching: true # Bool (default false)
    #------ OR ------ outlet reflection
    Rt: <reflection coefficient> # 1.0 ≤ Rt ≥ -1.0; Float
```
