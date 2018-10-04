# Configuration file

The first line in the configuration file contains the `project name` variable. This is a string that will be used to name the results folder. At the end of the simulation, all the results files will be found in `project name_results/` directory.

Following, the configuration file has three main sections: `solver`, `blood`, and `network`.

## solver

It contains the values for the numerical scheme parameters:

- `Ccfl` is the Courant's number used to compute the $\Delta t$ and it is usually taken equal to $0.9$;

```@docs
openBF.calculateDeltaT
```

- `cycles` is the maximum number of cardiac cycles to be simulated. This is used to stop openBF execution in case of non-converging simulations. An openBF simulation usually takes less than 20 cardiac cycles to converge, therefore $100$ cycles is a good threshold value;

- `convergence tolerance` is the maximum percentage error allowed between two consecutive cardiac cycle to claim convergence. During openBF development and validation, this was set to $5$%;

```@docs
openBF.checkConvergence
```

- `jump` is the number of time-points to be saved in the result files; by default this is $100$ and can be omitted in the `.yml` file.

## blood

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
- `R0` or `Rp` and `Rd` describe the lumen radius. If `R0` is defined, the vessel is assumed to have a constant reference lumen radius; it `Rp` and `Rd` are specified, the vessel is set to taper linearly from the proximal lumen radius (i.e. `Rp`, the lumen radius at `sn`) to the distal lumen radius (i.e. `Rd`, the lumen radius at `tn`)

```@docs
openBF.computeRadiusSlope
```

where

- `M` is the number of division along the vessels used to compute the artery $\Delta x$. The value for `M` is __optional__ as by default openBF automatically meshes the arteries so that $\Delta x$ is at least $1 mm$.

```@docs
openBF.meshVessel
```

Other __optional__ parameters are:

- `h0` constant value for wall thickness in $m$. If this is not specified, the thickness is automatically computed w.r.t. the lumen radius

```@docs
openBF.computeThickness
```

- `Pext` vessel external pressure in $Pa$. If this is not specified, the external pressure is set to $0.0 Pa$;

- `gamma profile` is the radial velocity profile parameter used in the calculation of the viscous losses term.

```@docs
openBF.computeViscousTerm
```

## Boundary conditions

The system boundary conditions (BCs) are applied to inlet vessel(s) and outlet vessel(s) only.
The inlet BC is given trough an ASCII file containing a list of values (pressure or volumetric flow rate) in time. To indicate a vessel as an inlet vessel, add the following parameters to its description:

- `inlet` specifies the type of inlet and it can be either `Q` (volumetric flow rate) or `P`;
- `inlet number` is a progressive integer taking count of the number of inlets. The first inlet will have `inlet number: 1`, the second (if any) `inlet number: 2` and so on;
- `inlet file` the path to the inlet time function file in ASCII format. The inlet file values should be arranged in two columns separated by a blank space. The former contains the time values and the latter contains the pressure/flow values in SI units.

In the case a vessel outlet is not connected to any other vessel, an outlet BC must be assigned by imposing a reflection coefficient `Rt` or by coupling a _windkessel_ model. The outlet BC is defined by

- `outlet` specifies the type of inlet: `reflection`, `wk3` (three-element windkessel), or `wk2` (two-element windkessel);

In case of `outlet: reflection` you need to add also
- `Rt` 1.0 ≤ Rt ≥ -1.0

In case of `outlet: wk2`
- `R1` peripheral resistance
- `Cc` peripheral compliance

In case of `outlet: wk3`
- `R1` first peripheral resistance
- `R2` second peripheral resistance
- `Cc` peripheral compliance

The windkessel equation solved is the same for `wk2` and `wk3` cases. In the former case, openBF will automatically assign a value to the second resistance so that the total peripheral resistance will be equal to the one assigned in the configuration file.

```@docs
openBF.computeWindkesselInletImpedance
```

## Template

The configuration file template is reported below. This contains all the mandatory and optional paramenter that can be specified in openBF. The parameter type is reported for each line after the #.

```yml
project name: <project name> # String

solver:
  Ccfl: <Courant's number> # 0.0 < Ccfl ≤ 1.0; Float
  cycles: <Max number of cardiac cycles> # Int
  convergence tolerance: <Max %error between two cycles> # Float
  jump: <Number of timepoints in result files; default 100> # Int

blood:
  mu:  <dynamic viscosity> # [Pa⋅s]; Float
  rho: <density> # [kg/m^3]; Float

network:
  - label: <vessel name> # String
    sn: <source node> # Int
    tn: <target node> # Int

    L: <length> # [m]; Float
    E: <wall Young's modulus> # [Pa]; Float

    R0: <constant lumen radius> # [m]; Float
    #------ OR ------ to assign a linear change of radius along the vessel
    Rp: <proximal lumen radius, i.e. lumen radius at sn>
    Rd: <distal lumen radius, i.e. lumen radius at tn>

    h0: <constant wall thickness; if not specified, computed as h0(R0)> # Float
    M: <number of divisions along the vessel; default M: Δx = L/M = 1.0mm> # Int
    Pext: <external pressure; default 0.0 Pa> # [Pa]; Float

    inlet: <type of inlet> # Q or P; String
    inlet number: <inlet ID> # ≥ 1; Int
    inlet file: <path to inlet time function file> # String

    gamma profile: <radial velocity profile parameter; default 9> # ≥ 2.0; Float

    outlet: <type of outlet> # String
    # outlet: wk2
    R1: <windkessel peripheral resistance> # [Pa⋅s⋅m^-3]; Float
    Cc: <compliance> # [m^3⋅Pa^-1]; Float
    #------ OR ------ outlet: wk3
    R1: <windkessel inlet impedance>
    R2: <peripheral resistance>
    Cc: <compliance>
    #------ OR ------ outlet: reflection
    Rt: <reflection coefficient> # 1.0 ≤ Rt ≥ -1.0; Float
```
