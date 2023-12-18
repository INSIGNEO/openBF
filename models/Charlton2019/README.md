# Reproducing simulations from Charlton et al. [1]
Code in this folder enables to reproduce the whole-body hemodynamics simulations from [1].
It produces configuration files for openBF from the ones used by Charlton et al.

## Usage:
A configuration file and inlet flow have been pre-generated, so the following instructions are only necessary if you aim to generate simulations for patient different from the reference one.
In that case, you need Python (>=3.0) with `pandas`, `numpy`, `matplotlib`, `argparse` and, `yaml` installed. Then `cd` into the `whole-body` folder and execute 
```
python create_template.py --patient_id $ID_OF_PATIENT
```
where `$ID_OF_PATIENT` is a value between 1 and 4374 that corresponds to one simulation from [1].

In addition to the standard usage, the code can be slightly adapted to generate simulations with new heart function parameters and different patient's height by calling the function `create_template_and_config` with the appropriate arguments (see documentation in the code).

## Limitations:
The code and openBF in their current implementations do not reproduce exactly the simulations of [1], as can be found [here](https://peterhcharlton.github.io/pwdb/).

## References:
[1]: Charlton, P. H., Mariscal Harana, J., Vennin, S., Li, Y., Chowienczyk, P., & Alastruey, J. (2019). Modeling arterial pulse waves in healthy aging: a database for in silico evaluation of hemodynamics and pulse wave indexes. American Journal of Physiology-Heart and Circulatory Physiology, 317(5), H1062-H1085.