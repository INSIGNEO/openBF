# openBF

openBF is an open-source 1D blood flow solver based on MUSCL finite-volume numerical scheme, written in [Julia](https://julialang.org/downloads/) and released under [Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0) free software license.

openBF has been developed by [Dr. A. Melis](https://alemel.is) under the supervision of [Dr. A. Marzo](https://www.sheffield.ac.uk/mecheng/academic-staff/alberto-marzo) as part of the PhD project on cardiovascular modelling at [INSIGNEO](https://www.sheffield.ac.uk/insigneo) Institute, Department of Mechanical Engineering of The University of Sheffield. All the contributors are listed on [github](https://github.com/INSIGNEO/openBF/graphs/contributors).

The accompaigning [PhD Thesis](https://etheses.whiterose.ac.uk/19175/) contains the relevant theoretical background and the explanation of the numerical scheme implemented. A shorter backgound introduction is given in the [Overview](man/overview.md) page.

Head over to the [Quickstart](man/quickstart.md) section for installation and running instructions.

## Publications 

openBF has been used in the following works:

- Benemerito I, Mustafa A, Wang N, Narata AP, Narracott A, Marzo A. [A multiscale computational framework to evaluate flow alterations during mechanical thrombectomy for treatment of ischaemic stroke](https://www.frontiersin.org/articles/10.3389/fcvm.2023.1117449/full), _Frontiers in Cardiovascular Medicine_, 2023. DOI: 10.3389/fcvm.2023.1117449

- Ning W, Sharma K, Sourbron SP, Benemerito I, Marzo A. [Distinguishing hypertensive renal injury from diabetic nephropathy using MR imaging and computational modelling of renal blood flow](https://vph-conference.org/), _VPH2022_ September 2022, Porto, PO _In proceedings_

- Benemerito I, Narata AP, Narracott A, Marzo A. [Determining clinically-viable biomarkers for ischaemic stroke through a mechanistic and machine learning approach](https://link.springer.com/article/10.1007/s10439-022-02956-7), _Annals of Biomedical Engineering_, 2022. DOI: 10.1007/s10439-022-02956-7

- Mustafa A. [An efficient computational approach to guide intervention in treatment of stroke](https://etheses.whiterose.ac.uk/29992/), _PhD Thesis_, 2021

- Benemerito I, Jordan B, Mustafa A, Marzo A. [Quantification of the effects of ageing, hypertension and atherosclerosis on flow reversal during a mechanical thrombectomy procedure](https://www.sheffield.ac.uk/insigneo/overview/events/biomedeng-2021-conference), _BioMedEng21_ September 2021, Sheffield, UK _In proceedings_

- Benemerito I, Narata AP, Narracott A, Marzo A. [Pulsatility indices can inform on distal perfusion following ischaemic stroke](https://www.cmbbe-symposium.com/2021/wp-content/uploads/sites/2/2021/09/Program-CMBBE21-A4.qxp_Detailed.pdf), _CMMBE_ September 2021, Online event, _In proceedings_

- Benemerito I, Narata AP, Narracott A, Marzo A. [Identification of biomarkers for perfusion following an ischaemic event](https://cbmc21.vfairs.com/), _CMBE_ September 2021, Online event, _In proceedings_

- Mercuri M, Wustmann K, von Tengg-Kobligk H, Goksu C, Hose DR, Narracott A. [Subject-specific simulation for non-invasive assessment of aortic coarctation: towards a translational approach](https://www.sciencedirect.com/science/article/pii/S1350453319302449), _Medical Engineering & Physics_, 2020. DOI: 10.1016/j.medengphy.2019.12.003

- Melis A, Moura F, Larrabide I, Janot K, Clayton RH, Narata AP, Marzo A. [Improved biomechanical metrics of cerebral vasospasm identified via sensitivity analysis of a 1D cerebral circulation model](https://www.sciencedirect.com/science/article/pii/S0021929019302830), _Journal of Biomechanics_, 2019. DOI: 10.1016/j.biomech.2019.04.019

- Melis A. [Gaussian process emulators for 1D vascular models](http://etheses.whiterose.ac.uk/19175/), _PhD Thesis_, 2017.

- Melis A, Clayton RH, Marzo A. [Bayesian sensitivity analysis of a 1D vascular model with Gaussian process emulators](http://rdcu.be/AqLm). _International Journal for Numerical Methods in Biomedical Engineering_, 2017. DOI: 10.1002/cnm.2882

- Melis A, Clayton RH, Marzo A. [A more efficient approach to perform sensitivity analyses in 0D/1D cardiovascular models](http://www.compbiomed.net/2015/cmbe-proceedings.htm), _CMBE_ July 2015, Cachan, FR. _In proceedings_

Have you used openBF for your research? Let us know!

## Cite

If you find openBF useful in your work, we kindly request that you cite the following figshare repository

```bibtex
@misc{openBF.jl-2018,
title={openBF: Julia software for 1D blood flow modelling}, 
url={https://figshare.com/articles/openBF_Julia_software_for_1D_blood_flow_modelling/7166183/1}, 
DOI={10.15131/shef.data.7166183}, 
abstractNote={
openBF is an open-source 1D blood flow solver based on MUSCL finite-volume numerical scheme, written in Julia and released under Apache 2.0 free software license.

See https://github.com/INSIGNEO/openBF for the git repository and https://insigneo.github.io/openBF/ for the documentation.
}, 
publisher={figshare}, 
author={Melis, Alessandro}, 
year={2018}, 
month={Oct}}
```

and the accompaining PhD thesis

```bibtex
@phdthesis{wreo19175,
           month = {August},
           title = {Gaussian process emulators for 1D vascular models},
          school = {University of Sheffield},
          author = {Alessandro Melis},
       publisher = {University of Sheffield},
            year = {2017},
             url = {https://etheses.whiterose.ac.uk/19175/},
}
```
