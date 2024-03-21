var documenterSearchIndex = {"docs":
[{"location":"man/overview/","page":"Overview","title":"Overview","text":"The cardiovascular system is a complex network of elastic vessels through which blood is pumped by contraction of the heart. This pulsatile regime and vessel elasticity cause pressure to propagate along the arterial circulation as waves. Mechanical discontinuities caused by bifurcation, bends or cardiovascular pathology cause pressure waves to be reflected in all directions.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"(Image: )","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"Pressure waves measured at a specific location can be seen as the result of a superimposition of incident (forward) and reflected (backward) waves. The analysis of this superposition mechanism allows for the study of mechanical properties upstream and downstream of the measurement point and can be a rich source of diagnostic information about the system through which these waves propagate. However, given the vascular system complexity, it is at present difficult to ascribe a particular waveform feature to a specific trait of the arterial circulation","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"The understanding of the cause-effect mechanisms governing wave propagation in the cardiovascular system was used at INSIGNEO to develop a computer model that can reproduce this behaviour in a quantitative way. openBF is a Julia package for the solution of one-dimensional (1D) blood flow in networks of elastic arteries.","category":"page"},{"location":"man/overview/#Why-1D?","page":"Overview","title":"Why 1D?","text":"","category":"section"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"The 1D model is derived from the three-dimensional formulation of the viscous fluid motion equations, the Navier-Stokes equations. These are simplified by assuming the arteries to be long, narrow, and elastic tubes and the blood to be incompressible (i.e., its density is constant over time) and behaving as a Newtonian fluid (i.e., its viscosity does not depend on its velocity).","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"(Image: )","category":"page"},{"location":"man/overview/#Does-it-work?","page":"Overview","title":"Does it work?","text":"","category":"section"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"Despite all the above assumptions, the resulting partial differential equation system is complex enough to simulate the physics of pulse wave transmission and reflection. This is at the expense of results accuracy in areas where the flow structure is inherently three-dimensional, e.g., in the proximity of heart valves and bifurcations. Conversely, the computational requirements are small compared to 3D simulations, allowing the simulation of the entire cardiovascular system rather than only circumscribed parts.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"(Image: )","category":"page"},{"location":"man/overview/#Acknowledgements","page":"Overview","title":"Acknowledgements","text":"","category":"section"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"The software openBF was developed at INSIGNEO Institute for in silico Medicine at The University of Sheffield (UK) as part of A. Melis’ PhD thesis under the supervision of Dr. A. Marzo. We gratefully acknowledge funding from the UK Engineering and Physical Sciences Research Council (Grant Number EP/K037145/1). openBF is based on MUSCL finite-volume numerical scheme, written in Julia, and released under Apache 2.0 free software license on GitHub.","category":"page"},{"location":"man/overview/","page":"Overview","title":"Overview","text":"The solver is currently used in several HPC systems within the H2020 CompBioMed Centre of Excellence in Computational Biomedicine (Grant Agreement №675451). In SURFSara’s HPC-Cloud, it is being used for large-scale sensitivity analysis and uncertainty quantification studies, towards clinical application. The simulations reported in the thesis and in this paper where run in TUoS’s ShARC tier-3 cluster. More information can be found on the official page and in CompBioMed Software Hub.","category":"page"},{"location":"man/examples/","page":"Examples","title":"Examples","text":"The repository models folder contains configuration and inlet files for:","category":"page"},{"location":"man/examples/","page":"Examples","title":"Examples","text":"In-vitro model from Matthys KS, Alastruey J, Peiró J, Khir AW, Segers P, Verdonck PR, Parker KH, Sherwin SJ. Pulse wave propagation in a model human arterial network: assessment of 1-D numerical simulations against in vitro measurements. Journal of biomechanics. 2007 Dec 31;40(15):3476-86.\nCircle of Willis model from Alastruey J, Parker KH, Peiró J, Byrd SM, Sherwin SJ. Modelling the circle of Willis to assess the effects of anatomical variations and occlusions on cerebral flows. Journal of biomechanics. 2007 Dec 31;40(8):1794-805\nBenchmark models from Boileau E, Nithiarasu P, Blanco PJ, Müller LO, Fossan FE, Hellevik LR, Donders WP, Huberts W, Willemet M, Alastruey J. A benchmark study of numerical schemes for one‐dimensional arterial blood flow modelling. International journal for numerical methods in biomedical engineering. 2015 Oct 1;31(10).","category":"page"},{"location":"#openBF","page":"Home","title":"openBF","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"openBF is an open-source 1D blood flow solver based on MUSCL finite-volume numerical scheme, written in Julia and released under Apache 2.0 free software license.","category":"page"},{"location":"","page":"Home","title":"Home","text":"openBF has been developed by Dr. A. Melis under the supervision of Dr. A. Marzo as part of the PhD project on cardiovascular modelling at INSIGNEO Institute, Department of Mechanical Engineering of The University of Sheffield. All the contributors are listed on github.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The accompaigning PhD Thesis contains the relevant theoretical background and the explanation of the numerical scheme implemented. A shorter backgound introduction is given in the Overview page.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Head over to the Quickstart section for installation and running instructions.","category":"page"},{"location":"#Publications","page":"Home","title":"Publications","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"openBF has been used in the following works:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Ning W, Sharma K, Sourbron SP, Benemerito I, Marzo A. Distinguishing hypertensive renal injury from diabetic nephropathy using MR imaging and computational modelling of renal blood flow, VPH2022 September 2022, Porto, PO In proceedings\nBenemerito I, Narata AP, Narracott A, Marzo A. Determining clinically-viable biomarkers for ischaemic stroke through a mechanistic and machine learning approach. Annals of Biomedical Engineering, 2022. DOI: 10.1007/s10439-022-02956-7\nMustafa A. An efficient computational approach to guide intervention in treatment of stroke. PhD Thesis, 2021\nBenemerito I, Jordan B, Mustafa A, Marzo A. Quantification of the effects of ageing, hypertension and atherosclerosis on flow reversal during a mechanical thrombectomy procedure, BioMedEng21 September 2021, Sheffield, UK In proceedings\nBenemerito I, Narata AP, Narracott A, Marzo A. Pulsatility indices can inform on distal perfusion following ischaemic stroke, CMMBE September 2021, Online event, In proceedings\nBenemerito I, Narata AP, Narracott A, Marzo A. Identification of biomarkers for perfusion following an ischaemic event, CMBE September 2021, Online event, In proceedings\nMelis A. Gaussian process emulators for 1D vascular models. PhD Thesis, 2017.\nMelis A, Clayton RH, Marzo A. Bayesian sensitivity analysis of a 1D vascular model with Gaussian process emulators. International Journal for Numerical Methods in Biomedical Engineering, 2017. DOI: 10.1002/cnm.2882\nMelis A, Clayton RH, Marzo A. A more efficient approach to perform sensitivity analyses in 0D/1D cardiovascular models, CMBE July 2015, Cachan, FR. In proceedings","category":"page"},{"location":"","page":"Home","title":"Home","text":"Have you used openBF for your research? Let us know!","category":"page"},{"location":"#Cite","page":"Home","title":"Cite","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you find openBF useful in your work, we kindly request that you cite the following figshare repository","category":"page"},{"location":"","page":"Home","title":"Home","text":"@misc{openBF.jl-2018,\ntitle={openBF: Julia software for 1D blood flow modelling}, \nurl={https://figshare.com/articles/openBF_Julia_software_for_1D_blood_flow_modelling/7166183/1}, \nDOI={10.15131/shef.data.7166183}, \nabstractNote={\nopenBF is an open-source 1D blood flow solver based on MUSCL finite-volume numerical scheme, written in Julia and released under Apache 2.0 free software license.\n\nSee https://github.com/INSIGNEO/openBF for the git repository and https://insigneo.github.io/openBF/ for the documentation.\n}, \npublisher={figshare}, \nauthor={Melis, Alessandro}, \nyear={2018}, \nmonth={Oct}}","category":"page"},{"location":"","page":"Home","title":"Home","text":"and the accompaining PhD thesis","category":"page"},{"location":"","page":"Home","title":"Home","text":"@phdthesis{wreo19175,\n           month = {August},\n           title = {Gaussian process emulators for 1D vascular models},\n          school = {University of Sheffield},\n          author = {Alessandro Melis},\n       publisher = {University of Sheffield},\n            year = {2017},\n             url = {https://etheses.whiterose.ac.uk/19175/},\n}","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"The first line in the configuration file contains the project_name variable. ","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"project_name: \"my_simulation_name\"","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"This is a string that will be used to name the results folder. At the end of the simulation, all the results files will be found in project_name_results/ directory.","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"inlet_file: \"my_inlet.dat\"","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"The inlet BC is given trough an ASCII file containing a list of values (pressure or volumetric flow rate) in time. For now, the inlet vessel is the one whose sn node is 1.","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"Following the write_results list is defined.","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"write_results: [\"P\"] # [\"P\", \"Q\", \"u\", \"A\"]","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"Here you can specify which quantities to write in the output files (Pressure, Q flow, u velocity and Area). These are all optional but at least one is recommended.","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"Then there are three main sections: solver, blood, and network.","category":"page"},{"location":"man/config/#solver","page":"Configuration","title":"solver","text":"","category":"section"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"solver:\n  Ccfl: 0.9\n  cycles: 100\n  convergence_tolerance: 5.0\n  jump: 100","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"It contains the values for the numerical scheme parameters:","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"Ccfl is the Courant's number used to compute the Delta t and it is usually taken equal to 09.\ncycles is the maximum number of cardiac cycles to be simulated. This is used to stop openBF execution in case of non-converging simulations. An openBF simulation usually takes less than 20 cardiac cycles to converge.\nconvergence_tolerance is the maximum error in mmHg allowed between two consecutive cardiac cycle to claim convergence.\njump is the number of time-points to be saved in the result files.","category":"page"},{"location":"man/config/#blood","page":"Configuration","title":"blood","text":"","category":"section"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"blood:\n  rho: 1060.0\n  mu: 0.004","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"Blood rheological properties:","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"mu dynamic viscosity in Pa cdot s;\nrho density in kg cdot m^-3.","category":"page"},{"location":"man/config/#network","page":"Configuration","title":"network","text":"","category":"section"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"This contains the list of vessels in the network. Each vessel has the following mandatory properties:","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"label the vessel name that will be used to name result files;\nsn segment source node;\ntn segment target node;\nL vessel lenght in m;\nE wall Young's modulus in Pa;\nR0 or Rp and Rd describe the lumen radius. If R0 is defined, the vessel is assumed to have a constant reference lumen radius; it Rp and Rd are specified, the vessel is set to taper linearly from the proximal lumen radius (i.e. Rp, the lumen radius at sn) to the distal lumen radius (i.e. Rd, the lumen radius at tn).","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"Optional parameters are:","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"to_save a boolean flag (default true) to tell openBF to save results for the current vessel;\nM is the number of division along the vessels used to compute the artery Delta x. When not specified, openBF automatically meshes the arteries so that Delta x is at least 1 mm;\nPext vessel external pressure in Pa, default 00 Pa;\ngamma_profile is the radial velocity profile parameter used in the calculation of the viscous losses term, default 2 (parabolic profile).","category":"page"},{"location":"man/config/#Boundary-conditions","page":"Configuration","title":"Boundary conditions","text":"","category":"section"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"The system boundary conditions (BCs) are applied to inlet vessel(s) and outlet vessel(s).","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"In the case a vessel outlet is not connected to any other vessel, an outlet BC must be assigned by imposing a reflection coefficient Rt or by coupling a windkessel model. In case of Rt:","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"Rt 1.0 ≤ Rt ≥ -1.0","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"In case of two-element windkessel:","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"R1 peripheral resistance\nCc peripheral compliance","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"In case of three-element windkessel:","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"R1 first peripheral resistance\nR2 second peripheral resistance\nCc peripheral compliance","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"You can also enable inlet_impedance_matching to let openBF optimise R1 at runtime. This will match the windkessel inlet impedance and minimise artificial reflections (recommended).","category":"page"},{"location":"man/config/#Template","page":"Configuration","title":"Template","text":"","category":"section"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"The configuration file template is reported below. This contains all the mandatory and optional paramenter that can be specified in openBF. The parameter type is reported for each line after the #.","category":"page"},{"location":"man/config/","page":"Configuration","title":"Configuration","text":"project_name: <project name> # String\ninlet_file: <inlet name>.dat\nwrite_results: [\"P\"] # [\"P\", \"Q\", \"u\", \"A\"]\n\nsolver:\n  Ccfl: <Courant's number> # 0.0 < Ccfl ≤ 1.0; Float\n  cycles: <Max number of cardiac cycles> # Int\n  convergence_tolerance: <Max RMSE between two cycles> # Pressure, Float\n  jump: <Number of timepoints in result files; default 100> # Int\n\nblood:\n  mu:  <dynamic viscosity> # [Pa⋅s]; Float\n  rho: <density> # [kg/m^3]; Float\n\nnetwork:\n  - label: <vessel name> # String\n    sn: <source node> # Int\n    tn: <target node> # Int\n\n    L: <length> # [m]; Float\n    M: <number of divisions along the vessel; default so that Δx = L/M = 1.0mm; minimum M=5> # Int (optional)\n    E: <wall Young's modulus> # [Pa]; Float\n    h0: <wall thickness; default computed as h0(x) = f(R0)> # [m]; Float (optional)\n\n    R0: <constant lumen radius> # [m]; Float\n    #------ OR ------ to assign a linear tapering\n    Rp: <proximal lumen radius, i.e. lumen radius at sn>\n    Rd: <distal lumen radius, i.e. lumen radius at tn>\n\n    Pext: <external pressure; default 0.0 Pa> # [Pa]; Float (optional)\n\n    gamma_profile: <radial velocity profile parameter; default 9> # ≥ 2.0; Float\n\n    # outlet wk2\n    R1: <windkessel peripheral resistance> # [Pa⋅s⋅m^-3]; Float\n    Cc: <compliance> # [m^3⋅Pa^-1]; Float\n    #------ OR ------ outlet wk3\n    R1: <windkessel inlet impedance>\n    R2: <peripheral resistance>\n    Cc: <compliance>\n    inlet_impedance_matching: true\n    #------ OR ------ outlet reflection\n    Rt: <reflection coefficient> # 1.0 ≤ Rt ≥ -1.0; Float","category":"page"},{"location":"man/quickstart/#Installation","page":"Quickstart","title":"Installation","text":"","category":"section"},{"location":"man/quickstart/","page":"Quickstart","title":"Quickstart","text":"The latest Julia binary can be downloaded for all platforms from the official website.","category":"page"},{"location":"man/quickstart/","page":"Quickstart","title":"Quickstart","text":"openBF can be installed via Pkg.clone. Start julia in a terminal session and run","category":"page"},{"location":"man/quickstart/","page":"Quickstart","title":"Quickstart","text":"using Pkg\nPkg.clone(\"https://github.com/INSIGNEO/openBF.git\")","category":"page"},{"location":"man/quickstart/#Run-a-simulation","page":"Quickstart","title":"Run a simulation","text":"","category":"section"},{"location":"man/quickstart/","page":"Quickstart","title":"Quickstart","text":"openBF API consists of a single function","category":"page"},{"location":"man/quickstart/","page":"Quickstart","title":"Quickstart","text":"run_simulation(\n        yaml_config_path::String;\n        verbose::Bool = false,\n        out_files::Bool = false,\n        conv_ceil::Bool = false,\n        )","category":"page"},{"location":"man/quickstart/","page":"Quickstart","title":"Quickstart","text":"whose main argument is the name of a .yaml file (see here for an introduction to yaml) containining the description of the vascular network, blood properties, and numerical solver parameters (see Configuration).","category":"page"},{"location":"man/quickstart/","page":"Quickstart","title":"Quickstart","text":"Running a simulation consists in calling","category":"page"},{"location":"man/quickstart/","page":"Quickstart","title":"Quickstart","text":"using openBF\nrun_simulation(\"input.yml\", verbose=true)","category":"page"},{"location":"man/quickstart/","page":"Quickstart","title":"Quickstart","text":"This will create a <project_name>_results folder containing all the output files from the simulation.","category":"page"}]
}
