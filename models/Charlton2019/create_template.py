import argparse

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from yaml.loader import SafeLoader
import math
import pickle5 as pickle
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use('agg')

class Namespace():
    def __init__(self):
        return

# Code adapted by Antoine Wehenkel from the Matlab code of Charlton's dataset paper
# https://github.com/peterhcharlton/pwdb/blob/master/pwdb_v0.1/Additional%20Functions/AorticFlowWave.m
# Charlton et al: Charlton, P. H., Mariscal Harana, J., Vennin, S., Li, Y., Chowienczyk, P., & Alastruey, J. (2019). Modeling arterial pulse waves in healthy aging: a database for in silico evaluation of hemodynamics and pulse wave indexes. American Journal of Physiology-Heart and Circulatory Physiology, 317(5), H1062-H1085.
def create_waveform(HR: float, SV: float, LVET: float, PFT: float, RFV: float, **kwargs):
    '''
    Return: the inlet flow that corresponds to the given parameters (see Charlton et al. for details).
    Parameters:
        HR: The heart rate in beats per minute
        SV: The stroke volume in mL
        LVET: The left ventricular ejection time in ms
        PFT: The peak flow time in ms
        RFV: The reflected fraction of volume in mL
        kwargs: Only used to make the function resilient to being called with a dictionary with other parameters.
    '''
    with open('./config/Q_template.pickle', 'rb') as handle:
        Q = pickle.load(handle)
    Q = np.array(Q)
    Q[1:52] = 0
    Q[362:] = 0
    Q = np.append(Q[53:], Q[1:52]) / 10 ** 6  # from mL to m^3
    template = Q
    t = np.linspace(0, (Q.shape[0] - 1) / 1000, Q.shape[0])
    assert (len(template.shape) == 1)

    CO = HR * SV / 1000  # in l/min
    rev_flow = 1  # between 0 and 1 - the proportion of reverse flow to include
    contractility_const = 1  # a multiplication factor applied to the time of the systolic peak
    T = 60. / HR
    T_Peak_Flow = PFT  # (0.080/0.282)*LVET  # Mynard's
    T_sys_upslope_ip = (0.010 / 0.080) * T_Peak_Flow;  # Mynard's
    Q_sys_upslope_ip = 0.32  # Haven't measured Mynard's - don't think it's used
    T_Min_Flow = (0.298 / 0.282) * LVET  # Mynard's
    T_dia = (0.309 / 0.282) * LVET  # Mynard's
    Reg_Vol = RFV  # 1.2775 # Mynard's, ml

    # Proportionality constants
    prop = 1.2
    prop2 = 1.1
    prop3 = 1.1
    prop4 = 1.04
    prop5 = 1.05
    prop6 = 1.03

    # timings
    template_flow = Namespace()
    template_flow.T = 1
    template_flow.ti_0 = (1 / prop) * 0.024 * template_flow.T;  # Time of inflection point in systolic increase
    template_flow.tmax_0 = 0.08 * template_flow.T;  # Time of systolic peak
    template_flow.ti2_0 = (
                                      1 / prop5) * prop2 * 0.17 * template_flow.T;  # 0.150; Time of inflection point in systolic decrease [old only]
    template_flow.ti3_0 = prop4 * prop2 * 0.21 * template_flow.T;  # 0.150; Time of inflection point in systolic decrease [old only]
    template_flow.ts_0 = 0.290 * template_flow.T;  # Time of start of dicrotic notch
    template_flow.tmin_0 = 0.310 * template_flow.T;  # Time of minimum flow during dicrotic notch  ### THIS HAS BEEN ADJUSTED from 0.30
    template_flow.td_0 = 0.330 * template_flow.T;  # Time of start of diastole

    # flows
    template_flow.Qmax = 1;  # Systolic peak flow()
    template_flow.Qi = 0.367;  # Flow at inflection point on systolic upslope
    template_flow.Qi2 = prop5 * prop3 * (1 / prop2) * 0.647;  # Flow at inflection point in systolic decrease [old]
    template_flow.Qi3 = prop6 * (1 / prop3) * (
                1 / prop2) * 0.496;  # Flow at second inflection point in systolic decrease [old]

    # find approximate Qmin [during reverse flow]
    template_flow.Qmin = -0.1;  # Minimum flow during dicrotic notch

    # round timings to  right sizes of time vector
    # TODO check effect of rounding at 1e-3
    ti = template_flow.ti_0
    tmax = template_flow.tmax_0
    ti2 = template_flow.ti2_0
    ti3 = template_flow.ti3_0
    ts = template_flow.ts_0
    tmin = template_flow.tmin_0
    td = template_flow.td_0

    # time step
    dt = 1 / 1000

    # Get template
    template_flow.age = "young"
    template_flow.v = template
    template_flow.t = np.linspace(0, (template.shape[0] - 1) / 1000., template.shape[0])

    # Now creating a new waveform with desired properties:

    # setup
    fs = 1000  # Hz
    dt = 1 / fs  # seconds

    # Find sytolic upslope
    deb = 0
    fin = np.argmax(template_flow.v)
    sys_upslope = Namespace()
    sys_upslope.v = template_flow.v[deb:fin]
    sys_upslope.t = np.linspace(0, T_Peak_Flow, sys_upslope.v.shape[0])

    # Find systolic downslope
    deb = fin + 1
    fin = np.argmax(template_flow.v <= 0)  # Find first negative value.
    # fin = find(template_flow.v[1:end-1] .> 0 & template_flow.v[2:end] <= 0); fin = fin[1]
    sys_downslope = Namespace()
    sys_downslope.v = template_flow.v[deb:fin]
    sys_downslope.t = np.linspace(sys_upslope.t[-1] + dt, LVET, sys_downslope.v.shape[0])

    # Find reverse flow()
    deb = fin - 1
    fin = np.argwhere(template_flow.v < 0.)[-1][0]  # last negative value
    # fin = find(template_flow.v != 0, 1, "last")
    reverse = Namespace()
    reverse.v = template_flow.v[deb:fin]
    reverse.t = np.linspace(sys_downslope.t[-1] + dt, T_dia, reverse.v.shape[0])

    # Find diastolic flat line()
    no_els = T - T_dia
    diastolic = Namespace()
    diastolic.t = np.arange(T_dia + dt, T * fs + dt, 1)
    diastolic.v = np.zeros(diastolic.t.shape[0])

    # concatenate to give waveform
    mod_flow = Namespace()
    mod_flow.t = np.concatenate([sys_upslope.t, sys_downslope.t, reverse.t, diastolic.t])
    mod_flow.v = np.concatenate([sys_upslope.v, sys_downslope.v, reverse.v, diastolic.v])

    # resample to give constant fs; without irregular spacing at joins
    # inflow.t = (0 : floor(mod_flow.t[end]*inflow.input_params.fs))/inflow.input_params.fs
    # inflow.v = interp1(mod_flow.t, mod_flow.v, inflow.t)
    inflow = Namespace()
    inflow.t = mod_flow.t
    inflow.v = mod_flow.v

    # Scale to give desired reverse flow volume
    scale_factor = Reg_Vol / (inflow.v[inflow.v < 0] * dt).sum() / 1e6
    inflow.v[inflow.v < 0] = -inflow.v[inflow.v < 0] * scale_factor;  # flow is now in ml/s

    # Scale to give desired stroke volume
    # AW: I do not understand why the RFV is added to the SV to compute the scale_factor
    # scale_factor = (inflow.input_params.SV+inflow.input_params.Reg_Vol)/(sum(inflow.v(inflow.v>0))*dt);
    scale_factor = SV / (inflow.v[inflow.v > 0] * dt).sum() / 1e6
    inflow.v[inflow.v > 0] = inflow.v[inflow.v > 0] * scale_factor

    return inflow.t / fs, inflow.v


def trim_zeros(t: np.array, Q: np.array):
    '''
    Return the cleaned version of t and Q where more than 2 consecutive 0 are removed (as these make openBF fail to run simulations).
    Parameters:
        t: the time of measurements
        Q: the values of the measurements.
    '''
    start_zero = -1
    new_t = []
    new_Q = []
    for i in range(len(t)):
        if Q[i] == 0 and i < (len(t) - 1):
            if start_zero == -1:
                start_zero = i
        elif Q[i] == 0:
            if start_zero != -1:
                if start_zero < (i - 1):
                    new_t += [t[start_zero], t[i]]
                    new_Q += [1e-9, 1e-10]
                else:
                    new_t += [t[start_zero]]
                    new_Q += [1e-10]
            else:
                new_t.append(t[i])
                new_Q.append(1e-10)
        else:
            if start_zero != -1:
                if start_zero < (i - 1):
                    new_t += [t[start_zero], t[i - 1]]
                    new_Q += [1e-9, 1e-10]
                else:
                    new_t += [t[start_zero]]
                    new_Q += [1e-10]
                start_zero = -1
            new_t.append(t[i])
            new_Q.append(Q[i])

    return np.array(new_t), np.array(new_Q)

def generate_template(params: dict, savefolder: str, plot: bool=True):
    '''
    Return: None. It generates, save and plot the inlet flow.
    Parameters:
        params: the parameters to be given to the `create_waveform` function (HR, SV, LVET, PFT, RFV)
        savefolder: the folder where the plot and template (in the format expected by openBF) will be saved.
        plot: save a plot of the generated inlet flow.
    '''
    t, Q = create_waveform(**params)
    t, Q = trim_zeros(t, Q)
    pd.DataFrame(data=[t, Q]).T.to_csv(savefolder + '/inlet.dat', sep=' ', index=False, header=False)
    plt.figure()
    plt.title(str(params), fontsize=8)
    plt.ylabel("Blood flow [mL/s]")
    plt.xlabel("Time [ms]")
    plt.ylim(-.3, .12*10)
    plt.plot(t*1000, Q*1000) 
    plt.savefig(savefolder + '/template.png')
    plt.close()

def generate_config(save_path, param_heart, height, subject_number=1, tapering=False, split_length=.1, **kwargs):
    create_template_and_config(subject_number, param_heart, tapering, split_length, False, height, save_path)

def create_template_and_config(subject_number: int, param_heart: dict, tapering: bool=False, split_length: float=-1., 
                               viscosity: bool=False, height_patient: float=170., save_path: str=None):
    '''
    Return: None. Create the inlet flow and configuration file to produce a whole-body hemodynamics similar to Charlton et al.
    Parameters:
        subject_number: integer between 1 and 4374 that defines the configuration from Charlton et al. to be used.
        param_heart: If None, will use the heart function parameters from Charlton et al. 
                    Else a dict with HR, SV, LVET, PFT, RFV that defines the inlet flow.
        tapering: Recommend False as results with tapering are not consistent. 
                    If "discrete" will use tapered arteries by splitting them into fixed `split_length` segment of the corresponding radius.
                    If "openBF" will use Rp and Rd of openBF to define linear tappering.
                    if "both" will use Rp and Rd on segmented arteries.
        split_length: Length of arterial segment if tapering is "discrete" or "both"
        viscosity: If True arteries are visco-elastic, not well handled.
        height_patient: If value is different from 170 will rescale the length of the arteries by height_patient/170.
        save_path: Path where config file and inlet flow will be saved.
    '''
    print("Creating parameters from template of patient with ID %d." % subject_number)
    if tapering:
        print("Using tapering.")
    else:
        print("Not using tapering.")
        
    alastruey_config = pd.read_csv('./config/pwdb_model_configs.csv', index_col='Subject Number').loc[subject_number, :]

    if param_heart is None:
        param_heart = {'HR': float(alastruey_config[' hr [bpm]']),
                       'LVET': float(alastruey_config[' lvet [ms]']),
                       'SV': float(alastruey_config[' sv [ml]']),
                       'PFT': float(alastruey_config[' pft [ms]']),
                       'RFV': float(alastruey_config[' rfv [ml]'])
                        }
        print(param_heart)
        print("Using the heart function parameters from Charlton et al for patient %d." % subject_number)
    
    # Generate the inlet flow as a function of Heart Rate (HR), Left Ventricular Ejection Time (LVET), Stroke Volume (SV), Peak Flow Time (PFT), Reflected Fraction Volume (RFV).
    generate_template(param_heart, save_path)

    rho = 1060

    # Will increase length of arteries linearly with the height of the simulated patient.
    height_ratio = height_patient/170.
    

    # Now taking parameters from the config files in Charlton's paper and making them compatible with openBF:
    # Note: Elasticity, viscosity and tapering are not handled very well by openBF and the conversion here may be incorrect.
    # This code has shown to give realistic simulations when tapering is ignored or a
    k1 = alastruey_config[" k1 [g/s^2/cm]"]#3e6 # unit: g*s^-2*cm^-1
    k2 = alastruey_config[" k2 [/cm]"]#-13.5
    k3 = alastruey_config[" k3 [g/s^2/cm]"]#5.3849e+05#430.118 - 1871.3 * age + 244.11 * age**2
    Eh = lambda R_d: 100 * R_d * (k1 * np.exp(k2 * 100 * R_d) + k3)/1000 #Unit: 10^(-2) * g*s^-2 = 10^(-3) * 0.1Pa = 1N/cm^2

    ah = 0.2802
    bh = -5.053e2
    ch = 0.1324
    dh = -0.1114e2 
    h0 = lambda R_d: R_d *(ah * np.exp(bh * R_d) + ch * np.exp(dh * R_d)) #Unit: m

    E0 = lambda Eh, h0: Eh/h0 #Pa*m^-1

    b0 = alastruey_config[' b0 [g/s]']#600 # g*s^-1
    b1 = alastruey_config[' b1 [g cm/s]']#150 # g*cm*s^-1
    gamma = lambda Rd: (b1/(2*100*Rd) + b0)/(1000) # Wall viscosity
    #gamma = lambda Rd: (b1/(2*Rd*100) + b0)/(10 * 100**3) * 1 / (math.pi * Rd**2)
    phi = lambda Rd: 1.5 * gamma(Rd)/(h0(Rd) * math.sqrt(math.pi)) * (math.pi * Rd**2)#/1000.
    phi_bis = lambda gamma, Rd: 1.5 * gamma/(h0(Rd) * math.sqrt(math.pi)) / (math.pi * Rd**2)#/1000.
    
    alastruey_cleaned = pd.read_csv('./config/geo/pwdb_geo_{:04d}.csv'.format(subject_number))

    columns_julia = ['E', 'L', 'R0', 'inlet', 'inlet file', 'inlet number', 'label', 'sn', 'tn', 'Cc', 'R1', 'outlet']
    df_julia_network_as_alastruey = pd.DataFrame(columns=columns_julia)

    df_julia_network_as_alastruey[['L', 'Rp', 'Rd', 'sn', 'tn', 'seg_id', 'Cc', 'R1']] = alastruey_cleaned.loc[:, [' length', ' inlet_radius', ' outlet_radius', ' inlet_node', 
                                                                                                               ' outlet_node', 'seg_no', ' peripheral_c', ' peripheral_r']].values
    alastruey_network = pd.read_csv("./config/116_artery_model.txt", sep='\t')
    df_julia_network_as_alastruey.loc[:, ['label']] = alastruey_network.loc[:, 'Name']
    df_julia_network_as_alastruey.loc[df_julia_network_as_alastruey.label.str.contains('Right Superior Thyroid'), 'label'] = df_julia_network_as_alastruey.loc[df_julia_network_as_alastruey.label.str.contains('Right Superior Thyroid'), 'label'].values[0].replace('/', '').replace(' ', '_')
    df_julia_network_as_alastruey.loc[df_julia_network_as_alastruey.label.str.contains('Left Superior Thyroid'), 'label'] = df_julia_network_as_alastruey.loc[df_julia_network_as_alastruey.label.str.contains('Left Superior Thyroid'), 'label'].values[0].replace('/', '').replace(' ', '_')

    df_julia_network_as_alastruey.loc[0, 'inlet'] = 'Q'
    df_julia_network_as_alastruey.loc[0, 'inlet number'] = 1.0
    df_julia_network_as_alastruey.loc[0, 'inlet file'] = 'inlet.dat'

    # Only long arteries are visco-elastic:
    elmt_shortL_vw = 0.011
    df_julia_network_as_alastruey.loc[:, 'viscous'] = False
    df_julia_network_as_alastruey.loc[:, 'L'] = df_julia_network_as_alastruey.loc[:, 'L'] * height_ratio
    df_julia_network_as_alastruey.loc[df_julia_network_as_alastruey.loc[:, 'L'] > elmt_shortL_vw, 'viscous'] = True
    viscous_arteries = df_julia_network_as_alastruey.loc[:, 'L'] > elmt_shortL_vw
    R0 = 0.5*(df_julia_network_as_alastruey.loc[viscous_arteries, 'Rp'].values + df_julia_network_as_alastruey.loc[viscous_arteries, 'Rd'].values)
    df_julia_network_as_alastruey.loc[viscous_arteries, 'gamma_visco'] = gamma(R0)
    if tapering in ["discrete", "both"]:
        df_julia_network_as_alastruey = split_network(df_julia_network_as_alastruey, split_length)
    else:
        df_julia_network_as_alastruey.loc[:, ['label']] = df_julia_network_as_alastruey.loc[:, ['label']] + ' seg_0'

    factor_peripheral = 1.
    is_terminal = (df_julia_network_as_alastruey['R1'] > 0) | (~df_julia_network_as_alastruey['tn'].isin(df_julia_network_as_alastruey['sn'].values))
    R1 = lambda Rend: 1/(math.pi * Rend**2) * (2/3 * rho * Eh(Rend)/Rend)**.5 #blood.rho*waveSpeed(A0[end], gamma[end])/A0[end]

    df_julia_network_as_alastruey.loc[is_terminal, 'R2'] = df_julia_network_as_alastruey.loc[is_terminal, 'R1'].values  - R1(.5*(df_julia_network_as_alastruey.loc[is_terminal, 'Rp'].values + df_julia_network_as_alastruey.loc[is_terminal, 'Rd'].values))
    df_julia_network_as_alastruey.loc[is_terminal, 'R1'] = R1(.5*(df_julia_network_as_alastruey.loc[is_terminal, 'Rp'].values + df_julia_network_as_alastruey.loc[is_terminal, 'Rd'].values))#df_julia_network_as_alastruey.loc[is_terminal, 'R1'].values * factor_peripheral
    df_julia_network_as_alastruey.loc[is_terminal, 'Cc'] = df_julia_network_as_alastruey.loc[is_terminal, 'Cc'].values / factor_peripheral

    is_terminal[0] = False
    df_julia_network_as_alastruey.loc[is_terminal, 'outlet'] = 'wk3'
    df_julia_network_as_alastruey.loc[~is_terminal, 'R1'] = float('nan')
    df_julia_network_as_alastruey.loc[~is_terminal, 'R2'] = float('nan')
    df_julia_network_as_alastruey.loc[~is_terminal, 'Cc'] = float('nan')
    df_julia_network_as_alastruey.loc[is_terminal, 'Pout'] = 33.2 * 133.322 # Pa = N/m^2

    df_julia_network_as_alastruey.loc[:, 'Pext'] = alastruey_config[" dbp [mmHg]"] * 133.322 # Pa = N/m^2

    R0 = (df_julia_network_as_alastruey.loc[:, 'Rp'].values + df_julia_network_as_alastruey.loc[:, 'Rd'].values)/2.
    df_julia_network_as_alastruey.loc[:, 'h0'] = h0(R0)

    alpha_velocity_profile = alastruey_config[" alpha [-]"]
    df_julia_network_as_alastruey.loc[:, 'gamma_profile'] =  (2 - alpha_velocity_profile)/(alpha_velocity_profile - 1.) #R0 * 1000
    df_julia_network_as_alastruey.loc[df_julia_network_as_alastruey.loc[:, 'gamma_profile'] > 9, 'gamma_profile'] = 9.
    df_julia_network_as_alastruey.loc[df_julia_network_as_alastruey.loc[:, 'gamma_profile'] < 2, 'gamma_profile'] = 2.

    df_julia_network_as_alastruey.loc[:, 'E'] = E0(Eh(R0), h0(R0)) # Pa
    df_julia_network_as_alastruey.loc[:, 'phi'] =  0.
    viscous_arteries = df_julia_network_as_alastruey.loc[:, 'viscous'] == True
    R0_viscous = (df_julia_network_as_alastruey.loc[viscous_arteries, 'Rp'].values + df_julia_network_as_alastruey.loc[viscous_arteries, 'Rd'].values)/2.
    gamma_viscous = df_julia_network_as_alastruey.loc[viscous_arteries, 'gamma_visco'].values 
    df_julia_network_as_alastruey.loc[viscous_arteries, 'phi'] = phi(R0_viscous) if viscosity else 0.

    idx_no_lumen_radius_change = (df_julia_network_as_alastruey['Rp'] == df_julia_network_as_alastruey['Rd']) #| (df_julia_network_as_alastruey['L'] < .2) #| (df_julia_network_as_alastruey['Rp'] < .5)
    df_julia_network_as_alastruey.loc[idx_no_lumen_radius_change, 'R0'] = df_julia_network_as_alastruey['Rd']
    df_julia_network_as_alastruey.loc[idx_no_lumen_radius_change, 'Rp'] = float('nan')
    df_julia_network_as_alastruey.loc[idx_no_lumen_radius_change, 'Rd'] = float('nan')

    no_lumen_tappering = tapering not in ["openBF", "both"]
    if no_lumen_tappering:
        df_julia_network_as_alastruey.loc[~idx_no_lumen_radius_change, ['R0']] = (df_julia_network_as_alastruey.loc[~idx_no_lumen_radius_change, ['Rp']].values + df_julia_network_as_alastruey.loc[~idx_no_lumen_radius_change, ['Rd']].values)/2.
        df_julia_network_as_alastruey.loc[~idx_no_lumen_radius_change, 'Rd'] = float('nan')
        df_julia_network_as_alastruey.loc[~idx_no_lumen_radius_change, 'Rp'] = float('nan')

    def remove_nan(d):
        new_dic = {}
        for k, v in d.items():
            if not pd.isna(v):
                new_dic[k] = v
        return new_dic
    df_julia_network_as_alastruey.drop(columns=['seg_id', 'gamma_profile', 'gamma_visco', 'viscous'], inplace=True)
    df_julia_network_as_alastruey[['sn', 'tn']] = df_julia_network_as_alastruey[['sn', 'tn']].astype('int')
    new_network = [remove_nan(l) for l in df_julia_network_as_alastruey.to_dict('records')]

    data = {}
    data['network'] = new_network
    data['Heart function'] = param_heart
    data['blood'] = {}
    data['blood']['mu'] = 0.0025  # in poise = 0.1 x Pa
    data['blood']['rho'] = rho  # in poise = 0.1 x Pa
    data['project name'] = 'julia_patient_{:04d}'.format(1)
    data['solver'] = {}
    data['solver']['cycles'] = 50  # The maximum number of iterations.
    data['solver']['jump'] = 200  # The number of timesteps in the output signals.
    data['solver'][
        'convergence tolerance'] = 1.  # This value can be reduced to achieve more accurate but slower simulation.
    data['solver']['Ccfl'] = .75  # Courant-Friedrichs-Lewy (CFL) number.
    # Open the file and load the file
    save_folder = save_folder if save_path is None else save_path
    with open(save_folder + '/config.yml'.format(subject_number), 'w') as f:
        yaml.dump(data, f)


def split_artery(artery, seg_leng, offset=0):
    nb_segments = round(artery.L / seg_leng) if artery.Rp != artery.Rd and seg_leng < artery.L else 1
    #nb_segments = round(artery.L / seg_leng) if seg_leng < artery.L else 1
    dl = artery.L/nb_segments
    dr = (artery.Rd - artery.Rp)/nb_segments
    sub_arteries = [{'inlet': float('NaN'),
                     'inlet file': float('NaN'),
                     'inlet number': float('NaN'),
                     #'viscous': artery.viscous,
                     #'E': artery.E,
                     #'phi': artery.phi,
                     #'gamma_visco': artery.gamma_visco,
                     'label': artery.label + ' seg_%d' % i,
                     'L': dl, 
                     'Rp': artery.Rp + dr * i, 
                     'Rd': artery.Rp + dr * (i + 1), 
                     'sn': i + offset, 
                     'tn': (i + 1) + offset, 
                     'Cc': float('NaN'), 
                     'R1': float('NaN')} 
                    for i in range(nb_segments)]
    df = pd.DataFrame(data=sub_arteries)
    df.loc[0, 'sn'] = artery.sn
    df.loc[nb_segments-1, 'tn'] = artery.tn
    df.loc[nb_segments-1, 'Cc'] = artery.Cc
    df.loc[nb_segments-1, 'R1'] = artery.R1
    df.loc[0, 'inlet'] = artery['inlet']
    df.loc[0, 'inlet file'] = artery['inlet file']
    df.loc[0, 'inlet number'] = artery['inlet number']
    return df.drop(columns=['seg_id', 'gamma_profile', 'gamma_visco', 'viscous'])


def split_network(df, seg_len):
    new_df = pd.DataFrame(columns=df.columns)
    offset = df.shape[0] - 1
    for i, row in df.iterrows():
        new_artery = split_artery(row, seg_leng=seg_len, offset=offset)
        offset += new_artery.shape[0] - 1
        new_df = pd.concat([new_df, new_artery], axis=0, ignore_index=True)
    return new_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Template of Alastruey\'s dataset for the julia solver.')
    parser.add_argument('--patient_id', type=int, required=True)

    parser.add_argument('--tapering', choices=["False", "openBF", "discrete", "both"], default="False")
    parser.add_argument('--split_length', type=float, default=0.02)
    parser.add_argument('--viscosity', action="store_true")
    args = parser.parse_args()
    create_template_and_config(args.patient_id, tapering=args.tapering, split_length=args.split_length,
                               viscosity=args.viscosity, save_path='.', param_heart=None, height_patient=170)
    print("Configuration file successfully created and saved.")