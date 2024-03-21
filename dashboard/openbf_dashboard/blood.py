import streamlit as st


def get_blood_props():
    c3, c4 = st.columns([1, 1])
    y = {}
    with c3:
        y["rho"] = st.number_input(
            "$\\rho$ $(kg \cdot m^{-3})$",
            help="density",
            value=1060.0,
            step=10.0,
            min_value=1030.0,
            max_value=1090.0,
        )

    with c4:
        y["mu"] = st.number_input(
            "$\\mu$ $(Pa \cdot s)$",
            help="dynamic viscosity",
            value=0.004,
            step=0.001,
            min_value=0.001,
            max_value=0.009,
            format="%e",
        )
    return y
