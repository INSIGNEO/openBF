import streamlit as st


def get_solver_props():
    y = {}
    c1, c2 = st.columns([1, 1])
    with c1:
        y["Ccfl"] = st.number_input(
            "CFL",
            value=0.9,
            min_value=0.01,
            max_value=1.0,
            step=0.1,
            help="Courant's number",
        )
        y["cycles"] = st.number_input(
            "cycles",
            value=100,
            min_value=1,
            step=10,
            help="simulation length in cardiac cycles",
        )
    with c2:
        y["jump"] = st.number_input(
            "jumps",
            value=100,
            min_value=10,
            step=10,
            help="number of time steps in the output files",
        )
        y["convergence tolerance"] = st.number_input(
            "convergence tolerance (%)",
            value=1.0,
            min_value=0.1,
            max_value=100.0,
            help="stop simulation when error < tolerance",
        )

    return y
