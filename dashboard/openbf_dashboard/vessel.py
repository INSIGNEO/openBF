import streamlit as st


def get_vessel():
    v = {}
    v["label"] = st.text_input(
        "label",
        value="vessel_name",
    )

    c1, c2 = st.columns([1, 1])
    with c1:
        sn = st.number_input(
            label="sn",
            min_value=1,
            step=1,
            help="source node",
        )
        L = st.number_input(
            "$\ell$ $(m)$",
            min_value=0.001,
            max_value=1.0,
            step=0.01,
            value=0.01,
            format="%6.5f",
            help="length",
        )
        Rp = st.number_input(
            "$R_{0p}$ $(m)$",
            min_value=0.0001,
            max_value=1.0,
            step=0.0001,
            value=0.001,
            format="%7.6f",
            help="proximal radius",
        )
        E = st.number_input(
            "E $(Pa)$",
            min_value=1e3,
            value=125e3,
            help="wall Young's modulus",
            format="%e",
        )
    with c2:
        tn = st.number_input(
            label="tn",
            min_value=2,
            step=1,
            help="target node",
        )
        M = st.number_input(
            label="M",
            min_value=5,
            value=int(max(5, L * 1e3)),
            help="number of nodes along the vessel",
        )
        Rd = st.number_input(
            "$R_{0d}$ $(m)$",
            min_value=0.0001,
            max_value=1.0,
            step=0.0001,
            value=0.001,
            format="%7.6f",
            help="distal radius",
        )
        gamma_profile = st.slider(
            "$\gamma_v$",
            min_value=2,
            max_value=9,
            help="velocity profile gamma value",
        )

    Pext = st.number_input(
        "$P_{ext}$ $(Pa)$",
        value=0.0,
        step=1e3,
        format="%e",
        help="external pressure",
    )

    if sn == tn:
        raise ValueError("source node must be different than target node")

    v["sn"] = sn
    v["tn"] = tn
    v["L"] = L
    v["M"] = M
    v["Rp"] = Rp
    v["Rd"] = Rd
    v["E"] = E
    v["gamma profile"] = gamma_profile
    v["Pext"] = Pext

    if st.checkbox("inlet"):
        v["inlet"] = st.radio(
            "",
            options=["Q", "P"],
            index=0,
            horizontal=True,
            label_visibility="collapsed",
        )
        v["inlet file"] = st.text_input(
            "inlet file", help="path to project_name_inlet.dat file", value="inlet.dat"
        )

        # WARNING what's this???
        v["inlet number"] = st.number_input(
            "inlet nuber",
            value=1,
            min_value=1,
        )

    if st.checkbox("outlet"):
        v["outlet"] = st.radio(
            "",
            options=["Rt", "wk2", "wk3"],
            index=2,
            horizontal=True,
            label_visibility="collapsed",
        )

        if v["outlet"] == "Rt":
            v["Rt"] = st.number_input(
                "$R_t$",
                min_value=0.0,
                max_value=1.0,
                step=0.1,
                help="reflection coefficient",
            )
        elif v["outlet"] in ["wk2", "wk3"]:
            c3, c4 = st.columns([1, 1])
            with c3:
                v["R1"] = st.number_input(
                    "$R_1$ $(Pa \cdot s \cdot m^{-3}$)",
                    min_value=0.0,
                    value=1e9,
                    format="%e",
                    help="proximal resistance",
                )
                if v["outlet"] == "wk3":
                    v["R2"] = st.number_input(
                        "$R_2$ $(Pa \cdot s \cdot m^{-3}$)",
                        min_value=0.0,
                        value=1e6,
                        format="%e",
                        help="distal resistance",
                    )
            with c4:
                v["Cc"] = st.number_input(
                    "$C_c$ $(m^3 \cdot Pa$)",
                    min_value=0.0,
                    value=1e-10,
                    format="%e",
                    help="compliance",
                )

    return v
