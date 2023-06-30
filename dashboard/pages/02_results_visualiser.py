import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def strip_filename(filename):
    p = Path(filename)
    vessel, q = p.stem.rsplit("_", 1)
    return vessel, q


MULTIPLIERS = {"P": 1 / 133.332, "Q": 1e6, "A": 1, "u": 1, "c": 1}
UNITS = {
    "P": "$(mmHg)$",
    "Q": "$(ml \cdot s^{-1})$",
    "A": "$(m^2)$",
    "u": "$(m \cdot s^{-1})$",
    "c": "$(m \cdot s^{-1})$",
}

uploaded_file = st.file_uploader("Pick a file (.temp, .last or .out)")
if uploaded_file:
    vessel, q = strip_filename(uploaded_file.name)
    df = pd.read_csv(uploaded_file, sep=" ", header=None)

    col_names = ["time", "inlet", "2/5", "middle", "4/5", "outlet"]
    cols = dict(zip(range(7), col_names))
    df = df.rename(columns=cols)

    for c in col_names[1:]:
        df[c] = df[c].apply(lambda x: MULTIPLIERS[q] * x)

    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        inlet = st.checkbox("inlet")
    with c2:
        middle = st.checkbox("middle")
    with c3:
        outlet = st.checkbox("outlet")
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if inlet:
        sns.lineplot(data=df, x="time", y="inlet", ax=ax)
    if middle:
        sns.lineplot(data=df, x="time", y="middle", ax=ax)
    if outlet:
        sns.lineplot(data=df, x="time", y="outlet", ax=ax)

    plt.ylabel(f"{q} {UNITS[q]}")
    plt.title(vessel)
    plt.xlabel("time $(s)$")

    sns.despine()
    st.pyplot(fig)
