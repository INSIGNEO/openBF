import streamlit as st
from openbf_dashboard.blood import get_blood_props
from openbf_dashboard.solver import get_solver_props
from openbf_dashboard.vessel import get_vessel

def parse_project_name(s: str) -> str:
    if s:
        return s.replace(" ", "_").lower()
    else:
        raise ValueError("project_name must be defined")


st.set_page_config(
    page_title="network_builder",
    layout="wide",
)

c1, c2 = st.columns([1, 1])

with c1:
    y = {}
    y["project_name"] = parse_project_name(
        st.text_input(
            "project_name",
            value="my_simulation",
        )
    )

    with st.expander("Blood"):
        y["blood"] = get_blood_props()

    with st.expander("Solver"):
        y["solver"] = get_solver_props()

    with st.expander("Network"):
        st.text("network definition")
        v = get_vessel()
        if v:
            y["network"] = [v]

with c2:
    st.text("viz")

with c1:
    st.divider()

    if y:
        st.json(y)
