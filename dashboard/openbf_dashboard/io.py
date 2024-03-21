import streamlit as st

def parse_project_name(s: str) -> str:
    if s:
        return s.replace(" ", "_").lower()
    else:
        raise ValueError("project_name must be defined")

def get_io():
    i = {}
    i["project_name"] = parse_project_name(
        st.text_input(
            "project_name",
            value="my_simulation",
        )
    )

    i["inlet_file"] = st.text_input(
                "inlet_file",
                help="path to project_name_inlet.dat file",
                value="inlet.dat",
                key=f"inlet_file",
            )

    i["output_directory"] = st.text_input(
                "output_directory",
                help="where to save the output files",
                value=f"{i['project_name']}_results",
                key=f"output_directory",
            )

    i["write_results"] = []
    st.text("write_results")
    check_cols = st.columns(4)
    for j, q in enumerate(["P", "Q", "A", "u"]):
        with check_cols[j]:
            if st.checkbox(q, key=f"save_{q}", value=q in ["P", "Q"]):
                i["write_results"].append(q)

    return i
