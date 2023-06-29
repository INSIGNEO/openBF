import streamlit as st
import yaml
from streamlit_agraph import agraph, Node, Edge, Config
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

    st.divider()
    y["network"] = []
    vessels = st.number_input(
        "number of vessels",
        min_value=1,
        step=1,
    )
    for i in range(vessels):
        v = get_vessel(i)
        if v:
            y["network"].append(v)


with c2:
    st.download_button(
        "Download .yaml",
        data=yaml.dump(
            y,
            sort_keys=False,
            default_flow_style=False,
        ),
        file_name=f"{y['project_name']}.yaml",
    )

    nodes = []
    edges = []
    seen_nodes = []
    for v in y["network"]:
        if v["sn"] not in seen_nodes:
            nodes.append(
                Node(
                    id=v["sn"],
                    label=str(v["sn"]),
                    size=10 if "inlet" in v else 5,
                )
            )
            seen_nodes.append(v["sn"])

        if v["tn"] not in seen_nodes:
            nodes.append(
                Node(
                    id=v["tn"],
                    label=str(v["tn"]),
                    size=1 if "outlet" in v else 5,
                )
            )
            seen_nodes.append(v["tn"])

        edges.append(Edge(source=v["sn"], label=v["label"], target=v["tn"]))

    config = Config(
        width=750,
        height=950,
        directed=True,
        hierarchical=False,
        physics=True,
    )

    agraph(nodes=nodes, edges=edges, config=config)

with c1:
    st.divider()

    if y:
        st.json(y)
