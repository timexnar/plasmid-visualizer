"""
Streamlit web app tying together plasmid parsing, restriction site
detection, insert matching, and both map styles (static and interactive)
into one browser-based UI.

Run with: .venv/bin/streamlit run app.py
"""

import tempfile
from pathlib import Path

import streamlit as st

from insert_matching import build_insert_features, find_insert_matches
from plasmid_io import load_record
from plot_plasmid import build_static_figure
from plot_plasmid_interactive import build_interactive_figure
from restriction_sites import find_unique_cutters

DEFAULT_FILE = "addgene-plasmid-50005-sequence-513490.gbk"
ALLOWED_UPLOAD_TYPES = ["gb", "gbk", "genbank", "fasta", "fa", "fna", "txt"]

st.set_page_config(page_title="Plasmid Visualizer", layout="wide")
st.title("Plasmid Visualizer")

uploaded_file = st.sidebar.file_uploader(
    "Plasmid file (GenBank, FASTA, or plain-text sequence)",
    type=ALLOWED_UPLOAD_TYPES,
)
insert_seq = st.sidebar.text_area("Insert sequence to highlight (optional)", height=100)
map_style = st.sidebar.radio("Map style", ["Static", "Interactive"])

if uploaded_file is not None:
    suffix = Path(uploaded_file.name).suffix
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp_file:
        tmp_file.write(uploaded_file.getvalue())
        input_path = tmp_file.name
else:
    input_path = DEFAULT_FILE
    st.sidebar.caption(f"No file uploaded - using the bundled sample plasmid ({DEFAULT_FILE}).")

record = load_record(input_path)

insert_seq = insert_seq.strip()
if insert_seq:
    matches = find_insert_matches(record.seq, insert_seq)
    record.features.extend(build_insert_features(matches, len(record.seq)))
    if matches:
        st.sidebar.success(f"Insert found: {len(matches)} match(es)")
    else:
        st.sidebar.warning("No match found for the given insert sequence.")

map_col, info_col = st.columns([2, 1])

with map_col:
    st.subheader(record.name)
    st.caption(f"{len(record.seq)} bp")
    if map_style == "Static":
        st.pyplot(build_static_figure(record))
    else:
        st.plotly_chart(build_interactive_figure(record), use_container_width=True)

with info_col:
    st.subheader("Features")
    if record.features:
        st.dataframe(
            [
                {
                    "type": feature.type,
                    "label": feature.qualifiers.get("label", ["(unnamed)"])[0],
                    "start": int(feature.location.start) + 1,
                    "end": int(feature.location.end),
                    "strand": {1: "+", -1: "-"}.get(feature.location.strand, "?"),
                }
                for feature in record.features
            ],
            hide_index=True,
        )
    else:
        st.write("No annotated features found in this file.")

    st.subheader("Unique restriction sites")
    cut_sites = find_unique_cutters(record.seq)
    if cut_sites:
        st.dataframe(
            [{"enzyme": name, "position": position} for name, position in cut_sites],
            hide_index=True,
        )
    else:
        st.write("None found among the common cloning enzymes checked.")
