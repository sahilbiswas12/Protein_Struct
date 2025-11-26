# =====================================================
# 🧬 Protein Structure Prediction and Visualization App
# Streamlit Version (Stable, Tested on Python 3.12)
# =====================================================

import streamlit as st
import requests, gzip
from io import BytesIO
from collections import Counter
from matplotlib import pyplot as plt
import py3Dmol

# ----------------- Streamlit Page Setup ----------------- #
st.set_page_config(
    page_title="Protein Structure Prediction and Visualization",
    layout="wide",
    page_icon="🧬"
)

# Custom Styling
st.markdown("""
    <style>
    [data-testid="stAppViewContainer"] {
        background: linear-gradient(to bottom right, #c8e9f8, #eaf3fc);
        background-attachment: fixed;
    }
    h1, h2, h3 {
        color: #004d80;
    }
    .stTabs [data-baseweb="tab-list"] {
        justify-content: center;
    }
    </style>
""", unsafe_allow_html=True)

# ----------------- Helper Functions ----------------- #
AMINO_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

@st.cache_data(show_spinner=False)
def get_proteome_data(species, max_seq=200):
    """Fetch Swiss-Prot reviewed sequences from UniProt REST API."""
    proteome_ids = {
        "Human": "UP000005640",
        "Mouse": "UP000000589",
        "Fruit Fly": "UP000000803",
        "E. coli": "UP000000625",
        "Yeast": "UP000002311",
    }
    base_url = "https://rest.uniprot.org"
    pid = proteome_ids[species]
    url = f"{base_url}/uniprotkb/stream?format=fasta&query=(proteome:{pid})+AND+(reviewed:true)&compressed=true"

    try:
        r = requests.get(url, stream=True, timeout=120)
        r.raise_for_status()
        data = b"".join(r.iter_content(8192))
        with gzip.GzipFile(fileobj=BytesIO(data)) as gz:
            fasta = gz.read().decode("utf-8")
    except Exception as e:
        st.error(f"Error fetching data: {e}")
        return []

    seqs, header, sequence = [], None, []
    for line in fasta.splitlines():
        if line.startswith(">"):
            if header:
                seqs.append({"header": header, "sequence": "".join(sequence)})
            header, sequence = line[1:], []
        else:
            sequence.append(line.strip())
    if header:
        seqs.append({"header": header, "sequence": "".join(sequence)})

    data = []
    for s in seqs[:max_seq]:
        header = s["header"]
        seq = s["sequence"]
        if len(seq) < 20:
            continue
        p0 = header.split("|")
        info = {
            "uniprot_id": p0[1] if len(p0) >= 3 else "Unknown",
            "protein_name": header.split(" OS=")[0],
            "organism": (header.split("OS=")[1].split("OX=")[0].strip() if "OS=" in header else "Unknown"),
            "gene_name": (header.split("GN=")[1].split()[0] if "GN=" in header else "Unknown"),
            "sequence": seq,
            "length": len(seq)
        }
        data.append(info)
    return data


def aa_composition(seq):
    """Compute amino acid composition percentages."""
    c = Counter(seq)
    return [100 * c.get(a, 0) / len(seq) for a in AMINO_ORDER]


def show_3d_structure(uniprot_id):
    """Visualize protein 3D model from SWISS-MODEL."""
    url = f"https://swissmodel.expasy.org/repository/uniprot/{uniprot_id}.pdb"
    r = requests.get(url, timeout=60)
    if r.status_code != 200 or "ATOM" not in r.text:
        st.warning("No 3D structure available for this protein.")
        return
    view = py3Dmol.view(width=800, height=500)
    view.addModel(r.text, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=520, scrolling=False)


# ----------------- Main Layout ----------------- #
st.title("🧬 Protein Structure Prediction and Visualization")

tabs = st.tabs(["Overview", "Sequences", "3D Viewer", "About"])

# ----------------- Tab 1: Overview ----------------- #
with tabs[0]:
    st.header("Species Selection and Data Fetching")

    col1, col2, col3 = st.columns([1, 1, 2])
    species = col1.selectbox("Select Species", ["Human", "Mouse", "Fruit Fly", "E. coli", "Yeast"])
    max_seq = col2.number_input("Max Sequences", 10, 500, 100)
    fetch_btn = col3.button("Fetch Data", use_container_width=True)

    if fetch_btn:
        with st.spinner("Fetching data..."):
            data = get_proteome_data(species, max_seq)
        if data:
            st.session_state["proteins"] = data
            st.success(f"Loaded {len(data)} proteins for {species}.")

            avg_len = round(sum(p["length"] for p in data) / len(data), 1)
            max_len = max(p["length"] for p in data)
            st.metric("Total Proteins", len(data))
            st.metric("Average Length", f"{avg_len} aa")
            st.metric("Longest Protein", f"{max_len} aa")
        else:
            st.error("No data returned. Check internet connection or UniProt availability.")


# ----------------- Tab 2: Sequence Browser ----------------- #
with tabs[1]:
    st.header("Protein Sequence Prediction and Analysis")

    if "proteins" not in st.session_state:
        st.warning("Please fetch data first from the Overview tab.")
    else:
        proteins = st.session_state["proteins"]
        names = [f"{p['uniprot_id']} | {p['protein_name']} | {p['length']} aa" for p in proteins]
        selected = st.selectbox("Select a Protein", names)

        if selected:
            p = proteins[names.index(selected)]
            st.subheader(f"{p['protein_name']} ({p['uniprot_id']})")
            st.markdown(f"""
            **Organism:** {p['organism']}  
            **Gene:** {p['gene_name']}  
            **Length:** {p['length']} amino acids
            """)

            # Amino acid composition plot
            comp = aa_composition(p["sequence"])
            fig, ax = plt.subplots(figsize=(8, 3))
            ax.bar(AMINO_ORDER, comp, color="#007acc")
            ax.set_ylabel("%")
            ax.set_title(f"Amino Acid Composition — {p['uniprot_id']}")
            st.pyplot(fig)

            # Sequence text box
            seq_display = "\n".join(p["sequence"][i:i+80] for i in range(0, len(p["sequence"]), 80))
            st.text_area("Protein Sequence", seq_display, height=250)


# ----------------- Tab 3: 3D Viewer ----------------- #
with tabs[2]:
    st.header("3D Structure Viewer")

    if "proteins" not in st.session_state:
        st.warning("Fetch Swiss-Prot data first.")
    else:
        proteins = st.session_state["proteins"]
        ids = [p["uniprot_id"] for p in proteins]
        selected_id = st.selectbox("Select UniProt ID", ids)
        if st.button("Load 3D Structure"):
            with st.spinner("Loading 3D model..."):
                show_3d_structure(selected_id)


# ----------------- Tab 4: About ----------------- #
with tabs[3]:
    st.header("ℹ️ About the Application")
    st.markdown("""
    ### Protein Structure Prediction and Visualization (2025)
    **An interactive educational app** for exploring reviewed Swiss-Prot proteins from UniProt.

    #### 🔍 Features
    - Fetch predicted protein sequence data for Human, Mouse, Yeast, etc.
    - Analyze amino acid composition
    - Visualize 3D protein models 
    - Built for classroom learning and bioinformatics projects

    #### 🧪 Technologies Used
    - Python • Streamlit • Matplotlib • Py3Dmol • Custom ProteinGPT LLM

    #### 👨‍💻 Developers
    Sahil Biswas • Ayaan Bhoje • Arfat Fakih • Raghav Deshpande (2025)
    """)