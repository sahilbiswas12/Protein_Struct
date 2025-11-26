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
import base64
import os

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
        background: linear-gradient(to bottom right, #0f172a, #1e3a8a, #0f172a);
        background-attachment: fixed;
    }
    h1, h2, h3 {
        color: #ffffff;
    }
    p, div, span, label {
        color: #e0e0e0;
    }
    .stTabs [data-baseweb="tab-list"] {
        justify-content: center;
    }
    
    /* Navbar Styles */
    html, body {
        margin: 0 !important;
        padding: 0 !important;
        overflow-x: hidden;
    }
    [data-testid="stHeader"] {
        display: none !important;
        height: 0 !important;
        visibility: hidden !important;
    }
    header[data-testid="stHeader"] {
        display: none !important;
    }
    #MainMenu {
        visibility: hidden;
    }
    footer {
        visibility: hidden;
    }
    .stDeployButton {
        display: none;
    }
    [data-testid="stAppViewContainer"] {
        padding-top: 0 !important;
        margin-top: 0 !important;
    }
    .main .block-container {
        padding-top: 5.5rem !important;
        margin-top: 0 !important;
    }
    section[data-testid="stSidebar"] {
        top: 4.5rem !important;
    }
    .navbar {
        display: flex;
        justify-content: space-between;
        align-items: center;
        padding: 1rem 2rem;
        background: #002366;
        box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        width: 100%;
        z-index: 9999;
        margin: 0 !important;
        padding-top: 1rem !important;
        padding-bottom: 1rem !important;
        border-radius: 0;
        box-sizing: border-box;
    }
    .navbar-logo {
        display: flex;
        align-items: center;
        font-size: 1.5rem;
        font-weight: bold;
        color: white;
        text-decoration: none;
    }
    .navbar-logo-icon {
        font-size: 2rem;
        margin-right: 0.5rem;
    }
    .navbar-menu {
        display: flex;
        gap: 1.5rem;
        align-items: center;
    }
    .navbar-link {
        color: white !important;
        text-decoration: none;
        font-weight: 500;
        padding: 0.5rem 1rem;
        border-radius: 5px;
        transition: background-color 0.3s ease;
    }
    .navbar-link:hover {
        background-color: rgba(255, 255, 255, 0.2);
        color: white !important;
    }
    .hero-section {
        background: linear-gradient(to bottom right, #0f172a, #1e3a8a, #0f172a);
        padding: 4rem 2rem;
        margin: 0 !important;
        margin-top: -5.5rem !important;
        position: relative;
        left: 50%;
        right: 50%;
        margin-left: -50vw !important;
        margin-right: -50vw !important;
        width: 100vw;
        display: flex;
        align-items: center;
        min-height: 100vh;
        box-sizing: border-box;
        padding-top: calc(4.5rem + 4rem);
    }
    .hero-content {
        max-width: 1200px;
        margin: 0 auto;
        width: 100%;
        display: flex;
        gap: 3rem;
        align-items: center;
    }
    .hero-left {
        flex: 1;
    }
    .hero-right {
        flex: 1;
        display: flex;
        justify-content: center;
        align-items: center;
    }
    .hero-tag {
        display: inline-flex;
        align-items: center;
        gap: 0.5rem;
        padding: 0.5rem 1rem;
        background-color: rgba(15, 23, 42, 0.6);
        border: 1px solid #3b82f6;
        border-radius: 20px;
        color: white;
        font-size: 0.875rem;
        font-weight: 500;
        margin-bottom: 1.5rem;
    }
    .hero-tag-icon {
        font-size: 1rem;
        color: #60a5fa;
    }
    .hero-header {
        margin: 0;
        text-align: left;
        line-height: 1.2;
        margin-bottom: 1.5rem;
    }
    .hero-header-part1 {
        font-size: 3.5rem;
        font-weight: bold;
        color: white;
        display: block;
    }
    .hero-header-part2 {
        font-size: 4rem;
        font-weight: bold;
        color: #60a5fa;
        display: block;
    }
    .hero-description {
        color: white;
        font-size: 1.125rem;
        line-height: 1.6;
        margin-bottom: 2rem;
        max-width: 600px;
    }
    .hero-buttons {
        display: flex;
        gap: 1rem;
        margin-bottom: 3rem;
    }
    .hero-btn-primary {
        padding: 0.875rem 2rem;
        background: linear-gradient(to right, #3b82f6, #06b6d4);
        color: white;
        border: none;
        border-radius: 8px;
        font-size: 1rem;
        font-weight: 600;
        cursor: pointer;
        transition: transform 0.2s, box-shadow 0.2s;
    }
    .hero-btn-primary:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(59, 130, 246, 0.4);
    }
    .hero-btn-secondary {
        padding: 0.875rem 2rem;
        background-color: rgba(15, 23, 42, 0.8);
        color: white;
        border: 1px solid #64748b;
        border-radius: 8px;
        font-size: 1rem;
        font-weight: 600;
        cursor: pointer;
        transition: background-color 0.2s;
    }
    .hero-btn-secondary:hover {
        background-color: rgba(15, 23, 42, 1);
    }
    .hero-stats {
        display: flex;
        gap: 3rem;
        margin-top: 2rem;
    }
    .hero-stat {
        display: flex;
        flex-direction: column;
    }
    .hero-stat-number {
        font-size: 2.5rem;
        font-weight: bold;
        color: #60a5fa;
        line-height: 1;
        margin-bottom: 0.5rem;
    }
    .hero-stat-label {
        font-size: 1rem;
        color: white;
        font-weight: 500;
    }
    @media (max-width: 968px) {
        .hero-content {
            flex-direction: column;
            gap: 2rem;
        }
        .hero-left, .hero-right {
            flex: 1;
            width: 100%;
        }
        .hero-header-part1 {
            font-size: 2.5rem;
        }
        .hero-header-part2 {
            font-size: 3rem;
        }
    }
    .section-container {
        margin: 3rem 0;
        padding: 2rem 0;
        border-bottom: 2px solid rgba(96, 165, 250, 0.3);
    }
    .section-container:last-child {
        border-bottom: none;
    }
    .about-header-center {
        text-align: center;
    }
    .about-cards-container {
        display: flex;
        flex-direction: column;
        gap: 1.5rem;
        margin-top: 2rem;
        max-width: 1200px;
        margin-left: auto;
        margin-right: 0;
        padding-left: 1rem;
    }
    .about-card {
        background: rgba(15, 23, 42, 0.8);
        border-radius: 12px;
        padding: 2rem;
        display: flex;
        align-items: flex-start;
        gap: 1.5rem;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
        transition: transform 0.2s, box-shadow 0.2s;
        max-width: 600px;
    }
    .about-card-left {
        align-self: flex-start;
        margin-right: auto;
        margin-left: 1.5rem;
    }
    .about-card-right {
        align-self: flex-end;
        margin-left: auto;
        margin-right: 0;
    }
    .about-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 25px rgba(0, 0, 0, 0.4);
    }
    .about-card-icon {
        font-size: 2.5rem;
        color: #60a5fa;
        flex-shrink: 0;
        width: 64px;
        height: 64px;
        display: flex;
        align-items: center;
        justify-content: center;
        background: rgba(96, 165, 250, 0.15);
        border-radius: 16px;
    }
    .about-card-content {
        flex: 1;
    }
    .about-card-title {
        font-size: 1.5rem;
        font-weight: bold;
        color: #ffffff;
        margin: 0 0 0.75rem 0;
    }
    .about-card-description {
        color: #d1d5db;
        font-size: 1rem;
        line-height: 1.6;
        margin: 0;
    }
    .hero-carousel {
        position: relative;
        width: 100%;
        max-width: 500px;
        display: flex;
        flex-direction: column;
        align-items: center;
        gap: 1.5rem;
    }
    .carousel-container {
        position: relative;
        width: 100%;
        height: 400px;
        overflow: hidden;
        border-radius: 12px;
        box-shadow: 0 8px 32px rgba(0, 0, 0, 0.3);
    }
    .carousel-slide {
        position: absolute;
        width: 100%;
        height: 100%;
        opacity: 0;
        transition: opacity 0.5s ease-in-out;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    .carousel-slide.active {
        opacity: 1;
        z-index: 1;
    }
    .carousel-slide img {
        width: 100%;
        height: 100%;
        object-fit: cover;
        border-radius: 12px;
    }
    .carousel-indicators {
        display: flex;
        gap: 0.75rem;
        justify-content: center;
        align-items: center;
    }
    .carousel-indicator {
        width: 12px;
        height: 12px;
        border-radius: 50%;
        background-color: rgba(255, 255, 255, 0.4);
        border: 2px solid rgba(255, 255, 255, 0.6);
        cursor: pointer;
        transition: all 0.3s ease;
    }
    .carousel-indicator:hover {
        background-color: rgba(255, 255, 255, 0.6);
        transform: scale(1.2);
    }
    .carousel-indicator.active {
        background-color: #60a5fa;
        border-color: #60a5fa;
        width: 14px;
        height: 14px;
    }
    </style>
""", unsafe_allow_html=True)

# Helper function to load images as base64
def get_image_base64(image_filename):
    """Load image file and return as base64 string."""
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    image_path = os.path.join(script_dir, image_filename)
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode()
    except FileNotFoundError:
        # Try relative path from current working directory
        try:
            with open(image_filename, "rb") as img_file:
                return base64.b64encode(img_file.read()).decode()
        except FileNotFoundError:
            st.error(f"Image not found: {image_filename}")
            return ""

# Navbar HTML
st.markdown("""
    <nav class="navbar">
        <div class="navbar-logo">
            <span class="navbar-logo-icon">🧬</span>
            <span>ProteinStruct</span>
        </div>
        <div class="navbar-menu">
            <a href="#overview" class="navbar-link">Overview</a>
            <a href="#sequences" class="navbar-link">Sequences</a>
            <a href="#3d-viewer" class="navbar-link">3D Viewer</a>
            <a href="#about" class="navbar-link">About</a>
        </div>
    </nav>
""", unsafe_allow_html=True)

# Hero Section - Load images first
img1_base64 = get_image_base64('img_1.png')
img2_base64 = get_image_base64('img_2.png')

hero_html = f"""
    <div class="hero-section">
        <div class="hero-content">
            <div class="hero-left">
                <div class="hero-tag">
                    <span class="hero-tag-icon">🧬</span>
                    <span>Bioinformatics Platform</span>
                </div>
                <h1 class="hero-header">
                    <span class="hero-header-part1">Explore Protein</span>
                    <span class="hero-header-part2">Structures</span>
                </h1>
                <p class="hero-description">
                    Visualize and analyze curated protein sequences from Swiss-Prot. Discover amino acid compositions and explore 3D molecular structures with our interactive platform.
                </p>
                <div class="hero-buttons">
                    <button class="hero-btn-primary">Get Started</button>
                    <button class="hero-btn-secondary">Learn More</button>
                </div>
                <div class="hero-stats">
                    <div class="hero-stat">
                        <div class="hero-stat-number">5+</div>
                        <div class="hero-stat-label">Species</div>
                    </div>
                    <div class="hero-stat">
                        <div class="hero-stat-number">200+</div>
                        <div class="hero-stat-label">Proteins</div>
                    </div>
                    <div class="hero-stat">
                        <div class="hero-stat-number">3D</div>
                        <div class="hero-stat-label">Visualization</div>
                    </div>
                </div>
            </div>
            <div class="hero-right">
                <div class="hero-carousel">
                    <div class="carousel-container">
                        <div class="carousel-slide active" id="slide-0">
                            <img src="data:image/png;base64,{img1_base64}" alt="Protein Structure 1">
                        </div>
                        <div class="carousel-slide" id="slide-1">
                            <img src="data:image/png;base64,{img2_base64}" alt="Protein Structure 2">
                        </div>
                    </div>
                    <div class="carousel-indicators">
                        <div class="carousel-indicator active" onclick="changeSlide(0)"></div>
                        <div class="carousel-indicator" onclick="changeSlide(1)"></div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    <script>
        function changeSlide(index) {{
            const slides = document.querySelectorAll('.carousel-slide');
            const indicators = document.querySelectorAll('.carousel-indicator');
            
            slides.forEach((slide, i) => {{
                if (i === index) {{
                    slide.classList.add('active');
                }} else {{
                    slide.classList.remove('active');
                }}
            }});
            
            indicators.forEach((indicator, i) => {{
                if (i === index) {{
                    indicator.classList.add('active');
                }} else {{
                    indicator.classList.remove('active');
                }}
            }});
        }}
        
        // Auto-rotate carousel every 5 seconds
        let currentSlide = 0;
        setInterval(() => {{
            currentSlide = (currentSlide + 1) % 2;
            changeSlide(currentSlide);
        }}, 5000);
    </script>
"""
st.markdown(hero_html, unsafe_allow_html=True)

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
# ----------------- Section 1: About ----------------- #
st.markdown('<div id="about"></div>', unsafe_allow_html=True)
st.markdown("""
    <div class="section-container">
""", unsafe_allow_html=True)

st.markdown("""
    <div class="about-header-center">
        <h1>ℹ️ About the Application</h1>
    </div>
""", unsafe_allow_html=True)
st.markdown("""
    <div class="about-cards-container">
        <div class="about-card about-card-left">
            <div class="about-card-icon">🔍</div>
            <div class="about-card-content">
                <h3 class="about-card-title">Features</h3>
                <p class="about-card-description">Fetch predicted protein sequence data for Human, Mouse, Yeast, etc. Analyze amino acid composition. Visualize 3D protein models. Built for classroom learning and bioinformatics projects.</p>
            </div>
        </div>
        <div class="about-card about-card-right">
            <div class="about-card-icon">🧪</div>
            <div class="about-card-content">
                <h3 class="about-card-title">Technologies Used</h3>
                <p class="about-card-description">Python • Streamlit • Matplotlib • Py3Dmol • Custom ProteinGPT LLM</p>
            </div>
        </div>
        <div class="about-card about-card-left">
            <div class="about-card-icon">👨‍💻</div>
            <div class="about-card-content">
                <h3 class="about-card-title">Developers</h3>
                <p class="about-card-description">Sahil Biswas • Ayaan Bhoje • Arfat Fakih • Raghav Deshpande (2025)</p>
            </div>
        </div>
    </div>
""", unsafe_allow_html=True)

st.markdown("</div>", unsafe_allow_html=True)

# ----------------- Section 2: Overview ----------------- #
st.markdown('<div id="overview"></div>', unsafe_allow_html=True)
st.markdown("""
    <div class="section-container">
""", unsafe_allow_html=True)

st.header("📊 Overview - Species Selection and Data Fetching")

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

st.markdown("</div>", unsafe_allow_html=True)

# ----------------- Section 3: Sequence Browser ----------------- #
st.markdown('<div id="sequences"></div>', unsafe_allow_html=True)
st.markdown("""
    <div class="section-container">
""", unsafe_allow_html=True)

st.header("🧬 Sequences - Protein Sequence Prediction and Analysis")

if "proteins" not in st.session_state:
    st.warning("Please fetch data first from the Overview section above.")
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

st.markdown("</div>", unsafe_allow_html=True)

# ----------------- Section 4: 3D Viewer ----------------- #
st.markdown('<div id="3d-viewer"></div>', unsafe_allow_html=True)
st.markdown("""
    <div class="section-container">
""", unsafe_allow_html=True)

st.header("🎯 3D Viewer - 3D Structure Viewer")

if "proteins" not in st.session_state:
    st.warning("Fetch Swiss-Prot data first from the Overview section above.")
else:
    proteins = st.session_state["proteins"]
    ids = [p["uniprot_id"] for p in proteins]
    selected_id = st.selectbox("Select UniProt ID", ids)
    if st.button("Load 3D Structure"):
        with st.spinner("Loading 3D model..."):
            show_3d_structure(selected_id)

st.markdown("</div>", unsafe_allow_html=True)