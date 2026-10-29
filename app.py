import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO

# RDKit
from rdkit import Chem
from rdkit.Chem import (
    Draw,
    AllChem,
    Descriptors,
    Crippen,
    Lipinski,
    rdMolDescriptors
)

# 3D Visualization
import py3Dmol
import streamlit.components.v1 as components

# PDF
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.colors import HexColor
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    PageBreak,
    Image,
    Preformatted
)
from reportlab.lib.units import inch
from reportlab.lib.utils import ImageReader
import tempfile
import os

st.set_page_config(
    page_title="CHEM𖦹CHAOS",
    layout="wide"
)

# Initialize theme
if "dark_mode" not in st.session_state:
    st.session_state.dark_mode = False

# Toggle button
if st.sidebar.button(
    "☀️ Switch to Light Mode"
    if st.session_state.dark_mode
    else "🌙 Switch to Dark Mode"
):
    st.session_state.dark_mode = not st.session_state.dark_mode

st.sidebar.markdown("---")
st.sidebar.markdown("### Navigation")
page = st.sidebar.radio(
    "Go to",
    [
        "★ Home",
        "✴︎ Theory",
        "⟡ Molecular Visualization",
        "✽ Molecular Properties Report",
        "༄ Structure–Property Analysis",
        "𓍢ִ໋🀦 Assessment",
        "🗫 Final Project Report"
    ]
)

# Theme colors
if st.session_state.dark_mode:
    background = "#12121F"
    sidebar = "#202035"
    text = "#00FFFF"
    heading = "#39FF14"
    input_background = "#292940"
else:
    background = "#FFFFFF"
    sidebar = "#F0F2F6"
    text = "#222222"
    heading = "#154360"
    input_background = "#FFFFFF"

# Apply theme
st.markdown(
    f"""
    <style>
    .stApp {{
        background-color: {background};
        color: {text};
    }}

    [data-testid="stSidebar"] {{
        background-color: {sidebar};
    }}

    /* Main text */
    .stApp p,
    .stApp label,
    .stApp li,
    .stApp span {{
        color: {text} !important;
    }}

    /* Headings */
    .stApp h1,
    .stApp h2,
    .stApp h3,
    .stApp h4 {{
        color: {heading} !important;
    }}

    /* Sidebar text */
    [data-testid="stSidebar"] * {{
        color: {text} !important;
    }}

    /* Text and number inputs */
    .stTextInput input,
    .stNumberInput input,
    .stTextArea textarea {{
        background-color: {input_background} !important;
        color: {text} !important;
        -webkit-text-fill-color: {text} !important;
    }}

    /* Dropdowns */
    [data-baseweb="select"] {{
        background-color: {input_background} !important;
        color: {text} !important;
    }}

    /* Buttons */
    .stButton button {{
        border: 1px solid {heading};
        border-radius: 8px;
        font-weight: bold;
    }}

    /* Alerts */
    [data-testid="stAlert"] {{
        color: {text} !important;
    }}
    </style>
    """,
    unsafe_allow_html=True
)

st.sidebar.markdown("---")
st.sidebar.markdown(
    "<div style='text-align: center; font-size: 12px; opacity: 0.8;'>"
    "<b>Developed by:</b> CHEM𖦹CHAOS Team"
    "</div>",
    unsafe_allow_html=True
)

# ============================================================
# MAIN PAGE ROUTING
# ============================================================
if page == "★ Home":
    logo_col1, logo_col2, logo_col3 = st.columns([1, 4, 1])
    with logo_col2:
        st.image("1000051916.png", use_container_width=True)
        
    st.title("★ Home")
    st.subheader("Welcome to CHEM𖦹CHAOS - The Virtual Cheminformatics Laboratory")
    st.write(
        "This interactive virtual laboratory helps students understand "
        "important concepts in cheminformatics and computational chemistry."
    )
    st.markdown("### ✴ What can students do? ✴")
    st.markdown(
        "Compare molecular formulas, molecular weights, "
        "hydrogen-bond donors, and acceptors."
    )

    # Sample Data for the Molecular Chart
    data = {
        "Molecule": ["Water", "Ethanol", "Glucose", "Aspirin", "Caffeine"],
        "Formula": ["H2O", "C2H5OH", "C6H12O6", "C9H8O4", "C8H10N4O2"],
        "Molecular Weight (g/mol)": [18.015, 46.07, 180.16, 180.16, 194.19]
    }

    df = pd.DataFrame(data)
    st.dataframe(df, use_container_width=True)
    st.bar_chart(df.set_index("Molecule"), y="Molecular Weight (g/mol)")
    st.caption("Molecular weights are approximate values in g/mol.")

elif page == "✴︎ Theory":
    st.title("✴︎ Theory")
    st.write(
        "Welcome to the theoretical foundation of **CHEM𖦹CHAOS**. "
        "Here we explore how computers represent and analyze chemical structures."
    )

    st.markdown("### 1. SMILES Notation")
    st.write(
        "Simplified Molecular-Input Line-Entry System (SMILES) is a specification in "
        "form of a line notation for describing the structure of chemical species using short ASCII strings."
    )
    st.info("Example: Ethanol is represented as `CCO`, and Benzene as a closed ring `c1ccccc1`.")

    st.markdown("### 2. Key Molecular Descriptors")
    st.markdown(
        "- **LogP (Partition Coefficient):** Measures the lipophilicity of a compound. High LogP values indicate fat-solubility, while negative values indicate water-solubility.\n"
        "- **TPSA (Topological Polar Surface Area):** Measures the surface sum over all polar atoms (usually oxygens and nitrogens) and is a key indicator of membrane permeability.\n"
        "- **Hydrogen Bond Donors/Acceptors:** Important parameters for predicting drug-likeness and intermolecular binding."
    )

    st.markdown("### 3. Lipinski's Rule of Five")
    st.write(
        "A rule of thumb to evaluate if a chemical compound with a certain pharmacological "
        "or biological activity has chemical properties that would make it a likely orally active drug in humans:"
    )
    st.markdown(
        "1. Molecular weight under 500 Daltons\n"
        "2. No more than 5 Hydrogen bond donors\n"
        "3. No more than 10 Hydrogen bond acceptors\n"
        "4. LogP not higher than 5"
    )

    st.markdown("### 4. Fraction $C_{sp^3}$ (Saturation & 3D Character)")
    st.write(
        "Fraction $C_{sp^3}$ ($F_{csp3}$) is the ratio of $sp^3$-hybridized carbons to the total carbon count. "
        "It is a key indicator of a molecule's three-dimensional complexity:"
    )
    st.markdown(
        "- **Low $F_{csp3}$ (near 0):** Indicates a flat, highly aromatic, or unsaturated molecular framework.\n"
        "- **High $F_{csp3}$ (near 1):** Indicates a saturated, more complex 3D shape, which often correlates with improved clinical success and solubility."
    )

elif page == "⟡ Molecular Visualization":
    st.title("⟡ Molecular Visualization")
    st.write("Visualization tools go here.")

elif page == "✽ Molecular Properties Report":
    st.title("✽ Molecular Properties Report")
    st.write("Properties calculation goes here.")

elif page == "༄ Structure–Property Analysis":
    st.title("༄ Structure–Property Analysis")
    st.write("Structure-property analysis goes here.")

elif page == "𓍢ִ໋🀦 Assessment":
    st.title("𓍢ִ໋🀦 Assessment")
    st.write("Test your knowledge on cheminformatics principles and molecular descriptors!")

    st.markdown("#### Question 1: What does LogP primarily measure in cheminformatics?")
    q1 = st.radio(
        "Select your answer for Q1:",
        [
            "A) Molecular weight and mass",
            "B) Lipophilicity (fat vs. water solubility)",
            "C) Number of aromatic rings",
            "D) Exact boiling point"
        ],
        key="q1_ans"
    )

    if st.button("Submit Q1 Answer"):
        if "B) Lipophilicity" in q1:
            st.success("Correct! LogP measures the partition coefficient, indicating how lipophilic a molecule is.")
        else:
            st.error("Incorrect. Remember that LogP relates to partition and lipophilicity.")

elif page == "🗫 Final Project Report":
    st.title("🗫 Final Project Report")
    st.write("PDF report generator goes here.")

# ============================================================
# HELPER FUNCTIONS & STATE
# ============================================================
EXAMPLE_MOLECULES = {
    "Ethanol": "CCO",
    "Benzene": "c1ccccc1",
    "Phenol": "Oc1ccccc1",
    "Acetic Acid": "CC(=O)O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}

if "active_smiles" not in st.session_state:
    st.session_state.active_smiles = "CCO"

if "selected_example" not in st.session_state:
    st.session_state.selected_example = "Ethanol"

def load_shared_example():
    st.session_state.active_smiles = EXAMPLE_MOLECULES[
        st.session_state.selected_example
    ]

def get_molecule(smiles):
    if not smiles or not smiles.strip():
        return None
    try:
        return Chem.MolFromSmiles(smiles.strip())
    except Exception:
        return None

def calculate_bond_information(mol):
    heavy_atom_bonds = mol.GetNumBonds()
    kekule_mol = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule_mol, clearAromaticFlags=True)
    except Exception:
        pass

    pi_bonds = 0
    for bond in kekule_mol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.DOUBLE:
            pi_bonds += 1
        elif bond_type == Chem.BondType.TRIPLE:
            pi_bonds += 2
        elif bond_type == Chem.BondType.QUADRUPLE:
            pi_bonds += 3
        elif bond.GetIsAromatic():
            pi_bonds += 0.5

    mol_with_h = Chem.AddHs(Chem.Mol(mol))
    sigma_bonds = mol_with_h.GetNumBonds()
    total_bonds = sigma_bonds + pi_bonds

    return {
        "Heavy Atom Bonds": heavy_atom_bonds,
        "Sigma Bonds": sigma_bonds,
        "Pi Bonds": pi_bonds,
        "Total Bonds": total_bonds,
        "Total Atoms Including H": mol_with_h.GetNumAtoms(),
    }

def calculate_properties(mol):
    bond_info = calculate_bond_information(mol)
    molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
    properties = {
        "Molecular Formula": molecular_formula,
        "Molecular Weight": round(Descriptors.MolWt(mol), 3),
        "Exact Molecular Weight": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "Hydrogen Bond Donors": Lipinski.NumHDonors(mol),
        "Hydrogen Bond Acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Ring Count": Lipinski.RingCount(mol),
        "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "Aliphatic Rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "Number of Atoms": mol.GetNumAtoms(),
        "Heavy Atoms": mol.GetNumHeavyAtoms(),
        "Heavy Atom Bonds": bond_info["Heavy Atom Bonds"],
        "Sigma (σ) Bonds": bond_info["Sigma Bonds"],
        "Pi (π) Bonds": bond_info["Pi Bonds"],
        "Total Bonds (σ + π)": bond_info["Total Bonds"],
        "Fraction Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "Molar Refractivity": round(Crippen.MolMR(mol), 3)
    }
    return properties

def generate_2d_image(mol):
    return Draw.MolToImage(mol, size=(700, 500))

def generate_3d_molecule(mol):
    try:
        mol3d = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(mol3d, params)
        if status != 0:
            return None
        try:
            AllChem.MMFFOptimizeMolecule(mol3d)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol3d)
            except Exception:
                pass
        return mol3d
    except Exception:
        return None

def display_3d_structure(mol3d, style, show_surface=False):
    mol_block = Chem.MolToMolBlock(mol3d)
    view = py3Dmol.view(width=900, height=600)
    view.addModel(mol_block, "mol")

    if style == "Ball and Stick":
        view.setStyle({"stick": {"radius": 0.15}, "sphere": {"scale": 0.30}})
    elif style == "Stick":
        view.setStyle({"stick": {"radius": 0.20}})
    elif style == "Space Filling":
        view.setStyle({"sphere": {"scale": 1.0}})
    elif style == "Wireframe":
        view.setStyle({"line": {"linewidth": 2}})
    elif style == "Cartoon":
        view.setStyle({"stick": {}, "sphere": {"scale": 0.20}})

    if show_surface:
        view.addSurface(py3Dmol.VDW, {"opacity": 0.7})

    view.setBackgroundColor("white")
    view.zoomTo()
    html = view._make_html()
    components.html(html, height=620)

def generate_molecule_pdf(input_smiles, mol, properties):
    buffer = BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        leftMargin=45,
        rightMargin=45,
        topMargin=45,
        bottomMargin=45
    )
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "TitleStyle",
        parent=styles["Title"],
        fontSize=20,
        alignment=TA_CENTER,
        textColor=HexColor("#17365D")
    )
    heading_style = ParagraphStyle(
        "HeadingStyle",
        parent=styles["Heading1"],
        fontSize=14,
        textColor=HexColor("#17365D")
    )
    body_style = ParagraphStyle(
        "BodyStyle",
        parent=styles["BodyText"],
        fontSize=10,
        leading=14,
        alignment=TA_JUSTIFY
    )

    story = []
    story.append(Paragraph("MOLECULAR PROPERTIES REPORT", title_style))
    story.append(Spacer(1, 0.3 * inch))

    story.append(Paragraph("1. Molecular Information", heading_style))
    canonical_smiles = Chem.MolToSmiles(mol)
    basic_data = [
        ["Property", "Value"],
        ["Input SMILES", input_smiles],
        ["Canonical SMILES", canonical_smiles],
        ["Molecular Formula", properties["Molecular Formula"]]
    ]
    basic_table = Table(basic_data, colWidths=[2 * inch, 4 * inch])
    basic_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#17365D")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("PADDING", (0, 0), (-1, -1), 7)
        ])
    )
    story.append(basic_table)
    story.append(Spacer(1, 0.3 * inch))

    story.append(Paragraph("2. Molecular Properties", heading_style))
    property_data = [["Property", "Value"]]
    for key, value in properties.items():
        property_data.append([key, str(value)])

    property_table = Table(property_data, colWidths=[3 * inch, 3 * inch], repeatRows=1)
    property_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#2F5597")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("PADDING", (0, 0), (-1, -1), 6)
        ])
    )
    story.append(property_table)
    story.append(Spacer(1, 0.3 * inch))

    story.append(Paragraph("3. Interpretation", heading_style))
    interpretation = "This report was generated automatically using CHEM𖦹CHAOS."
    story.append(Paragraph(interpretation, body_style))

    doc.build(story)
    pdf_data = buffer.getvalue()
    buffer.close()
    return pdf_data

def get_property_interpretation(property_name, value, properties):
    if property_name == "Molecular Formula":
        return "Shows the elemental composition of the molecule."
    elif property_name == "Molecular Weight":
        return "Represents the average molecular mass."
    elif property_name == "LogP":
        return "Indicates lipophilicity."
    else:
        return "Property calculated successfully."

