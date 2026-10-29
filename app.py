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
    .stApp p, .stApp label, .stApp li, .stApp span {{
        color: {text} !important;
    }}
    .stApp h1, .stApp h2, .stApp h3, .stApp h4 {{
        color: {heading} !important;
    }}
    [data-testid="stSidebar"] * {{
        color: {text} !important;
    }}
    .stTextInput input, .stNumberInput input, .stTextArea textarea {{
        background-color: {input_background} !important;
        color: {text} !important;
        -webkit-text-fill-color: {text} !important;
    }}
    [data-baseweb="select"] {{
        background-color: {input_background} !important;
        color: {text} !important;
    }}
    .stButton button {{
        border: 1px solid {heading};
        border-radius: 8px;
        font-weight: bold;
    }}
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
    interpretation = "This report was generated automatically using the CHEM𖦹CHAOS application via RDKit toolkit."
    story.append(Paragraph(interpretation, body_style))

    doc.build(story)
    pdf_data = buffer.getvalue()
    buffer.close()
    return pdf_data

def get_property_interpretation(property_name, value, properties):
    if property_name == "Molecular Formula":
        return "Shows the elemental composition of the molecule and is useful for understanding its basic chemical composition."
    elif property_name == "Molecular Weight":
        return "Represents the average molecular mass. Molecular size and mass can influence diffusion, formulation and many physicochemical behaviours."
    elif property_name == "LogP":
        if value < 0:
            return "A negative value indicates a preference toward the aqueous phase and relatively high hydrophilicity."
        elif value < 3:
            return "Indicates moderate lipophilicity, reflecting a balance between hydrophilic and hydrophobic character."
        else:
            return "Indicates relatively high lipophilicity and a stronger preference for non-polar environments."
    elif property_name == "TPSA":
        if value < 60:
            return "A relatively low polar surface area, indicating limited overall molecular polarity."
        else:
            return "A high polar surface area, indicating substantial polarity and strong hydrogen-bonding potential."
    elif property_name == "Hydrogen Bond Donors":
        return f"The molecule has {value} donor site(s) capable of donating hydrogen in hydrogen-bond interactions."
    elif property_name == "Hydrogen Bond Acceptors":
        return f"The molecule has {value} acceptor site(s) capable of accepting hydrogen bonds."
    elif property_name == "Fraction Csp3":
        return "Ratio of sp3-hybridized carbons indicating 3D structural complexity."
    else:
        return "Standard physicochemical descriptor evaluated via RDKit."

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
        "Enter any SMILES notation\n\n"
        "☆ Validate molecular structures\n\n"
        "☆ Generate 2D structures\n\n"
        "☆ Generate interactive 3D structures\n\n"
        "☆ Select different molecular visualization styles\n\n"
        "☆ Calculate complete molecular properties\n\n"
        "☆ Download molecular reports\n\n"
        "☆ View molecular charts\n\n"
        "☆ Compare multiple molecules\n\n"
        "☆ Complete an assessment"
    )

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
        "- **TPSA (Topological Polar Surface Area):** Measures the surface sum over all polar atoms and is a key indicator of membrane permeability.\n"
        "- **Hydrogen Bond Donors/Acceptors:** Important parameters for predicting drug-likeness."
    )

elif page == "⟡ Molecular Visualization":
    st.title("⟡ Molecular Visualization Studio")
    st.write("Enter any SMILES notation to inspect its 2D topology and explore its 3D spatial conformation.")
    
    col1, col2 = st.columns([2, 1])
    with col1:
        smiles_input = st.text_input("Enter SMILES String", value=st.session_state.active_smiles)
    with col2:
        st.selectbox("Or select an example", list(EXAMPLE_MOLECULES.keys()), key="selected_example")
        if st.button("Load Example Molecule"):
            load_shared_example()
            st.rerun()

    mol = get_molecule(smiles_input)
    if mol:
        st.session_state.active_smiles = smiles_input
        st.success("Molecule validated successfully via RDKit!")
        
        tab1, tab2 = st.tabs(["2D Structure View", "Interactive 3D Visualizer"])
        
        with tab1:
            st.subheader("2D Graph Representation")
            img = generate_2d_image(mol)
            st.image(img, caption="2D Chemical Structure Render", use_container_width=True)
            
        with tab2:
            st.subheader("3D Conformation Visualizer")
            style_col, surf_col = st.columns(2)
            with style_col:
                style = st.selectbox("Render Style", ["Ball and Stick", "Stick", "Space Filling", "Wireframe", "Cartoon"])
            with surf_col:
                show_surf = st.checkbox("Enable Van der Waals Surface", value=False)
                
            mol3d = generate_3d_molecule(mol)
            if mol3d:
                display_3d_structure(mol3d, style, show_surf)
            else:
                st.error("Failed to generate 3D coordinates for this structure.")
    else:
        st.error("Invalid SMILES notation entered. Please check your syntax.")

elif page == "✽ Molecular Properties Report":
    st.title("✽ Molecular Properties & Report Generator")
    st.write("Calculate comprehensive cheminformatics properties and download formatted PDF laboratory reports.")
    
    smiles_input = st.text_input("Enter SMILES for Property Analysis", value=st.session_state.active_smiles, key="prop_input")
    mol = get_molecule(smiles_input)
    
    if mol:
        st.session_state.active_smiles = smiles_input
        properties = calculate_properties(mol)
        
        st.subheader("Calculated Physicochemical Properties")
        prop_df = pd.DataFrame(list(properties.items()), columns=["Property Name", "Value"])
        st.dataframe(prop_df, use_container_width=True)
        
        st.subheader("Detailed Property Interpretations")
        for k, v in properties.items():
            interpretation = get_property_interpretation(k, v, properties)
            st.markdown(f"- **{k} ({v}):** {interpretation}")
            
        st.markdown("---")
        st.subheader("Export Laboratory Report")
        pdf_data = generate_molecule_pdf(smiles_input, mol, properties)
        st.download_button(
            label="📄 Download PDF Molecular Report",
            data=pdf_data,
            file_name="CHEMOCHAOS_Molecular_Report.pdf",
            mime="application/pdf"
        )
    else:
        st.warning("Please provide a valid SMILES string above.")

elif page == "༄ Structure–Property Analysis":
    st.title("༄ Structure–Property Comparative Analysis")
    st.write("Compare multiple chemical structures simultaneously to evaluate structure-activity relationships.")
    
    default_batch = "CCO\nCC(=O)Oc1ccccc1C(=O)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C\nc1ccccc1"
    batch_input = st.text_area("Enter SMILES strings (one per line)", value=default_batch, height=120)
    
    lines = [line.strip() for line in batch_input.split("\n") if line.strip()]
    if lines:
        rows = []
        for line in lines:
            m = get_molecule(line)
            if m:
                p = calculate_properties(m)
                p["SMILES"] = line
                rows.append(p)
                
        if rows:
            comparison_df = pd.DataFrame(rows)
            cols = ["SMILES"] + [c for c in comparison_df.columns if c != "SMILES"]
            comparison_df = comparison_df[cols]
            st.dataframe(comparison_df, use_container_width=True)
            st.bar_chart(comparison_df.set_index("SMILES"), y="Molecular Weight")
        else:
            st.warning("No valid molecules could be parsed from the input list.")

elif page == "𓍢ִ໋🀦 Assessment":
    st.title("🗫 Cheminformatics Knowledge Assessment")
    st.write("Test your understanding of molecular descriptors, cheminformatics rules, and bond properties through this 5-question evaluation.")

    # Question 1
    st.markdown("#### Question 1: What does LogP primarily measure in cheminformatics?")
    q1 = st.radio(
        "Select your answer for Q1:",
        [
            "A) Molecular weight and total atomic mass",
            "B) Lipophilicity (fat vs. water solubility partition coefficient)",
            "C) Number of conjugated aromatic rings",
            "D) Exact boiling point temperature"
        ],
        key="q1_ans"
    )
    if st.button("Submit Q1 Answer"):
        if "B) Lipophilicity" in q1:
            st.success("Correct! LogP measures the partition coefficient, indicating how lipophilic a compound is.")
        else:
            st.error("Incorrect. Remember that LogP relates directly to partition and lipophilicity.")

    st.markdown("---")

    # Question 2
    st.markdown("#### Question 2: Which SMILES string correctly represents Ethanol?")
    q2 = st.radio(
        "Select your answer for Q2:",
        ["A) H2O", "B) CCO", "C) c1ccccc1", "D) CC(=O)O"],
        key="q2_ans"
    )
    if st.button("Submit Q2 Answer"):
        if "CCO" in q2:
            st.success("Correct! `CCO` is the standard SMILES notation for Ethanol.")
        else:
            st.error("Incorrect. `CCO` stands for a two-carbon chain attached to an alcohol group.")

    st.markdown("---")

    # Question 3
    st.markdown("#### Question 3: What does TPSA (Topological Polar Surface Area) primarily help estimate?")
    q3 = st.radio(
        "Select your answer for Q3:",
        [
            "A) Membrane permeability and drug absorption potential",
            "B) Total number of aliphatic carbon atoms",
            "C) Optical rotation angle in a polarimeter",
            "D) Melting point of the crystalline lattice"
        ],
        key="q3_ans"
    )
    if st.button("Submit Q3 Answer"):
        if "A) Membrane permeability" in q3:
            st.success("Correct! TPSA is widely used as a predictor of drug transport properties like intestinal absorption.")
        else:
            st.error("Incorrect. TPSA measures surface contributions of polar atoms (O and N).")

    st.markdown("---")

    # Question 4
    st.markdown("#### Question 4: Under Lipinski's Rule of Five, what is the maximum recommended molecular weight threshold?")
    q4 = st.radio(
        "Select your answer for Q4:",
        [
            "A) 200 Daltons",
            "B) 500 Daltons",
            "C) 1000 Daltons",
            "D) 2000 Daltons"
        ],
        key="q4_ans"
    )
    if st.button("Submit Q4 Answer"):
        if "B) 500 Daltons" in q4:
            st.success("Correct! Lipinski's rule states that most orally active drugs have a molecular weight under 500 Daltons.")
        else:
            st.error("Incorrect. The threshold for molecular weight in Lipinski's Rule of Five is 500 Daltons.")

    st.markdown("---")

    # Question 5
    st.markdown("#### Question 5: In chemical bonding analysis, how are Pi ($\pi$) bonds typically formed?")
    q5 = st.radio(
        "Select your answer for Q5:",
        [
            "A) Head-on overlap of atomic orbitals",
            "B) Lateral (side-by-side) overlap of unhybridized p-orbitals",
            "C) Electrostatic attraction between ions",
            "D) Sharing of a single electron pair"
        ],
        key="q5_ans"
    )
    if st.button("Submit Q5 Answer"):
        if "B) Lateral" in q5:
            st.success("Correct! Pi bonds arise from the lateral overlap of parallel p-orbitals, common in double and triple bonds.")
        else:
            st.error("Incorrect. Head-on overlap forms sigma ($\sigma$) bonds, whereas lateral overlap forms pi ($\pi$) bonds.")

elif page == "🗫 Final Project Report":
    st.title("🗫 CHEM𖦹CHAOS Final Project Report")
    st.markdown(
        "### Abstract & Overview\n"
        "**CHEM𖦹CHAOS** is an integrated virtual cheminformatics laboratory designed to bridge theoretical chemistry "
        "and computational analysis. By leveraging RDKit, Py3Dmol, and ReportLab, students can interactively explore molecular structures, "
        "calculate physicochemical descriptors, analyze bonding configurations, and generate professional PDF documentation."
    )
    st.markdown(
        "### Key Modules Included:\n"
        "1. **Home & Theory:** Foundational guidance on SMILES notation, descriptors, and Lipinski's rules.\n"
        "2. **Molecular Visualization Studio:** Real-time 2D rendering and interactive 3D conformational analysis.\n"
        "3. **Properties & Report Generator:** Comprehensive property calculations with automated PDF export.\n"
        "4. **Structure-Property Analysis:** Batch comparative analysis across multiple chemical entities.\n"
        "5. **Interactive Assessment:** 5-question evaluation module to test student comprehension."
    )
    st.success("All application modules are fully operational and verified.")
