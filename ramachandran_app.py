
import streamlit as st
from Bio.PDB import PDBParser, PPBuilder
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO, BytesIO
import requests

st.set_page_config(page_title="Ramachandran Plot App", layout="wide")
st.title("📊 Ramachandran Plot Generator by Abijeet and Yathish")

# --- [YOUR CONSTANTS AND HELPER FUNCTIONS - UNCHANGED] ---

# CORE (Favored) Regions for non-Glycine, non-Proline residues
CORE_REGIONS = [
    (-100, -35, -65, -20),
    (-140, -90, 135, 180),
    (-140, -90, -180, -160),
    (-140, -90, 100, 125)
]

# ALLOWED (Outer boundary) Regions for non-Glycine, non-Proline residues
ALLOWED_REGIONS = [
    (-140, -20, -100, 15),
    (-180, -60, 90, 180),
    (-180, -60, -180, -150),
    (35, 100, 20, 80)
]

def get_regions(residue):
    res = residue.upper()
    if res == "GLY":
        gly_core = [(-180, -60, -180, 180), (60, 180, -180, 180)]
        gly_allowed = [(-180, 180, -180, 180)]
        return gly_core, gly_allowed
    if res == "PRO":
        pro_core = [(-80, -35, -30, 30)]
        pro_allowed = [(-90, -30, -50, 50)]
        return pro_core, pro_allowed
    return CORE_REGIONS, ALLOWED_REGIONS

def is_allowed(residue, phi, psi):
    core_regions, allowed_regions = get_regions(residue)
    for phi_min, phi_max, psi_min, psi_max in core_regions:
        if phi_min <= phi <= phi_max and psi_min <= psi <= psi_max:
            return "Core"
    for phi_min, phi_max, psi_min, psi_max in allowed_regions:
        if phi_min <= phi <= phi_max and psi_min <= psi <= psi_max:
            return "Allowed"
    return "Disallowed"

def plot_regions(ax, regions, color, label):
    first_segment = True
    original_label = label
    for phi_min, phi_max, psi_min, psi_max in regions:
        current_label = original_label if first_segment else ""
        ax.fill_between([phi_min, phi_max], psi_min, psi_max, color=color, alpha=0.2, zorder=0, label=current_label)
        first_segment = False
    return label

# Ramachandran Plotting Function
def ramachandran_plot(pdb_file, chain_id="A", source_name="PDB"):
    st.markdown("---")
    st.subheader(f"Results for Structure: **{source_name}**, Chain: **{chain_id}**")
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        model = structure[0]
        if chain_id not in model:
            st.error(f"Chain '{chain_id}' not found in structure '{source_name}'.")
            return
        chain = model[chain_id]
        ppb = PPBuilder()
        phi_psi = []
        residues_list = []
        for pp in ppb.build_peptides(chain):
            for i, res in enumerate(pp):
                ph, ps = pp.get_phi_psi_list()[i]
                if ph is not None and ps is not None:
                    phi_psi.append((np.degrees(ph), np.degrees(ps)))
                    residues_list.append(res.get_resname())
        phi, psi = zip(*phi_psi) if phi_psi else ([], [])
        residues = residues_list
        if not phi:
            st.warning(f"No valid phi/psi angles found for chain '{chain_id}'.")
            return
            
        categories = [is_allowed(res, ph, ps) for res, ph, ps in zip(residues, phi, psi)]
        core_count = categories.count("Core")
        allowed_count = categories.count("Allowed") + core_count
        disallowed_count = categories.count("Disallowed")
        total_count = len(phi)
        
        if total_count == 0:
            st.warning("No residues analyzed.")
            return
            
        total_allowed_percent = allowed_count / total_count * 100
        disallowed_percent = disallowed_count / total_count * 100

        fig, ax = plt.subplots(figsize=(8,8))
        std_core, std_allowed = get_regions("ALA")
        plot_regions(ax, std_allowed, 'skyblue', 'Allowed')
        plot_regions(ax, std_core, 'dodgerblue', 'Favored (Core)')
        data_df = pd.DataFrame({'phi': phi, 'psi': psi, 'category': categories})
        color_map = {"Core": "black", "Allowed": "darkorange", "Disallowed": "red"}
        for category, color in color_map.items():
            subset = data_df[data_df['category'] == category]
            if not subset.empty:
                ax.scatter(subset['phi'], subset['psi'], c=color, s=20, alpha=0.7, zorder=5, label=f'{category} Points')
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_yticks(np.arange(-180, 181, 60))
        ax.set_xlabel("Phi (φ) [degrees]")
        ax.set_ylabel("Psi (ψ) [degrees]")
        ax.set_title(f"Ramachandran Plot: {source_name} Chain {chain_id}", fontsize=14)
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.axhline(0, color='gray', linewidth=0.5)
        ax.axvline(0, color='gray', linewidth=0.5)
        ax.legend(loc='lower left', frameon=True)
        st.pyplot(fig)
        
        st.subheader("Phi (φ) and Psi (ψ) Dihedral Angles Along Sequence")
        sequence_index = np.arange(1, total_count + 1)
        fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        ax1.plot(sequence_index, phi, marker='o', linestyle='-', markersize=3, color='dodgerblue', linewidth=1, label="Phi (φ)")
        ax1.set_ylabel("Phi (φ) [degrees]")
        ax1.set_title(f"Dihedral Angle Sequence for {source_name} Chain {chain_id}")
        ax1.set_ylim(-180, 180)
        ax1.axhline(0, color='gray', linestyle='--', linewidth=0.5)
        ax1.grid(True, linestyle=':', alpha=0.6)
        ax1.legend(loc='upper right')
        ax2.plot(sequence_index, psi, marker='o', linestyle='-', markersize=3, color='darkorange', linewidth=1, label="Psi (ψ)")
        ax2.set_ylabel("Psi (ψ) [degrees]")
        ax2.set_xlabel("Residue Index")
        ax2.set_ylim(-180, 180)
        ax2.axhline(0, color='gray', linestyle='--', linewidth=0.5)
        ax2.grid(True, linestyle=':', alpha=0.6)
        ax2.legend(loc='upper right')
        step = max(1, total_count // 10)
        ax2.set_xticks(np.arange(1, total_count + 1, step))
        plt.tight_layout()
        st.pyplot(fig2)

        st.subheader("Analysis Summary")
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Residues Checked", total_count)
        col2.metric("✅ Total Allowed Percentage", f"{total_allowed_percent:.2f}%", f"{allowed_count} points (Core + Allowed)")
        col3.metric("⚠️ Outlier/Disallowed", f"{disallowed_percent:.2f}%", f"{disallowed_count} points")
        
        st.subheader("Phi/Psi Data Table")
        df = pd.DataFrame({"Residue": [f"{res} ({i+1})" for i, res in enumerate(residues)], "Phi (°)": phi, "Psi (°)": psi, "Classification": categories})
        st.dataframe(df, use_container_width=True)
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(label="📥 Download Phi/Psi angles as CSV", data=csv, file_name=f'{source_name}_{chain_id}_phi_psi.csv', mime='text/csv')

    except Exception as e:
        st.error(f"An error occurred during plot generation: {e}")

# -----------------------------
# Execution Logic (Fetch PDB or use uploaded file)
# -----------------------------
st.subheader("1. Input Protein Structure")

# Create a form to wrap the inputs and the submit button
with st.form(key='my_form'):
    # PDB ID input
    col_id, col_chain = st.columns([3, 1])
    pdb_id = col_id.text_input("Enter PDB ID (e.g., 1UBQ) or leave blank to upload file").upper()
    chain_id = col_chain.text_input("Enter Chain ID (default = A)", "A")

    # File upload input
    uploaded_file = st.file_uploader("OR Upload a PDB file (.pdb)", type=["pdb"])

    # Create the submit button inside the form
    submit_button = st.form_submit_button(label='Generate Plot')

# Check if the submit button was pressed
if submit_button:
    # This entire block will only execute if the button is clicked
    if uploaded_file is not None:
        st.info(f"Using uploaded file: {uploaded_file.name}")
        try:
            bytes_data = uploaded_file.getvalue()
            string_data = bytes_data.decode("utf-8")
            pdb_file = StringIO(string_data)
            
            file_name = uploaded_file.name.split('.')[0]
            ramachandran_plot(pdb_file, chain_id, source_name=file_name)
        except Exception as e:
            st.error(f"Error reading uploaded file: {e}")
            
    elif pdb_id:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        st.info(f"Fetching PDB file from RCSB: {url}...")
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                pdb_file = StringIO(response.text)
                ramachandran_plot(pdb_file, chain_id, source_name=pdb_id)
            else:
                st.error(f"PDB ID '{pdb_id}' not found or inaccessible. Status code: {response.status_code}")
        except requests.exceptions.RequestException as e:
            st.error(f"Network error while fetching PDB ID: {e}")
    else:
        st.warning("Please enter a PDB ID or upload a PDB file to generate the plot.")
