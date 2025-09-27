import streamlit as st
from Bio.PDB import PDBParser, PPBuilder
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO, BytesIO
import requests

st.set_page_config(page_title="Ramachandran Plot App", layout="wide")
st.title("📊 Ramachandran Plot Generator (High Accuracy)")

# -----------------------------
# INPUT SECTION
# -----------------------------

st.subheader("1. Input Protein Structure")

# PDB ID input
col_id, col_chain = st.columns([3, 1])
pdb_id = col_id.text_input("Enter PDB ID (e.g., 1UBQ) or leave blank to upload file").upper()
chain_id = col_chain.text_input("Enter Chain ID (default = A)", "A")

# File upload input
uploaded_file = st.file_uploader("OR Upload a PDB file (.pdb)", type=["pdb"])

# ----------------------------------------------------
# 1. ACCURATE RAMACHANDRAN REGIONS (CORE & ALLOWED)
#    These boundaries are approximated based on statistical data (e.g., MolProbity/PROCHECK)
#    Each region is defined by (phi_min, phi_max, psi_min, psi_max)
# ----------------------------------------------------

# CORE (Favored) Regions for non-Glycine, non-Proline residues
CORE_REGIONS = [
    # Right-handed Alpha Helix (R_alpha)
    (-100, -35, -65, -20),
    # Parallel/Antiparallel Beta Sheet (P_beta/A_beta)
    (-140, -90, 135, 180),
    (-140, -90, -180, -160),
    (-140, -90, 100, 125)
]

# ALLOWED (Outer boundary) Regions for non-Glycine, non-Proline residues
ALLOWED_REGIONS = [
    # Alpha Helix region (wider)
    (-140, -20, -100, 15),
    # Beta Sheet region (wider)
    (-180, -60, 90, 180),
    (-180, -60, -180, -150),
    # Left-handed Alpha Helix (L_alpha)
    (35, 100, 20, 80)
]

# RESIDUE-SPECIFIC OVERRIDES
def get_regions(residue):
    """Returns (core_regions, allowed_regions) for a given residue name."""
    res = residue.upper()
    
    if res == "GLY":
        # Glycine is symmetrical and much more flexible
        gly_core = [
            (-180, -60, -180, 180),
            (60, 180, -180, 180)
        ]
        gly_allowed = [
            (-180, 180, -180, 180)
        ]
        return gly_core, gly_allowed
    
    if res == "PRO":
        # Proline is highly restricted due to cyclic side chain
        pro_core = [(-80, -35, -30, 30)]
        pro_allowed = [(-90, -30, -50, 50)]
        return pro_core, pro_allowed
    
    # Default for all other standard residues
    return CORE_REGIONS, ALLOWED_REGIONS

def is_allowed(residue, phi, psi):
    """Checks if a phi/psi pair falls within the 'Allowed' region boundary."""
    core_regions, allowed_regions = get_regions(residue)
    
    # 1. Check against the tighter CORE regions first
    for phi_min, phi_max, psi_min, psi_max in core_regions:
         if phi_min <= phi <= phi_max and psi_min <= psi <= psi_max:
            return "Core"

    # 2. If not in Core, check against the wider ALLOWED regions
    for phi_min, phi_max, psi_min, psi_max in allowed_regions:
        if phi_min <= phi <= phi_max and psi_min <= psi <= psi_max:
            return "Allowed"

    # 3. Otherwise, it's an Outlier
    return "Disallowed"

# -----------------------------
# PLOTTING HELPER FUNCTION
# -----------------------------

def plot_regions(ax, regions, color, label):
    """Draws rectangular Ramachandran regions on the matplotlib axis."""
    for phi_min, phi_max, psi_min, psi_max in regions:
        ax.fill_between(
            [phi_min, phi_max], 
            psi_min, 
            psi_max, 
            color=color, 
            alpha=0.2, 
            zorder=0,
            label=label if label else ""
        )
    if label:
        # Clear label after first drawing to avoid duplicate legends
        label = None
    return label


# -----------------------------
# Ramachandran function
# -----------------------------
def ramachandran_plot(pdb_file, chain_id="A", source_name="PDB"):
    st.markdown("---")
    st.subheader(f"Results for Structure: **{source_name}**, Chain: **{chain_id}**")
    
    try:
        # Note: BioPython PDBParser can handle file-like objects (StringIO/BytesIO)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        model = structure[0]

        if chain_id not in model:
            st.error(f"Chain '{chain_id}' not found in structure '{source_name}'. Please check the PDB entry for valid chains.")
            return

        chain = model[chain_id]
        ppb = PPBuilder()
        
        # Extract Phi/Psi and Residue Names
        phi_psi = []
        residues_list = []
        
        for pp in ppb.build_peptides(chain):
            # Iterate through residues to align names with angles
            for i, res in enumerate(pp):
                # Use phi/psi from the builder's peptide list
                ph, ps = pp.get_phi_psi_list()[i] 
                
                # Check for None (usually N-term, C-term, or proline before Pro)
                if ph is not None and ps is not None:
                    phi_psi.append((np.degrees(ph), np.degrees(ps)))
                    residues_list.append(res.get_resname())
        
        phi, psi = zip(*phi_psi) if phi_psi else ([], [])
        residues = residues_list
        
        if not phi:
            st.warning(f"No valid phi/psi angles found for chain '{chain_id}' in '{source_name}'. It might be too short or broken.")
            return
            
        # -----------------------------
        # Calculate Stats and Categorize Points
        # -----------------------------
        
        categories = [is_allowed(res, ph, ps) for res, ph, ps in zip(residues, phi, psi)]
        
        core_count = categories.count("Core")
        allowed_count = categories.count("Allowed") + core_count # Core points are also allowed
        disallowed_count = categories.count("Disallowed")
        total_count = len(phi)
        
        if total_count == 0:
            st.warning("No residues analyzed.")
            return
            
        total_allowed_percent = allowed_count / total_count * 100

        # -----------------------------
        # Plotting the Regions and Points
        # -----------------------------
        fig, ax = plt.subplots(figsize=(8,8))

        # 1. Draw Background Regions (Uses ALA/Standard regions for visual guide)
        std_core, std_allowed = get_regions("ALA")
        
        # Plot Allowed Region (Light Blue)
        plot_regions(ax, std_allowed, 'skyblue', 'Allowed')

        # Plot Core/Favored Region (Darker Blue)
        plot_regions(ax, std_core, 'dodgerblue', 'Favored (Core)')
        
        # 2. Plot the Actual Phi/Psi Points
        
        # Group points for colored plotting
        data_df = pd.DataFrame({'phi': phi, 'psi': psi, 'category': categories})
        
        # Define colors for better visualization
        color_map = {
            "Core": "black",
            "Allowed": "darkorange",
            "Disallowed": "red"
        }
        
        # Plot each category separately to control z-order and legend
        for category, color in color_map.items():
            subset = data_df[data_df['category'] == category]
            if not subset.empty:
                 ax.scatter(
                    subset['phi'], 
                    subset['psi'], 
                    c=color, 
                    s=20, 
                    alpha=0.7, 
                    zorder=5, 
                    label=f'{category} Points'
                )


        # 3. Final Plot Settings
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

        # -----------------------------
        # Display Stats
        # -----------------------------
        st.subheader("Analysis Summary")
        col1, col2, col3 = st.columns(3)
        
        col1.metric("Total Residues Checked", total_count)
        
        # Display the Total Allowed Percentage (Core + Allowed)
        col2.metric(
            "✅ Total Allowed Percentage", 
            f"{total_allowed_percent:.2f}%", 
            f"{allowed_count} points (Core + Allowed)"
        )
        
        # Display Outlier/Disallowed
        col3.metric(
            "⚠️ Outlier/Disallowed", 
            f"{disallowed_count}", 
            f"({disallowed_count/total_count*100:.2f}%) points"
        )
        
        # -----------------------------
        # Table & CSV
        # -----------------------------
        st.subheader("Phi/Psi Data Table")
        
        df = pd.DataFrame({
            "Residue": [f"{res} ({i+1})" for i, res in enumerate(residues)], # Add index for clarity
            "Phi (°)": phi,
            "Psi (°)": psi,
            "Classification": categories
        })
        st.dataframe(df, use_container_width=True)

        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Phi/Psi angles as CSV",
            data=csv,
            file_name=f'{source_name}_{chain_id}_phi_psi.csv',
            mime='text/csv'
        )

    except Exception as e:
        st.error(f"An error occurred during plot generation: {e}")

# -----------------------------
# Execution Logic (Fetch PDB or use uploaded file)
# -----------------------------

if uploaded_file is not None:
    # 1. Handle File Upload
    st.info(f"Using uploaded file: {uploaded_file.name}")
    try:
        # Read the uploaded file contents as bytes
        bytes_data = uploaded_file.getvalue()
        # Decode bytes to string, then wrap in StringIO for PDBParser
        string_data = bytes_data.decode("utf-8")
        pdb_file = StringIO(string_data)
        
        # Determine source name for display and file saving
        file_name = uploaded_file.name.split('.')[0]
        ramachandran_plot(pdb_file, chain_id, source_name=file_name)
    except Exception as e:
        st.error(f"Error reading uploaded file: {e}")

elif pdb_id:
    # 2. Handle PDB ID Fetch
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
