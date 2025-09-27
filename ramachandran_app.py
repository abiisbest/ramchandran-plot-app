import streamlit as st
from Bio.PDB import PDBParser, PPBuilder
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO
import requests

st.set_page_config(page_title="Ramachandran Plot App", layout="wide")
st.title("📊 Ramachandran Plot Generator (High Accuracy)")

# -----------------------------
# PDB ID input
# -----------------------------
pdb_id = st.text_input("Enter PDB ID (e.g., 1UBQ)").upper()
chain_id = st.text_input("Enter chain ID (default = A)", "A")

# -----------------------------
# Residue-specific allowed regions (simplified)
# -----------------------------
# Format: (phi_min, phi_max, psi_min, psi_max)
allowed_regions = {
    "ALL": [(-180, 0, -180, 180)],  # Generic fallback
    "GLY": [(-180, 0, -180, 180)],  # Glycine flexible
    "PRO": [(-80, -40, 130, 180)],  # Proline restricted
    # Right-handed helix typical
    "ALA": [(-160, -30, -70, -5)],
    "VAL": [(-160, -30, -70, -5)],
    "LEU": [(-160, -30, -70, -5)],
    "ILE": [(-160, -30, -70, -5)],
    "MET": [(-160, -30, -70, -5)],
    "CYS": [(-160, -30, -70, -5)],
    "SER": [(-160, -30, -70, -5)],
    "THR": [(-160, -30, -70, -5)],
    "TRP": [(-160, -30, -70, -5)],
    "PHE": [(-160, -30, -70, -5)],
    "TYR": [(-160, -30, -70, -5)],
    "HIS": [(-160, -30, -70, -5)],
    "ASN": [(-160, -30, -70, -5)],
    "GLN": [(-160, -30, -70, -5)],
    "ASP": [(-160, -30, -70, -5)],
    "GLU": [(-160, -30, -70, -5)],
    "LYS": [(-160, -30, -70, -5)],
    "ARG": [(-160, -30, -70, -5)],
}

def is_allowed(residue, phi, psi):
    residue = residue.upper()
    regions = allowed_regions.get(residue, allowed_regions["ALL"])
    for phi_min, phi_max, psi_min, psi_max in regions:
        if phi_min <= phi <= phi_max and psi_min <= psi <= psi_max:
            return True
    return False

# -----------------------------
# Ramachandran function
# -----------------------------
def ramachandran_plot(pdb_file, chain_id="A"):
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        model = structure[0]

        if chain_id not in model:
            st.error(f"Chain '{chain_id}' not found.")
            return

        chain = model[chain_id]
        ppb = PPBuilder()
        phi_psi = []
        residues_list = []

        for pp in ppb.build_peptides(chain):
            phi_psi.extend(pp.get_phi_psi_list())
            residues_list.extend([res.get_resname() for res in pp])

        phi, psi, residues = [], [], []
        for i, (ph, ps) in enumerate(phi_psi):
            if ph is not None and ps is not None:
                phi.append(np.degrees(ph))
                psi.append(np.degrees(ps))
                residues.append(residues_list[i] if i < len(residues_list) else "UNK")

        # -----------------------------
        # Plot
        # -----------------------------
        fig, ax = plt.subplots(figsize=(6,6))
        ax.scatter(phi, psi, c="blue", s=25, alpha=0.6, label="Residues")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xlabel("Phi (φ)")
        ax.set_ylabel("Psi (ψ)")
        ax.set_title(f"Ramachandran Plot: {pdb_id} Chain {chain_id}")

        ax.axhline(0, color='gray', linewidth=0.5)
        ax.axvline(0, color='gray', linewidth=0.5)

        st.pyplot(fig)

        # -----------------------------
        # Stats with residue-specific regions
        # -----------------------------
        allowed_count = sum(is_allowed(res, ph, ps) for res, ph, ps in zip(residues, phi, psi))

        st.write(f"📊 Total residues checked: {len(phi)}")
        st.write(f"✅ Allowed residues: {allowed_count} ({allowed_count/len(phi)*100:.2f}%)")
        st.write(f"❌ Disallowed residues: {len(phi)-allowed_count}")

        # -----------------------------
        # Table & CSV
        # -----------------------------
        df = pd.DataFrame({
            "Residue": residues,
            "Phi (°)": phi,
            "Psi (°)": psi
        })
        st.dataframe(df)

        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Phi/Psi angles as CSV",
            data=csv,
            file_name=f'{pdb_id}_phi_psi.csv',
            mime='text/csv'
        )

    except Exception as e:
        st.error(f"Error processing PDB file: {e}")

# -----------------------------
# Fetch PDB and run
# -----------------------------
if pdb_id:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_file = StringIO(response.text)
        ramachandran_plot(pdb_file, chain_id)
    else:
        st.error("PDB ID not found!")
