import streamlit as st
from Bio.PDB import PDBParser, PPBuilder
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO
import requests

st.set_page_config(page_title="Ramachandran Plot App", layout="wide")
st.title("📊 Ramachandran Plot Generator (Mobile-Friendly)")

# -----------------------------
# PDB ID input
# -----------------------------
pdb_id = st.text_input("Enter PDB ID (e.g., 1UBQ)").upper()
chain_id = st.text_input("Enter chain ID (default = A)", "A")

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

        allowed_regions = [
            {"phi": (-160, -40), "psi": (-80, 50)},  # Beta sheet
            {"phi": (-90, -30), "psi": (-70, 10)},  # Right-handed helix
        ]

        for reg in allowed_regions:
            ax.add_patch(
                plt.Rectangle(
                    (reg["phi"][0], reg["psi"][0]),
                    reg["phi"][1]-reg["phi"][0],
                    reg["psi"][1]-reg["psi"][0],
                    fill=False, edgecolor="red", linewidth=2, linestyle="--"
                )
            )

        ax.legend()
        st.pyplot(fig)

        # -----------------------------
        # Stats
        # -----------------------------
        allowed_count = 0
        for ph, ps in zip(phi, psi):
            for reg in allowed_regions:
                if reg["phi"][0] <= ph <= reg["phi"][1] and reg["psi"][0] <= ps <= reg["psi"][1]:
                    allowed_count += 1
                    break

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
