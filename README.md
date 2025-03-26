# X-SPAN---Crosslink-Structural-Pattern-Analyzer-

```

   _  __     _____ ____  ___    _   __
  | |/ /    / ___// __ \/   |  / | / /
  |   /_____\__ \/ /_/ / /| | /  |/ / 
 /   /_____/__/ / ____/ ___ |/ /|  /  
/_/|_|    /____/_/   /_/  |_/_/ |_/   

```

**Version:** 1.0.0  
**Authors:** A. Vetrano, A. Di Ianni, N. Di Fonte, G. Dell'Orletta, S. Reale, I. Daidone, and C. Iacobucci 
**Language:** Python 3  G
**Interface:** PyQt5 GUI

---

## 📌 Description

**X-SPAN** is an interactive tool with a graphical interface that enables automated, integrated analysis of **crosslinking mass spectrometry data**. It supports `.zhrm`, `.xlsx`, and `.mzid` files, downloads AlphaFold structures, extracts structural information (secondary structure and pLDDT), and identifies **continuous structural motifs** (helices, sheets, loops) between crosslinked regions. It outputs a series of processed CSV files and plots.


> ❗ **Why only intra-protein cross-links within 20 amino acids?**  
This threshold is based on structural insights presented in the X-SPAN publication.
Short-range cross-links (< 20 residues apart) are often dismissed as trivial, but our work shows they contain 
valuable information about local secondary structures. See the associated publication for further details.


---

## 🧰 Features

- ✅ User-friendly PyQt5 GUI
- ✅ Supports **.zhrm**, **.xlsx**, and **.mzid** inputs
- ✅ Automatic conversion of `.zhrm` to `.xlsx`
- ✅ AlphaFold **PDB downloader** for UniProt IDs
- ✅ Structural parsing via **PyMOL** (pymol2)
- ✅ Identification of continuous structural motifs
- ✅ Output as structured `.csv` files
- ✅ Automatic generation of **plots** and **frequency summaries**

---

## 🗂️ Requirements

- Python ≥ 3.7  
- PyQt5  
- pandas  
- matplotlib  
- tqdm  
- requests  
- pymol2  
- openpyxl  
- lxml (for `.mzid` parsing)

Install all dependencies with:

```bash
pip install -r requirements.txt
```

> ⚠️ You need a **PyMOL installation with pymol2 support** (e.g., via conda).

---

## 🚀 Usage

1. **Launch the script**:
   ```bash
   python xspan.py
   ```

2. **Select paths** in the GUI:
   - Input file (`.zhrm`, `.xlsx`, or `.mzid`)
   - Output folder
   - Output folder name

3. The script will:
   - Convert the input if needed
   - Normalize and clean the Excel dataset
   - Download AlphaFold PDBs based on UniProt IDs
   - Analyze secondary structure and pLDDT using PyMOL
   - Annotate crosslinks with structural motifs
   - Export:
     - Annotated CSVs
     - Frequency tables
     - PNG plots

---

## 📄 Expected .xlsx Format

| PDB_ID (UniProt) | Residue1 | Residue2 |
|------------------|----------|----------|
| Q9Y6N9           | 42       | 58       |

---

## 📦 Output Files

Inside the output folder (`your_output_dir/`), you'll find:

- `*_proteins.txt`: UniProt IDs list
- `*_XL.xlsx`: Cleaned dataset
- `output.txt`: Secondary structure and pLDDT data
- `*_results.csv`: Final annotated results
- `filtered*.csv`: Filtered motif-based crosslinks
- `counts.csv`, `counts_updated.csv`: Distance distributions and frequencies
- PNG plots: One for each motif type and frequency

---

## 🛑 Limitations

- Only **intra-protein** crosslinks are analyzed
- Distance threshold between residues is < 20
- Requires internet for downloading AlphaFold structures
- Requires PyMOL with pymol2 support

---

## 📬 Contact

For support or collaboration, contact: **claudio.iacobucci@univaq.it**
