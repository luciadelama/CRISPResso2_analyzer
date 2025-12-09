import os
import re
from sys import prefix
import zipfile
import tempfile
import glob
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from modules import parser
from typing import Optional, Tuple

def extract_zip_to_tmp(zip_path: str):
    """Extracts the ZIP to a temporary directory and returns the TemporaryDirectory object."""
    tmp_dir = tempfile.TemporaryDirectory()     
    with zipfile.ZipFile(zip_path, 'r') as z:
        z.extractall(tmp_dir.name)
    return tmp_dir

def parse_fastq_metadata(name: str) -> Optional[Tuple[str, str]]:
    m = re.match(
        r"^Day(?P<day>\d+)_"
        r"(?:(?P<gen>[^_]+)_)?"
        r"(?P<variant>[^_]+)_"
        r"(?P<treat>[^_]+)_"
        r"(?P<rep>\d+)"
        r"(?P<rerun>_rerun)?"
        r".*\.fastq(?:\.gz)?$", 
        name, 
        re.IGNORECASE
        )
    if not m:
        return None
    day = m.group("day")
    rep = m.group("rep")
    treat = (m.group("treat") or "untreated").lower()
    rerun = m.group("rerun")
    return f"Replica-{rep}", f"Day{day}_{treat}_rerun" if rerun else f"Day{day}_{treat}"

def reorganize_extracted_fastqs(tmp_path: Path):
    for p in Path(tmp_path).rglob("*.fastq.gz"):
        meta = parse_fastq_metadata(p.name)
        if meta is None:
            dest = Path(tmp_path) / "Replica-unknown" / "Day-unknown"
        else:
            rep, treat = meta
            dest = Path(tmp_path) / rep / treat
        dest.mkdir(parents=True, exist_ok=True)
        newp = dest / p.name
        if p.resolve() != newp.resolve():
            try:
                p.rename(newp)
            except Exception:
                pass

# def ensure_tmpdir(prefix="crispr_tmp_"):
#     """Creates and returns a TemporaryDirectory with a given prefix."""
#     d = tempfile.TemporaryDirectory(prefix=prefix)
#     # Keep a strong reference by attaching attribute so it doesn't GC; caller should hold it
#     return d

def list_replicates_from_outputs(outputs_root: Path):
    root = Path(outputs_root)
    reps = sorted([p.name for p in root.iterdir() if p.is_dir() and p.name.startswith("Replica-")])
    if reps:
        return reps
    return sorted([p.name for p in root.iterdir() if p.is_dir()])

def list_treatments_from_outputs(outputs_root: Path):
    root = Path(outputs_root)
    reps = [p for p in root.iterdir() if p.is_dir() and p.name.startswith("Replica-")]
    names = set()
    for r in reps:
        for t in r.iterdir():
            if t.is_dir():
                names.add(t.name)
    if names:
        return sorted(names)
    names = set(p.name for p in root.iterdir() if p.is_dir())
    return sorted(names)

def process_treatment_folder(treatment_path, wt_seq, mut_seq, amplicon_seq):
    """Process a single CRISPResso output folder and extracts key editing metrics.

    CRISPResso may output files directly under `treatment_path` or within a nested
    directory like `CRISPResso_on_*`. This function adapts to both layouts.
    """

    # Determine actual directory containing CRISPResso outputs
    actual_dir = treatment_path
    direct_log = os.path.join(actual_dir, "CRISPResso_RUNNING_LOG.txt")
    if not os.path.exists(direct_log):
        found = None
        try:
            found = next((p.parent for p in Path(treatment_path).rglob("CRISPResso_RUNNING_LOG.txt")), None)
        except Exception:
            found = None
        if not found:
            return None
        actual_dir = str(found)

    # Define the expected CRISPResso output files paths (under actual_dir)
    map_stat_path = os.path.join(actual_dir, "CRISPResso_mapping_statistics.txt")
    quant_edit_path = os.path.join(actual_dir, "CRISPResso_quantification_of_editing_frequency.txt")
    frameshift_path = os.path.join(actual_dir, "Frameshift_analysis.txt")
    # amplicon_path = os.path.join(actual_dir, "CRISPResso_RUNNING_LOG.txt")

    # Locate HTML report path under the sample folder
    sample_dir = actual_dir
    base_name = os.path.basename(actual_dir)
    if base_name.startswith("CRISPResso_on_"):
        sample_dir = os.path.dirname(actual_dir)
    html_matches = sorted(glob.glob(os.path.join(sample_dir, "CRISPResso_on_*.html")))
    report_html = html_matches[0] if html_matches else None

    # Look for the allele frequency table (filename may vary by sgRNA sequence)
    allele_candidates = [
        os.path.join(actual_dir, "Alleles_frequency_table_around_sgRNA_*.txt"),
        os.path.join(actual_dir, "Alleles_frequency_table_around_cut_site_*.txt"),
        os.path.join(actual_dir, "Alleles_frequency_table*.txt"),
    ]
    allele_matches = []
    for pat in allele_candidates:
        allele_matches.extend(glob.glob(pat))
    if not allele_matches:
        return None  # Allele file not found
    alleles_path = sorted(allele_matches)[0]

    # Verify that all required files exist before processing
    required_paths = [map_stat_path, quant_edit_path, frameshift_path, alleles_path]
    if not all(os.path.exists(p) for p in required_paths):
        return None

    # Use the functions to extract the relevant data
    total_reads = parser._get_total_reads(map_stat_path)
    modif_data = parser._get_modif_reads(quant_edit_path, total_reads)
    frame_data = parser._get_frame_reads(frameshift_path, total_reads)
    # amplicon_data = parser._get_input_run(amplicon_path)
    mut_wt_data = parser._get_mut_wt_reads(
        alleles_path, 
        wt_seq, 
        mut_seq, 
        amplicon_seq, 
        total_reads
    )

    # Return all the processed information in a dictionary
    return {
        "Treatment": os.path.basename(treatment_path),
        "WT/Unmodified%": modif_data["wt_unmodified_percent"],
        "Reads (WT/Unmodified)": modif_data["unmodified_reads"],
        "MUT%": mut_wt_data["mut_percent"],
        "Reads (MUT)": mut_wt_data["mut_reads"],
        "WT*%": mut_wt_data["wt_percent"],
        "Reads (WT*)": mut_wt_data["wt_reads"],
        "Frameshift%": frame_data["frameshift_percent"],
        "Reads (Frameshift)": frame_data["frameshift_reads"],
        "In-frame%": frame_data["inframe_percent"],
        "Reads (In-frame)": frame_data["inframe_reads"],
        "Indel%": modif_data["indel_percent"],
        "Reads (Modified)": modif_data["modified_reads"],
        "MUT/WT*%": mut_wt_data["mut_wt_percent"],
        "Total Reads": total_reads,
        "ReportPath": report_html,
    }

def calculate_sensitivity(df, column="MUT/WT*%", reference_day=None, reference_treatment=None):
    if df.empty or column not in df.columns:
        return pd.Series(dtype=float)
    if reference_treatment is not None:
        sel = df[df["Treatment"] == reference_treatment]
        if sel.empty:
            return pd.Series([None] * len(df))
        base_value = float(sel[column].iloc[0])
    elif reference_day is None:
        base_value = float(df[column].iloc[0])
    else:
        pref = f"Day{int(reference_day)}_"
        sel = df[df["Treatment"].astype(str).str.startswith(pref, na=False)]
        if "Day{}_untreated".format(int(reference_day)) in df["Treatment"].values:
            sel = df[df["Treatment"] == "Day{}_untreated".format(int(reference_day))]
        if sel.empty:
            return pd.Series([None] * len(df))
        base_value = float(sel[column].iloc[0])
    if base_value == 0:
        return pd.Series([None] * len(df))
    return round(df[column].astype(float) / base_value * 100, 2)


# def find_root_folder(tmp_dir_path: str):
    """Devuelve la primera carpeta dentro del tmp dir (raíz del ZIP)."""
    roots = [
        f for f in os.listdir(tmp_dir_path)
        if os.path.isdir(os.path.join(tmp_dir_path, f))
    ]
    if not roots:
        return None
    return os.path.join(tmp_dir_path, roots[0])


def build_dfs_by_replicate(out_root, samples, wt_seq, mut_seq, amplicon_seq, reference_day=None, reference_treatment=None):
    dfs = {}
    for rep in samples:
        rep_path = Path(out_root) / rep
        if not rep_path.is_dir():
            dfs[rep] = pd.DataFrame()
            continue
        subdirs = [d for d in rep_path.iterdir() if d.is_dir()]
        if not subdirs:
            data = process_treatment_folder(str(rep_path), wt_seq, mut_seq, amplicon_seq)
            if not data:
                dfs[rep] = pd.DataFrame()
                continue
            df = pd.DataFrame([data])
            df["Sensitivity"] = calculate_sensitivity(df, "MUT/WT*%", reference_day, reference_treatment)
            df["Replicate"] = rep
            cols = ["Replicate"] + [c for c in df.columns if c != "Replicate"]
            df = df.loc[:, cols]
            if "ReportPath" in df.columns:
                df = df.drop(columns=["ReportPath"])
            dfs[rep] = df
            continue
        rows = []
        for t in subdirs:
            d = process_treatment_folder(str(t), wt_seq, mut_seq, amplicon_seq)
            if d:
                rows.append(d)
        if not rows:
            dfs[rep] = pd.DataFrame()
            continue
        df = pd.DataFrame(rows)
        df["Sensitivity"] = calculate_sensitivity(df, "MUT/WT*%", reference_day, reference_treatment)
        df.insert(0, "Replicate", rep)
        if "ReportPath" in df.columns:
            df = df.drop(columns=["ReportPath"])
        dfs[rep] = df
    return dfs