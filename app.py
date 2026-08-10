import io
import os
import re
import shutil
import tempfile
import zipfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
from Bio import SeqIO
from anarci import anarci
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


st.set_page_config(page_title="VHH-ANARCI Explorer", layout="wide")
st.title("VHH-ANARCI Explorer")
st.caption("Upload a FASTA file to annotate VH/VHH regions and analyze CDR3 sequences.")

REGION_RANGES = {
    "imgt": {
        "FR1": (1, 26), "CDR1": (27, 38), "FR2": (39, 55),
        "CDR2": (56, 65), "FR3": (66, 104), "CDR3": (105, 117), "FR4": (118, 128),
    },
    "kabat": {
        "FR1": (1, 30), "CDR1": (31, 35), "FR2": (36, 49),
        "CDR2": (50, 65), "FR3": (66, 94), "CDR3": (95, 102), "FR4": (103, 113),
    },
    "chothia": {
        "FR1": (1, 25), "CDR1": (26, 32), "FR2": (33, 51),
        "CDR2": (52, 56), "FR3": (57, 94), "CDR3": (95, 102), "FR4": (103, 113),
    },
    "martin": {
        "FR1": (1, 25), "CDR1": (26, 35), "FR2": (36, 49),
        "CDR2": (50, 58), "FR3": (59, 94), "CDR3": (95, 102), "FR4": (103, 113),
    },
}

OUTPUT_COLUMNS = ["title", "FR1", "FR2", "FR3", "FR4", "CDR1", "CDR2", "CDR3"]


def safe_job_id(name):
    name = name.strip()
    if not name:
        raise ValueError("Please enter a job name.")
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")
    if not value:
        raise ValueError("Job name must contain at least one letter or number.")
    return value


def clean_sequence(seq):
    seq = "".join(str(seq).split()).upper().replace("*", "")
    start = seq.find("QVQ")
    return seq[start:] if start >= 0 else ""


def clean_title(description):
    # Keep candidate_N when present; otherwise keep text before the first |.
    m = re.search(r"(candidate[_-]?\d+)", description, flags=re.I)
    if m:
        return m.group(1).replace("-", "_")
    return description.split("|")[0].strip()


def pick_heavy_domain(numbered_entry, detail_entry):
    if not numbered_entry:
        return None, None
    for i, domain in enumerate(numbered_entry):
        details = detail_entry[i] if detail_entry and i < len(detail_entry) else {}
        if details.get("chain_type") == "H":
            return domain[0], details
    return None, None


def split_regions(domain_numbering, scheme_name):
    regions = {k: [] for k in ["FR1", "FR2", "FR3", "FR4", "CDR1", "CDR2", "CDR3"]}
    ranges = REGION_RANGES[scheme_name]
    for position, aa in domain_numbering:
        if aa == "-":
            continue
        base_pos = int(position[0])
        for region, (start, end) in ranges.items():
            if start <= base_pos <= end:
                regions[region].append(aa)
                break
    return {k: "".join(v) for k, v in regions.items()}


def levenshtein(a, b):
    if len(a) < len(b):
        a, b = b, a
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        current = [i]
        for j, cb in enumerate(b, 1):
            current.append(min(
                current[-1] + 1,
                previous[j] + 1,
                previous[j - 1] + (ca != cb)
            ))
        previous = current
    return previous[-1]


def normalized_levenshtein(a, b):
    denom = max(len(a), len(b))
    return 0.0 if denom == 0 else levenshtein(a, b) / denom


def analyze(uploaded_file, job_name, receptor_type, scheme, cluster_threshold):
    job_id = safe_job_id(job_name)
    workdir = tempfile.mkdtemp(prefix=f"{job_id}_")
    output_dir = os.path.join(workdir, f"{job_id}_results")
    os.makedirs(output_dir, exist_ok=True)

    fasta_path = os.path.join(workdir, uploaded_file.name)
    with open(fasta_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError("No FASTA records were found.")

    seq_tuples, failures = [], []
    for record in records:
        title = clean_title(record.description)
        seq = clean_sequence(record.seq)
        if not seq:
            failures.append({"title": title, "reason": "QVQ motif not found"})
            continue
        seq_tuples.append((title, seq))

    if not seq_tuples:
        raise ValueError("No sequences containing the QVQ motif were found.")

    numbered, alignment_details, hit_tables = anarci(
        seq_tuples, scheme=scheme, output=False, assign_germline=False
    )

    rows = []
    for i, (title, seq) in enumerate(seq_tuples):
        domain_numbering, details = pick_heavy_domain(numbered[i], alignment_details[i])
        if domain_numbering is None:
            failures.append({"title": title, "reason": "No ANARCI heavy-chain (H) domain detected"})
            continue
        rows.append({"title": title, **split_regions(domain_numbering, scheme)})

    df = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if df.empty:
        raise ValueError("ANARCI did not detect any heavy-chain sequences.")

    regions_csv = os.path.join(output_dir, f"{job_id}_anarci_regions.csv")
    df.to_csv(regions_csv, index=False)

    if failures:
        pd.DataFrame(failures).to_csv(
            os.path.join(output_dir, f"{job_id}_anarci_failures.csv"), index=False
        )

    # CDR3 length analysis
    cdr3_lengths = df["CDR3"].str.len()
    length_counts = (
        cdr3_lengths.value_counts().sort_index()
        .rename_axis("CDR3_length").reset_index(name="count")
    )
    length_counts.to_csv(
        os.path.join(output_dir, f"{job_id}_cdr3_length_counts.csv"), index=False
    )

    fig1, ax1 = plt.subplots(figsize=(9, 5))
    ax1.bar(length_counts["CDR3_length"], length_counts["count"])
    ax1.set_xlabel("CDR3 length (aa)")
    ax1.set_ylabel("Number of sequences")
    ax1.set_title(f"{job_name} — CDR3 length distribution — {receptor_type}, {scheme.upper()}")
    ax1.set_xticks(length_counts["CDR3_length"].tolist())
    fig1.tight_layout()
    length_png = os.path.join(output_dir, f"{job_id}_cdr3_length_distribution.png")
    fig1.savefig(length_png, dpi=300, bbox_inches="tight")

    # CDR3 clustering
    cluster_df = df.loc[
        df["CDR3"].notna() & (df["CDR3"].astype(str).str.len() > 0),
        ["title", "CDR3"]
    ].copy()
    cluster_df["CDR3_length"] = cluster_df["CDR3"].str.len()

    fig2 = None
    if len(cluster_df) >= 2:
        seqs = cluster_df["CDR3"].tolist()
        n = len(seqs)
        distance_matrix = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i + 1, n):
                d = normalized_levenshtein(seqs[i], seqs[j])
                distance_matrix[i, j] = d
                distance_matrix[j, i] = d

        Z = linkage(squareform(distance_matrix, checks=False), method="average")
        raw_ids = fcluster(Z, t=cluster_threshold, criterion="distance")
        unique_ids = sorted(set(raw_ids))
        remap = {old: new for new, old in enumerate(unique_ids, start=1)}
        cluster_df["cluster_id"] = [remap[x] for x in raw_ids]
        cluster_df["cluster"] = cluster_df["cluster_id"].map(lambda x: f"Cluster {x}")

        cluster_df[["title", "CDR3", "CDR3_length", "cluster"]].to_csv(
            os.path.join(output_dir, f"{job_id}_cdr3_clusters.csv"), index=False
        )

        fig_height = max(7, 0.35 * n)
        fig2, ax2 = plt.subplots(figsize=(12, fig_height))
        dendrogram(
            Z,
            labels=cluster_df["title"].tolist(),
            orientation="left",
            color_threshold=0,
            above_threshold_color="gray",
            link_color_func=lambda k: "gray",
            ax=ax2,
        )
        ax2.axvline(cluster_threshold, linestyle="--", linewidth=1.5, color="black",
                    label=f"Cutoff = {cluster_threshold}")

        cmap = plt.get_cmap("tab20")
        cluster_ids = sorted(cluster_df["cluster_id"].unique())
        cluster_colors = {cid: cmap((cid - 1) % 20) for cid in cluster_ids}
        title_to_cluster = dict(zip(cluster_df["title"], cluster_df["cluster_id"]))
        for tick in ax2.get_yticklabels():
            cid = title_to_cluster.get(tick.get_text())
            if cid is not None:
                tick.set_color(cluster_colors[cid])
                tick.set_fontweight("bold")

        ax2.set_title(f"{job_name} — CDR3 hierarchical clustering")
        ax2.set_xlabel("Normalized Levenshtein distance")
        ax2.set_ylabel("Sequence")
        ax2.legend(loc="best")
        fig2.tight_layout()
        fig2.savefig(
            os.path.join(output_dir, f"{job_id}_cdr3_clustering.png"),
            dpi=300, bbox_inches="tight"
        )

    zip_path = shutil.make_archive(
        os.path.join(workdir, f"{job_id}_anarci_analysis_results"),
        "zip",
        output_dir
    )
    with open(zip_path, "rb") as f:
        zip_bytes = f.read()

    return df, length_counts, cluster_df, fig1, fig2, zip_bytes, job_id


with st.sidebar:
    st.header("Analysis settings")
    job_name = st.text_input("Job name", value="my_project")
    receptor_type = st.selectbox("Antibody type", ["VHH", "VH"])
    scheme = st.selectbox("Numbering scheme", ["imgt", "kabat", "chothia", "martin"], index=0)
    cluster_threshold = st.number_input(
        "CDR3 cluster cutoff",
        min_value=0.01, max_value=1.0, value=0.35, step=0.01
    )

st.subheader("1. Upload FASTA")
uploaded_file = st.file_uploader(
    "Choose a FASTA file",
    type=["fasta", "fa", "faa", "fas", "txt"],
    help="Sequences are trimmed to start at the first QVQ motif."
)

if uploaded_file is None:
    st.info("Upload a FASTA file to begin.")
else:
    st.success(f"Selected: {uploaded_file.name}")

    if st.button("Run analysis", type="primary"):
        try:
            with st.spinner("Running ANARCI and CDR3 analysis..."):
                df, length_counts, cluster_df, fig1, fig2, zip_bytes, job_id = analyze(
                    uploaded_file, job_name, receptor_type, scheme, cluster_threshold
                )

            st.success(f"Analysis complete: {len(df)} heavy-chain sequences")

            st.subheader("2. FR/CDR regions")
            st.dataframe(df, use_container_width=True)

            st.subheader("3. CDR3 length distribution")
            st.pyplot(fig1)
            st.dataframe(length_counts, use_container_width=True)

            if fig2 is not None:
                st.subheader("4. CDR3 clustering")
                st.pyplot(fig2)
                st.dataframe(
                    cluster_df[["title", "CDR3", "CDR3_length", "cluster"]]
                    .sort_values(["cluster", "title"]),
                    use_container_width=True
                )
            else:
                st.warning("At least 2 valid CDR3 sequences are required for clustering.")

            st.download_button(
                "Download all results (.zip)",
                data=zip_bytes,
                file_name=f"{job_id}_anarci_analysis_results.zip",
                mime="application/zip",
            )

        except Exception as e:
            st.exception(e)
