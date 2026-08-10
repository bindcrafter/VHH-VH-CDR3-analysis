# VH-VHH-ANARCI Explorer

VH-VHH-ANARCI Explorer is a Streamlit application for analyzing antibody
VH and VHH sequences from FASTA files using ANARCI.

## Features

- Upload antibody sequences in FASTA format
- Analyze either VH or VHH sequences
- Supports IMGT, Kabat, Chothia, and Martin numbering schemes
- IMGT is used as the default numbering scheme
- Automatically trims sequences to start from the first `QVQ`
- Extracts:
  - FR1
  - CDR1
  - FR2
  - CDR2
  - FR3
  - CDR3
  - FR4
- Calculates CDR3 length distribution
- Performs hierarchical clustering of CDR3 sequences using normalized Levenshtein distance
- Exports results as CSV and PNG files
- Packages all results into a single downloadable ZIP file
- Output files are automatically prefixed with the user-defined job name

# Installation

## Windows users: WSL installation

This application is easiest to run on Windows through Windows Subsystem for Linux (WSL) with Ubuntu.

### 1. Open Ubuntu / WSL

Open the Ubuntu terminal from Windows.

### 2. Install HMMER and Python virtual-environment support

Run:

```bash
sudo apt update
sudo apt install -y hmmer python3-venv
```

Confirm that HMMER is installed:

```bash
hmmscan -h
```

### 3. Download the project

Download or clone this repository.

For example, if the project is located in your Windows Downloads folder:

```bash
cd /mnt/c/Users/YOUR_WINDOWS_USERNAME/Downloads/VH_VHH_CDR3_analysis
```

Replace `YOUR_WINDOWS_USERNAME` with your Windows username.

If the downloaded repository contains another project directory, enter it:

```bash
cd VH_VHH_CDR3_analysis-main
```

### 4. Create a Python virtual environment

Because recent Ubuntu versions protect the system Python environment,
install the Python packages inside a virtual environment.

It is recommended to create the virtual environment in the WSL Linux
home directory rather than under `/mnt/c/`.

```bash
python3 -m venv ~/vhvhh_env
```

Activate it:

```bash
source ~/vhvhh_env/bin/activate
```

After activation, your terminal should show something similar to:

```text
(vhvhh_env) user@computer:~$
```

### 5. Install Python dependencies

From the project directory:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### 6. Run the application

```bash
python -m streamlit run app.py
```

Streamlit should display:

```text
Local URL: http://localhost:8501
```

Open a Windows web browser and go to:

```text
http://localhost:8501
```

The VH-VHH-ANARCI Explorer should now be running.

# Running the App Again

You do not need to reinstall everything each time.

Open Ubuntu / WSL and go to the project directory:

```bash
cd /mnt/c/Users/YOUR_WINDOWS_USERNAME/Downloads/VH_VHH_CDR3_analysis/VH_VHH_CDR3_analysis-main
```

Activate the existing environment:

```bash
source ~/vhvhh_env/bin/activate
```

Then start the application:

```bash
python -m streamlit run app.py
```

# Using the Application

1. Open the Streamlit interface in your browser.
2. Enter a job name.
3. Upload a FASTA file containing VH or VHH sequences.
4. Select VH or VHH.
5. Select a numbering scheme.
6. Run the analysis.
7. Review the CDR3 analysis and clustering results.
8. Download the ZIP file containing the output files.

# Stopping the Application

Return to the WSL terminal and press:

```text
Ctrl + C
```

# Notes

- HMMER must be installed because it is required by ANARCI.
- Keep the WSL terminal open while the Streamlit application is running.
- A message such as `gio: http://localhost:8501: Operation not supported` does not necessarily indicate an application error under WSL. Open `http://localhost:8501` manually in your Windows browser.
- Avoid installing the Python dependencies globally with `--break-system-packages`. Using a virtual environment is safer.

# Repository Files

The repository should contain at least:

- `app.py`
- `requirements.txt`
- `README.md`
