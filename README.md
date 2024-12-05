# Genome Engineering Analysis using TIDE

This repository contains a Python implementation of the TIDE (Tracking of Indels by Decomposition) method for analyzing genome editing outcomes using Sanger sequencing data. The project aims to recreate the functionality of the TIDE web tool, enabling users to quantify indels (insertions and deletions) in edited genomes using custom Python code.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Input Files](#input-files)
- [Example Output](#example-output)
- [Code Overview](#code-overview)
- [Requirements](#requirements)
- [License](#license)

## Introduction
Genome engineering technologies like CRISPR/Cas9 allow for targeted genetic modifications. TIDE is a tool designed to estimate the frequency and spectrum of small insertions and deletions (indels) in a genome after editing. This Python implementation aims to recreate the TIDE method using linear algebra techniques such as non-negative least squares (NNLS) to determine the spectrum and efficiency of gene editing.

## Features
- Loads sequencing data from `.ab1` files and extracts nucleotide signals.
- Aligns sequences to determine the break site using a provided gRNA sequence.
- Decomposes sequencing data into an indel spectrum using NNLS.
- Plots chromatograms, indel spectra, and aberrant sequence signals.
- Calculates inserted nucleotide probabilities for +1 insertions.

## Installation

1. **Clone the repository**:
    ```
    git clone https://github.com/username/genome-engineering-tide.git
    ```
2. **Navigate to the directory**:
    ```
    cd genome-engineering-tide
    ```
3. **Install the required Python packages**:
    You can install the required packages using `pip`:
    ```
    pip install -r requirements.txt
    ```
    Alternatively, if you do not have a `requirements.txt` file, manually install the dependencies:
    ```
    pip install numpy matplotlib biopython pandas
    ```

## Usage

1. **Prepare the Files**:
   - Place the `.ab1` files in the root directory or in a separate folder named `data`.
   - The repository includes Python scripts such as `tide_analysis.py`.

2. **Run the Script**:
   ```
   python tide_analysis.py
   ```

3. **Parameters**:
   - You will need to specify the paths to the control and edited `.ab1` files, as well as the gRNA sequence.
   - Default parameter settings are used for alignment, decomposition, and indel size range, but these can be adjusted in the script if needed.

## Input Files

The `.ab1` files are Sanger sequencing chromatogram files. You should have:
- `example1.ab1`: The control sequence (unedited).
- `example2.ab1`: The edited sequence.

These files will be used to extract the nucleotide sequence signals for subsequent analysis.

## Example Output

The Python code will generate the following output:
1. **Chromatogram Plot**: Comparison between control and edited sequences to visualize differences.
2. **Indel Spectrum Plot**: A bar plot representing the percentage of sequences with each indel size.
3. **Aberrant Sequence Signal Plot**: Used for quality control, showing regions of divergence between control and edited traces.
4. **Inserted Nucleotide Probabilities**: A dictionary indicating the percentage of A, C, G, T being inserted at the break site.

Example of inserted nucleotide probabilities:
```
Inserted Nucleotide Probabilities (%): {'A': 9.83, 'C': 13.00, 'G': 39.01, 'T': 38.16}
```

## Code Overview

- **`tide_analysis.py`**: The main script that performs the analysis. It includes functions to:
  - Load `.ab1` files.
  - Align the sequences to find the break site.
  - Construct the trace matrix and use NNLS to decompose the indel spectrum.
  - Plot chromatograms, indel spectrum, and aberrant sequence signals.
  - Calculate nucleotide probabilities for +1 insertions.

## Requirements

The following Python packages are required to run the code:
- `numpy`: Used for numerical operations and handling the sequence data.
- `matplotlib`: Used for plotting chromatograms, indel spectrum, and quality control plots.
- `biopython`: Used to read the `.ab1` sequencing files.
- `pandas`: Used to perform rolling averages on the trace data to smooth the aberrant sequence signal.
- `scipy`: Used for non-negative least squares (NNLS) analysis.

You can install all requirements using the command mentioned in the [Installation](#installation) section.

## License

I don't know

## Acknowledgments

- **TIDE**: This code is inspired by the TIDE web tool developed by Brinkman et al. (2014). The tool uses a similar approach to estimate indel frequencies and track editing efficiency.
- **CRISPR/Cas9 Genome Editing**: Special thanks to researchers who have contributed to the advancement of CRISPR and other gene editing technologies.
