import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.optimize import nnls
import pandas as pd

# Load .ab1 file and extract sequence trace data
def load_sequence_data(file_path):
    """
    Load sequencing data from an .ab1 file.
    
    Parameters:
        file_path (str): The file path to the .ab1 file containing sequencing data.
    
    Returns:
        tuple: A dictionary containing channel data (A, C, G, T) and the base sequence as a string.
    
    Function:
        This function reads an .ab1 sequencing file and extracts the raw sequencing trace data for the A, C, G, and T channels.
        It also extracts the base sequence, which provides the nucleotide sequence called from the trace data.
    """
    with open(file_path, "rb") as file:
        record = SeqIO.read(file, "abi")
        channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']  # A, C, G, T channels
        channel_data = {channel: np.array(record.annotations['abif_raw'][channel]) for channel in channels}
        base_sequence = record.annotations['abif_raw']['PBAS1'].decode("ascii")  # Extract base sequence
        return channel_data, base_sequence

# Align the sequences to find the break site
# Use a default alignment window to search for the gRNA
def align_sequences(base_sequence, gRNA, alignment_window=(100, 190)):
    """
    Align the base sequence to find the break site using the gRNA sequence.
    
    Parameters:
        base_sequence (str): The nucleotide sequence of the control genome.
        gRNA (str): The guide RNA sequence to locate the break site.
        alignment_window (tuple): A tuple specifying the start and end indices for searching the gRNA.
    
    Returns:
        int: The index position of the break site within the sequence.
    
    Function:
        This function searches for the gRNA sequence within a specified alignment window of the base sequence to find the break site.
        If the gRNA is not found in the specified window, it falls back to searching the entire sequence.
    """
    search_region = base_sequence[alignment_window[0]:alignment_window[1]]
    start_idx = search_region.find(gRNA)
    if start_idx == -1:
        # Fallback to searching the entire sequence if gRNA not found in the specified window
        start_idx = base_sequence.find(gRNA)
        if start_idx == -1:
            raise ValueError("gRNA not found in the control sequence")
    else:
        start_idx += alignment_window[0]
    return start_idx + len(gRNA)  # Break site is at the end of gRNA

# Plot chromatogram to visualize differences between sequences
def plot_chromatogram(control_trace, edited_trace):
    """
    Plot chromatogram to visualize differences between control and edited sequences.
    
    Parameters:
        control_trace (dict): Dictionary containing trace data for the control sequence (A, C, G, T channels).
        edited_trace (dict): Dictionary containing trace data for the edited sequence (A, C, G, T channels).
    
    Function:
        This function plots the chromatogram of both the control and edited sequences to visually compare the signal intensity of each nucleotide channel.
    """
    plt.figure(figsize=(10, 8))
    channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']  # A, C, G, T channels
    colors = {'DATA9': 'green', 'DATA10': 'blue', 'DATA11': 'black', 'DATA12': 'red'}
    labels = {'DATA9': 'A', 'DATA10': 'C', 'DATA11': 'G', 'DATA12': 'T'}
    
    for channel in channels:
        plt.plot(control_trace[channel], label=f"Control {labels[channel]}", color=colors[channel], alpha=0.5)
        plt.plot(edited_trace[channel], label=f"Edited {labels[channel]}", linestyle='--', color=colors[channel], alpha=0.7)
    
    plt.title('Chromatogram - Control vs Edited Sequence')
    plt.xlabel('Base Pair Position')
    plt.ylabel('Signal Intensity')
    plt.legend()
    plt.show()

# Construct the trace matrix A for indel analysis
def construct_trace_matrix(control_trace, max_indel_size=10):
    """
    Construct the trace matrix for indel analysis.
    
    Parameters:
        control_trace (dict): Dictionary containing trace data for the control sequence (A, C, G, T channels).
        max_indel_size (int): Maximum size of indels to be considered in the analysis.
    
    Returns:
        numpy.ndarray: A matrix where each row represents a possible sequence trace, including wild type and possible indel variants.
    
    Function:
        This function constructs a trace matrix where the first row represents the wild type gene trace and subsequent rows represent possible insertion and deletion outcomes.
    """
    indel_matrix = []
    # First row is the wild type gene trace with 100% representation of A, T, C, G
    for channel in control_trace:
        nucleotide_trace = np.zeros_like(control_trace[channel])
        nucleotide_trace += 100  # Representing 100% accuracy for each nucleotide
        indel_matrix.append(nucleotide_trace)
    
    # Subsequent rows represent possible mutation outcomes (insertions/deletions)
    for size in range(-max_indel_size, max_indel_size + 1):
        if size == 0:
            continue  # Skip the original wild type row
        for channel in control_trace:
            shifted_trace = np.roll(control_trace[channel], size)
            if size < 0:
                shifted_trace[size:] = 0  # Truncate values rolled from the left
            elif size > 0:
                shifted_trace[:size] = 0  # Truncate values rolled from the right
            indel_matrix.append(shifted_trace)

    return np.array(indel_matrix).T

# Decompose indel spectrum using non-negative least squares
# Apply a default decomposition window
def decompose_indels(control_trace, edited_trace, max_indel_size=10, decomposition_window=(220, 690)):
    """
    Decompose the indel spectrum using non-negative least squares (NNLS).
    
    Parameters:
        control_trace (dict): Dictionary containing trace data for the control sequence (A, C, G, T channels).
        edited_trace (dict): Dictionary containing trace data for the edited sequence (A, C, G, T channels).
        max_indel_size (int): Maximum size of indels to be considered in the analysis.
        decomposition_window (tuple): The start and end indices specifying the decomposition window.
    
    Returns:
        numpy.ndarray: The indel frequencies as a percentage.
    
    Function:
        This function constructs the trace matrix and then uses NNLS to solve for the indel spectrum, representing the percentage of sequences with each indel type.
    """
    # Construct matrix A
    indel_matrix = construct_trace_matrix(control_trace, max_indel_size)
    edited_trace_combined = np.concatenate([edited_trace[channel] for channel in edited_trace])
    edited_trace_combined = edited_trace_combined[decomposition_window[0]:decomposition_window[1]]
    indel_matrix = indel_matrix[decomposition_window[0]:decomposition_window[1], :]
    
    # Solve A * X = Y using non-negative least squares
    solution, _ = nnls(indel_matrix, edited_trace_combined)
    return solution / solution.sum() * 100  # Normalize to percentage

# Plot indel spectrum
# Ensure the range and frequencies match properly
def plot_indel_spectrum(indel_frequencies):
    """
    Plot the indel spectrum to visualize the distribution of indels.
    
    Parameters:
        indel_frequencies (numpy.ndarray): The calculated frequencies of each indel as a percentage.
    
    Function:
        This function plots a bar chart representing the percentage of each indel size, with annotations for values greater than zero.
    """
    indel_range = np.arange(-len(indel_frequencies) // 2, len(indel_frequencies) // 2)
    bars = plt.bar(indel_range[:len(indel_frequencies)], indel_frequencies, color="blue", alpha=0.7)
    plt.axhline(0, color='black', linewidth=0.8)
    plt.title("Indel Spectrum")
    plt.xlabel("Indel Size (bp)")
    plt.ylabel("% of Spectrum")
    plt.ylim(0, 100)
    
    # Annotate bars with their values if greater than 0
    for bar in bars:
        yval = bar.get_height()
        if yval > 0:
            # Annotating y value at the top of the bar
            plt.text(bar.get_x() + bar.get_width() / 2, yval + 1, f'{yval:.2f}', ha='center', va='bottom')
            # Annotating x value at the bottom of the bar
            if int(bar.get_x()) != 0:
                plt.text(bar.get_x() + bar.get_width() / 2, -1, f'{int(bar.get_x())}', ha='center', va='top')

    plt.show()

# Plot aberrant sequence signal for quality control
# Use a rolling window to better capture local variability
def plot_aberrant_signal(control_trace, edited_trace, window_size=10):
    """
    Plot aberrant sequence signal for quality control.
    
    Parameters:
        control_trace (dict): Dictionary containing trace data for the control sequence (A, C, G, T channels).
        edited_trace (dict): Dictionary containing trace data for the edited sequence (A, C, G, T channels).
        window_size (int): The size of the rolling window used to smooth the data.
    
    Function:
        This function calculates and plots the aberrant sequence signal, which represents the differences between the control and edited sequence traces.
        A rolling window is used to smooth the data to reduce noise and highlight regions of divergence between the two sequences.
    """
    # Smooth the control and edited traces to reduce noise
    aberrant_signals = {}
    for channel in control_trace:
        control_trace_smooth = pd.Series(control_trace[channel]).rolling(window=window_size, center=True).mean().fillna(0)
        edited_trace_smooth = pd.Series(edited_trace[channel]).rolling(window=window_size, center=True).mean().fillna(0)
        
        # Calculate the aberrant signal as the absolute difference between smoothed traces
        diff = np.abs(control_trace_smooth - edited_trace_smooth)
        aberrant_signal = diff * 100 / control_trace_smooth.max()  # Normalize to percentage
        aberrant_signals[channel] = aberrant_signal
        
        plt.plot(aberrant_signal, label=f"Aberrant Signal {channel}", alpha=0.7)
    
    # Plot vertical line at the break site
    plt.axvline(x=len(control_trace['DATA9']) // 2, color='blue', linestyle="--", label="Break Site")
    plt.title("Quality Control - Aberrant Sequence Signal")
    plt.xlabel("Base Pair Position")
    plt.ylabel("% Aberrant Sequence")
    plt.legend()
    plt.ylim(0, 100)  # Adjust the y-axis to better visualize aberrant sequences
    plt.show()

# Calculate nucleotide probabilities for +1 insertion
# Adjusted to focus specifically on the region following the expected cut site
def calculate_inserted_nucleotide_probabilities(control_trace, edited_trace, break_site, window_size=20):
    """
    Calculate nucleotide probabilities for +1 insertion.
    
    Parameters:
        control_trace (dict): Dictionary containing trace data for the control sequence (A, C, G, T channels).
        edited_trace (dict): Dictionary containing trace data for the edited sequence (A, C, G, T channels).
        break_site (int): The index of the break site within the sequence.
        window_size (int): The size of the window used to analyze the insertion site.
    
    Returns:
        dict: A dictionary containing the probabilities for each nucleotide (A, C, G, T) being inserted, expressed as a percentage.
    
    Function:
        This function calculates the probability of each nucleotide being inserted at the break site by comparing the differences in signal intensity between the control and edited sequences.
        It uses a window starting immediately after the break site to focus specifically on the region where insertions are most likely to occur.
    """
    # Use a window starting immediately after the break site to better estimate insertion probabilities
    window_start = break_site
    window_end = min(len(control_trace['DATA9']), break_site + window_size)
    
    # Focus on the insertion site to determine probabilities
    inserted_diff = {channel: np.sum(np.abs(control_trace[channel][window_start:window_end] - edited_trace[channel][window_start:window_end])) for channel in control_trace}
    total_diff = sum(inserted_diff.values())
    if total_diff == 0:
        return {"A": 0, "C": 0, "T": 0, "G": 0}
    
    # Calculate probabilities for each nucleotide and normalize to match the format of the web tool
    probabilities = {channel: (inserted_diff[channel] / total_diff) * 100 for channel in inserted_diff}
    return {
        "A": probabilities['DATA9'],
        "C": probabilities['DATA10'],
        "G": probabilities['DATA11'],
        "T": probabilities['DATA12']
    }

# MAIN
if __name__ == "__main__":
    control_file = "example1.ab1"
    edited_file = "example2.ab1"
    gRNA = "ATCACTCTCGGCATGGACGA"

    control_trace, control_sequence = load_sequence_data(control_file)
    edited_trace, _ = load_sequence_data(edited_file)

    break_site = align_sequences(control_sequence, gRNA)

    # Plot chromatogram to visualize differences
    plot_chromatogram(control_trace, edited_trace)

    # Decompose indel spectrum using the refined matrix
    indel_frequencies = decompose_indels(control_trace, edited_trace)
    plot_indel_spectrum(indel_frequencies)

    # Plot aberrant sequence signal for quality control
    plot_aberrant_signal(control_trace, edited_trace)

    # Calculate inserted nucleotide probabilities
    nucleotide_probs = calculate_inserted_nucleotide_probabilities(control_trace, edited_trace, break_site)
    print("Inserted Nucleotide Probabilities (%):", nucleotide_probs)