# Modal Analysis using Prony Method

This MATLAB script performs modal analysis on impact test data using the Prony method. It is designed to identify modal parameters, such as natural frequencies, damping ratios, and mode shapes, from impulse response functions (IRF) obtained during the impact test.

## Features

- **Data Loading**: Load H1 frequency response functions from a `.mat` file.
- **Frequency Response Visualization**: Plot the frequency response function (H1) on a semilogarithmic scale.
- **MDOF (Multi-Degree-of-Freedom) Processing**: Pre-process the data for multi-degree-of-freedom analysis, crucial for handling systems with multiple modes.
- **Prony Method Application**: Extract modal parameters such as eigenfrequencies and damping coefficients.
- **Mode Shape Visualization**: Plot the mode shapes for a specified number of modes.

## Usage

1. **Load Data**: When prompted, select the `.mat` file containing the H1 data to be analyzed.
   
2. **Set Parameters**:
   - `fsamp`: Sampling frequency in Hz.
   - `df`: Frequency resolution.
   - `m`: Number of samples for the Prony method.
   - `order`: Maximum number of poles used in the Prony method to check stability
   - `poles_num`: Number of poles to consider.

3. **Run the Script**:
   - Execute the script in MATLAB.
   - The script will load the data, process it using the Prony method, and plot the results.

4. **View Results**:
   - The script will generate plots for the frequency response, the identified mode shapes, and more.
   ![alt text](https://github.com/shahin1009/prony/blob/main/pics/Pronymethodpub_01.png?raw=true)
   ![alt text](https://github.com/shahin1009/prony/blob/main/pics/Pronymethodpub_02.png?raw=true)
   ![alt text](https://github.com/shahin1009/prony/blob/main/pics/Pronymethodpub_03.png?raw=true)
   ![alt text](https://github.com/shahin1009/prony/blob/main/pics/Pronymethodpub_05.png?raw=true)


## Functions

- **hankel(IRF, N, m)**: Constructs the Hankel matrix for frequency analysis.
- **V_matrix(roots, m)**: Computes the Vandermonde matrix for the Prony method.
- **Prony(IRF, fsamp, order, m, H1_1, freq)**: Core function implementing the Prony method.


## Notes

- Ensure the H1 data is properly formatted and pre-processed before using this script.
- The script is configured for a specific type of impact test data and may require modifications for other types of data.
