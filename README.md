# DeepMALDI spectral processing
This code is the official implementation of [Semi-Quantitative MALDI Meaurements of Blood-Based Samples for Molecular Diagnostics](https://doi.org/10.3390/molecules27030997)

## Installation
Download the entire repository to a single working folder. Relative file paths are used (not absolute), so the MATLAB code can run directly from the MATLAB folder. Modifying or moving contained folders will cause the code to fail.

## Important Notes
This processing code is written up for Windows. Directory structures will be different for Mac/Linux operating systems (need to change "\" in the file directories to "/" for this code and associated codes for reading in and writing files.

## Running the spectral processing code
- Open the DeepMALDI_SpectralProcessing.mlx live script located in the downloaded MATLAB folder.
- Ensure that lines 4, 5, and 6 point to the appropriate folders for testing
- For new spectra to analyze, change the file path of data_path to the appropriate folder
  - data_path is the location of the directory containing the batch folder(s) of DeepMALDI spectra. To run, this directory must only contain folders/files along the lines of 
    - Batch1\SampleA\SampleA.txt
    - Batch1\SampleB\SampleB.txt
    - ...
    - BatchN\SampleX\SampleX.txt
  - Where the Batch# directories only have sample folders, which only have sample files. See the ExampleSpectra folder as an example
- Change the results output folder (Line 5) to the appropriate output file location. 

Code is set up only to run with the Biodesix file type read in. This code works if the input spectra is a .txt file with the pectral data flagged by a "#" and input as m/z in the first column and intensity in the second column. It can contain any header as long as the "#" symbol is not used until the line immediately before the data. For different file types, alternative file read in methods will be needed.

Once the file paths are correctly input, the user can run the MATLAB live script from the start.

## Output of the code
The code will write several files to the res_path folder location.
- Aligned_samples.csv: A list of all files that were able to be aligned.
- FT_Amplitude+Bumps.xlsx and FT_Amplitude+Bumps_TIC.xlsx are feature tables showing the raw and TIC normalized feature values, respectively, of each feature using the Amplitude and Bumps method.
- FT_FineAmplitude.xlsx and FT_FineAmplitude_TIC.xlsx are the feature tables showing the raw and TIC normalized feature values, respectively, of each feature using only the Fine strucuture amplitude.
- MasterPeakList.csv is the resulting master peak list for the analyzed samples
- Batch# folder containing
  - Subfolders for each sample spectrum. Within each is
    - ALIGN_COEFFS.csv - Alignment coefficients for the 4x regions used for alignment
    - CHUNK_IDX.csv - An index of the starts of each chunk used in analysis
    - PEAKS.csv - A list of peaks found in the individual spectrum listing the relative index within a given chunk (index), the m/z location (m0), the peak amplitude (A), integrated peak area (int), signal to noise ratio (SN), and multiplicity (mult)
    - SPECTRUM.csv - The sub-components of the spectra as analyzed by the code with m/z (mzs), overall intensity (its), background (BG), background corrected intensity (BGSUB), Fine structure (FINE), and Bumps structure (BUMPS).
