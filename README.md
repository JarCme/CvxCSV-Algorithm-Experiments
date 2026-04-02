# CvxCSV: Code for "Dynamic Independent Component Extraction with Blending Mixing Vector: CvxCSV Algorithm"

This repository provides the MATLAB code for the audio experiments and simulations presented in the paper.

**Important Note:**
The experiments, particularly those described in sections IV.A, IV.B, and IV.C, are computationally intensive. To facilitate reproducibility, pre-computed results from the paper are included in this repository. You can therefore focus on running the analysis scripts to generate the figures and tables directly.

## Repository Structure

The repository is organized into two main directories:

-   `simulations/`: Code for the numerical simulations (Experiments IV.A, IV.B, IV.C).
-   `audio_experiment/`: Code for the audio source separation experiments (Experiment IV.D).


## Dependencies

-   **MATLAB**
-   **[MTIMESX](https://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support)** (precompiled MEX included in this repository)
-   **[OverIVA](https://github.com/onolab-tmu/overiva)** (required only for the audio experiment)
    - The method is implemented in Python, and the repository includes a MATLAB wrapper function to run it through `pyrunfile`.
    - R. Scheibler and N. Ono, "Independent Vector Analysis with More Microphones Than Sources," 2019 IEEE Workshop on Applications of Signal Processing to Audio and Acoustics (WASPAA), New Paltz, NY, USA, 2019, pp. 185-189, doi: 10.1109/WASPAA.2019.8937080.

## How to Reproduce the Results

### Simulations (Experiments IV.A, IV.B, IV.C)

The simulation workflow is a two-step process:
1.  Run the main simulation script (`mtlb_fun.m`) with different configurations.
    - the repository includes all precomputed results presented in the paper.
2.  Run post-processing scripts to generate figures and tables from the saved results.

#### Step 1: Running the Simulations (`cvxCSV_task`) - can be skipped in case use of precomputed results

The `simulations/cvxCSV_task/main_fun.m` function is the main entry point for running the simulations. It was designed for execution on a grid computing cluster (running PBS scheduler), which is reflected in its command-line interface.

syntax in MATLAB:
```matlab
main_fun(current_subtask_index, ...
         total_tasks, ...
         cpu_cores, ...
         config_file);
```
**Parameters:**
-   `current_subtask_index` (e.g., `1`): The index of the current task chunk. For local runs, this is typically `1`.
-   `total_tasks` (e.g., `1`): The total number of chunks the simulation is divided into. For local runs, this is typically `1`.
-   `cpu_cores` (e.g., `8`): The number of CPU cores to utilize for parallel computation.
-   `config_file` (e.g., `'cfg_alpha_eye.m'`): A string specifying the configuration file that sets the parameters for a specific experiment.

To run it on a local multi-core machine, use the following:
```matlab
% Example: running Experiment IV.A divided into 500 tasks and processed in parallel
% parpool(number_of_available_cores) % 

total_tasks = 500;
parfor i = 1:total_tasks
    mtlb_fun(i, total_tasks, 1, 'cfg_alpha_eye.m');
end
```


**Experiment Configurations:**

-   **Experiment IV.A:** `cfg_alpha_eye.m`, `cfg_gamma_eye.m`, `cfg_tau_eye.m`
-   **Experiment IV.B:** `cfg_ini_model_lap.m`
-   **Experiment IV.C:** `cfg_ini_N_test.m`

#### Step 2: Generating Figures and Tables

After running the simulations (or using the provided pre-computed results), you can generate the final figures and tables.

-   **Figure for Experiment IV.A:**
    -   Navigate to `simulations/algs_vs_CRIB/`.
    -   Run `main_generate_figure.m`.

-   **Histograms for Experiment IV.B:**
    -   Navigate to `simulations/SIR_SDR_wrt_ini_eps/`.
    -   Run `main_generate_histograms.m`.

-   **Table for Experiment IV.C:**
    -   Navigate to `simulations/SIR_SDR_wrt_N/`.
    -   Run `main_generate_tex_table.m`.


### Audio Experiment (Experiment IV.D)

1.  Open MATLAB and navigate to the `audio_experiment/` directory.
2.  Run `main.m` to execute the experiment. (This step is optional, as pre-computed results are included in this repository.)
    - All settings are included within the script.
3.  Run `create_table.m` to generate the LaTeX table with the results.