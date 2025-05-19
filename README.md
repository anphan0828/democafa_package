# democafa

## Overview

This project is a CAFA (Critical Assessment of protein Function Annotation) implementation focusing on protein function prediction. It represents a conversion from a script-based approach to a proper Python package structure, offering improved modularity, reusability, and maintainability.

## Key Features

*   **Efficient Matrix Handling:** Utilizes efficient storage and reuse of large sparse matrices to prevent redundant matrix creation, save computation time, and reduce memory usage.
*   **Modular Design:** Functionality is clearly separated into focused modules, making the project easy to maintain, extend, and reuse.
*   **Configuration Management:** Centralized configuration in Python modules provides type-safe configuration that is easy to import, use, and version control.

## Project Structure
### Directory Breakdown

*   **democafa:** Contains the main source code for the democafa package.
    *   **\_\_init\_\_.py:** Initializes the democafa package, including version information and top-level imports.
    *   **config.py:** Provides configuration settings and parameters for the package.
    *   **config.yaml:** Contains data sources and versions.
    *   **main.py:** Contains the main entry point for the application.
    *   **datacollection:** Contains modules for data retrieval and preparation.
        *   **\_\_init\_\_.py:** Initializes the data module.
        *   **create\_test\_set.py:** Creates the combined test set.
        *   **retrieve\_sequences.py:** Collects sequence data.
        *   **retrieve\_terms.py:** Retrieves GO terms.
    *   **baselines:** Contains modules implementing prediction algorithms.
        *   **\_\_init\_\_.py:** Initializes the baselines module.
        *   **blast.py:** Implements a BLAST-based prediction algorithm.
        *   **goa_nonexp.py:** Implements a non-experimental GO annotations baseline prediction.
        *   **naive.py:** Implements a naive prediction algorithm.
        *   **prott5.py:** Implements a ProtT5-embeddings-based prediction algorith.
    *   **utils:** Contains utility modules.
        *   **\_\_init\_\_.py:** Initializes the utilities module.
        *   **dask\_write.py:** Provides utilities for working with Dask for efficient data writing operations.
        *   **ontology.py:** Provides shared ontology processing functions.
        *   **run_blast.sh:** Runs external NCBI BLAST (shell script) from terminal with bioconda::blast.
        *   **run_prott5.sh:** Runs external ProtT5 embeddings (shell script) from terminal with hugging face model.
        *   **prott5-baseline:** Folder containing scripts for ProtT5 embeddings retrieval.
*   **data:** Contains raw and processed data.        
*   **tests:** Contains test files mirroring the package structure.
*   **build:** Contains build artifacts.

### Package Structure

```
democafa_package/
    data/
        raw/
        processed/
        release/
    democafa/
        __init__.py
        config.py
        config.yaml
        main.py
        datacollection/
            __init__.py
            retrieve_terms.py
            retrieve_sequences.py
            create_test_set.py
        baselines/
            __init__.py
            naive.py
            blast.py
            goa_nonexp.py
            prott5.py
        utils/
            __init__.py
            ontology.py
            dask_write.py
            run_blast.sh
            run_prott5.sh
            prott5-baseline/
    tests/
    build/
    democafaenv.yaml
    pyproject.toml
    setup.py
    README.md
```

## Usage

To use this project, follow these steps:

0.  **Conda environments:** 
*   Create main CAFA environment using `conda env create -f democafaenv.yaml`, then activate it using `conda activate democafaenv`. 
*   Create another environment to run prott5:
```
conda create --name prott5-env python=3.10 
conda activate prott5-env
conda install pytorch torchvision torchaudio -c pytorch
conda install h5py transformers tqdm joblib scikit-learn scipy sentencepiece -c conda-forge
```
1.  **Installation:** Install the package using `pip install .` (from the `democafa_package` directory).
2.  **Configuration:** Configure the necessary settings in the `config` modules.
3.  **Data Preparation:** Run the scripts in the `datacollection` module to retrieve and prepare the data.
4.  **Baseline Prediction:** Use the scripts in the `baselines` module to make baseline predictions.
5.  **Evaluation:** Evaluate the predictions using the appropriate evaluation scripts.

To rerun any step within the pipeline, run the following command from the root directory (`democafa_package/`)

```
python3 -m democafabaselines.goa_nonexp <provide arguments as needed>
```

## Additional notes
- Baselines are constructed using GO at submission deadline (t0)
    - To evaluate baselines: propagate using t0 obo, then intersect with obo (t-1), then propagate obo to ensure connectivity
    - To evaluate any prediction: propagate using t-1 (the pivot graph)
    - To make ground truth comparable with the predictions: propagate ground truth using t1 obo, then intersect with t-1 obo, then propagate using t-1 obo
    - To use GO-slim: propagate predictions using GO slim obo (with roots included), intersect with t-1 obo, then propagate using t-1 obo
    
    
# TODO:
- Allow the use of GAF file for predictors (in addition to annotated sparse matrix)
- Add example data
- Add unit test
- Use logging module
- Proteins will be removed from evaluation set if their sequences changed between release and evaluation
- Separate script for run_blast and run_prott5, but access to config paths
- Add rerun argument to each of the script
