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
    *   **predictors:** Contains modules implementing prediction algorithms.
        *   **\_\_init\_\_.py:** Initializes the predictors module.
        *   **blast.py:** Implements a BLAST-based prediction algorithm.
        *   **goa_nonexp.py:** Implements a non-experimental GO annotations baseline prediction.
        *   **naive.py:** Implements a naive prediction algorithm.
    *   **utils:** Contains utility modules.
        *   **\_\_init\_\_.py:** Initializes the utilities module.
        *   **dask\_write.py:** Provides utilities for working with Dask for efficient data writing operations.
        *   **ontology.py:** Provides shared ontology processing functions.
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
        predictors/
            __init__.py
            naive.py
            blast.py
            goa_nonexp.py
        utils/
            __init__.py
            ontology.py
            dask_write.py
    tests/
    build/
    pyproject.toml
    setup.py
    README.md
```

## Usage

To use this project, follow these steps:

1.  **Installation:** Install the package using `pip install .` (from the `democafa_package` directory).
2.  **Configuration:** Configure the necessary settings in the `config` modules.
3.  **Data Preparation:** Run the scripts in the `data` module to retrieve and prepare the data.
4.  **Prediction:** Use the scripts in the `predictors` module to make predictions.
5.  **Evaluation:** Evaluate the predictions using the appropriate evaluation scripts.

# TODO:
- Instructions & run BLAST
- Allow the use of GAF file for predictors (in addition to annotated sparse matrix)
- Sparse matrix is currently stored as 3 separate aspects for smaller files
- Add feature to rerun any step of choice
- Add evaluator as a subpackage
- Add installation instructions (via pip or via conda environment)
- Add example data
- Use logging module