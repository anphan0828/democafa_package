# democafa

## Overview

This project is a CAFA (Critical Assessment of protein Function Annotation) implementation focusing on protein function prediction. It represents a conversion from a script-based approach to a proper Python package structure, offering improved modularity, reusability, and maintainability.

## Project Structure
### Directory Breakdown

*   **democafa:** Contains the main source code for the democafa package.
    *   **\_\_init\_\_.py:** Initializes the democafa package, including version information and top-level imports.
    *   **config.py:** Provides configuration settings and parameters for the package.
    *   **configuration.py:** Handles the loading and management of configuration values.
    *   **main.py:** Contains the main entry point for the application.
    *   **data:** Contains modules for data retrieval and preparation.
        *   **\_\_init\_\_.py:** Initializes the data module.
        *   **create\_test\_set.py:** Creates the combined test set.
        *   **retrieve\_sequences.py:** Collects sequence data.
        *   **retrieve\_terms.py:** Retrieves GO terms.
        *   **processed:** Stores generated and intermediate data files.
        *   **raw:** Stores original input data files.
        *   **release:** Stores data release files.
    *   **predictors:** Contains modules implementing prediction algorithms.
        *   **\_\_init\_\_.py:** Initializes the predictors module.
        *   **blast.py:** Implements a BLAST-based prediction algorithm.
        *   **naive.py:** Implements a naive prediction algorithm.
    *   **utils:** Contains utility modules.
        *   **\_\_init\_\_.py:** Initializes the utilities module.
        *   **dask\_write.py:** Provides utilities for working with Dask for efficient data writing operations.
        *   **ontology.py:** Provides shared ontology processing functions.
    *   **tests:** Contains test files mirroring the package structure.
*   **build:** Contains build artifacts.


## Key Features

*   **Efficient Matrix Handling:** Utilizes efficient storage and reuse of large sparse matrices to prevent redundant matrix creation, save computation time, and reduce memory usage.
*   **Modular Design:** Functionality is clearly separated into focused modules, making the project easy to maintain, extend, and reuse.
*   **Configuration Management:** Centralized configuration in Python modules provides type-safe configuration that is easy to import, use, and version control.

## Project Structure

```
democafa_package/
    democafa/
        __init__.py
        config.py
        config.yaml
        main.py
        data/
            __init__.py
            raw/
            processed/
            release/
            retrieve_terms.py
            retrieve_sequences.py
            create_test_set.py
        predictors/
            __init__.py
            naive.py
            blast.py
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

### Files

*   **pyproject.toml:** Contains project metadata and build configuration.
*   **setup.py:** Package installation script.
*   **README.md:** Project documentation (this file).

## Usage

To use this project, follow these steps:

1.  **Installation:** Install the package using `pip install .` (from the `democafa_package` directory).
2.  **Configuration:** Configure the necessary settings in the `config` modules.
3.  **Data Preparation:** Run the scripts in the `data` module to retrieve and prepare the data.
4.  **Prediction:** Use the scripts in the `predictors` module to make predictions.
5.  **Evaluation:** Evaluate the predictions using the appropriate evaluation scripts.

