import yaml
import os

DATA_DIR = os.environ.get('DATA_DIR', 'data')

def load_config(config_file="config.yaml"):
    """Loads configuration from a YAML file."""
    current_dir = os.path.dirname(os.path.realpath(__file__))
    config_file = os.path.join(current_dir, config_file)
    with open(config_file, 'r') as f:
        try:
            config = yaml.safe_load(f)
            return config
        except yaml.YAMLError as e:
            print(f"Error loading config file: {e}")
            return None

config = load_config() # Load configuration at module import time

if config:
    VERSIONS = config.get('versions', {}) 
    GO_CODES = config.get('go_codes', {}) 
    RAW_FILE_PATHS = config.get('raw_file_paths', {})
    for key in RAW_FILE_PATHS:
        RAW_FILE_PATHS[key] = os.path.join(DATA_DIR, RAW_FILE_PATHS[key])
    PROCESSED_PATHS = config.get('processed_paths', {})
    for key in PROCESSED_PATHS:
        PROCESSED_PATHS[key] = os.path.join(DATA_DIR, PROCESSED_PATHS[key])
else:
    VERSIONS = {}
    GO_CODES = {}
    RAW_FILE_PATHS = {}
    PROCESSED_PATHS = {}
