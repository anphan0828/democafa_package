#!/bin/bash

# Required packages: yq for parsing yaml
if ! command -v yq &> /dev/null; then
    echo "Error: yq is not installed. Please install it using:"
    echo "  wget https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64 -O /usr/local/bin/yq && chmod +x /usr/local/bin/yq"
    exit 1
fi

CONFIG_FILE="democafa/config.yaml"
OUTPUT_DIR=$1

mkdir -p $OUTPUT_DIR
echo "Reading URLs from $CONFIG_FILE..."

# Loop through each key-value pair in the versions section
yq eval '.versions | to_entries | .[] | .key + "|" + .value' $CONFIG_FILE | while IFS="|" read -r key value; do
    # Check if the value is a URL (starts with http or https)
    if [[ $value == http* ]]; then
        # TODO: Check if file already exists
        if [ -f "$OUTPUT_DIR/$key" ]; then
            echo "Skipping $key: File already exists in $OUTPUT_DIR"
            continue
        fi
        echo "Downloading: $value, filename: $OUTPUT_DIR/$key"
        
        wget -q --show-progress --progress=bar:force:noscroll "$value" -O "$OUTPUT_DIR/$key"
        
        if [ $? -eq 0 ]; then
            echo "Downloaded $key"
        else
            echo "Failed to download $key"
        fi
        echo "----------------------------------------"
    else
        echo "Skipping $key: Not a URL"
    fi
done

echo "Download process completed!"