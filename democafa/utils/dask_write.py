#!/usr/bin/env python3

import pandas as pd
import dask.dataframe as dd
import os
import gzip  # Import the gzip module explicitly
import dask

def write_dask_dataframe_to_gzipped_tsv(ddf: dd.DataFrame, output_path: str, blocksize: str = "64MB"):
    """
    Writes a Dask DataFrame to a gzipped TSV file.  
    
    Args:
        ddf: The Dask DataFrame to write.
        output_path: The path to the output gzipped TSV file (e.g., 'output.tsv.gz').
        blocksize: The Dask blocksize for parallel processing (default: "64MB"). Adjust based on your memory.
    """
    
    # 1. Define helper functions for partition writing
    def write_partition_no_header(partition: pd.DataFrame, filename: str, compression='gzip'):
        """Writes a single DataFrame partition to a gzipped TSV file (no header)."""
        try:
            partition.to_csv(filename, sep='\t', index=False, header=False, compression=compression)
        except Exception as e:
            print(f"Error writing partition to {filename}: {e}")
            raise  # Re-raise the exception to signal failure
            
    # 2. Generate temporary filenames (on the same filesystem!)
    temp_prefix = os.path.join(os.path.dirname(output_path), "temp_partition")
    temp_files = [f"{temp_prefix}_{i}.tsv.gz" for i in range(ddf.npartitions)]
    
    # 3. Write the header to a separate file
    header_file = os.path.join(os.path.dirname(output_path), "header.tsv")
    first_partition = ddf.partitions[0].compute()  # Extract the *actual* DataFrame
    first_partition.to_csv(header_file, sep='\t', index=False, header=True)
    
    # 4. Write the remaining partitions in parallel (without headers)
    delayed_writes = []
    for i in range(1, ddf.npartitions):
        partition = ddf.partitions[i]
        filename = temp_files[i]
        delayed_write = dask.delayed(write_partition_no_header)(partition, filename)
        delayed_writes.append(delayed_write)
    
    dask.compute(delayed_writes)  # Execute the delayed writes in parallel
    
    # 5. Concatenate the gzipped partitions into a single gzipped file
    def concatenate_gzipped_files(input_files, output_file, header_file):
        """Concatenates gzipped files, prepending the header and removing temp files."""
        with open(output_file, 'wb') as outfile:
            # Add the header
            with open(header_file, 'rb') as header_infile:  # Read as binary!
                outfile.write(header_infile.read())  # No need to encode, already bytes
            
            # Concatenate gzipped data (skipping gzip header)
            for filename in input_files:
                with open(filename, 'rb') as infile:  # Read as binary!
                    # Skip the gzip header (first 10 bytes) and the extra flag (next 2 bytes)
                    infile.read(10)
                    outfile.write(infile.read())
                os.remove(filename)  # Clean up
    
    # Concatenate the files
    concatenate_gzipped_files(temp_files[1:], output_path, header_file)
    
    # Remove the header file
    os.remove(header_file)
    
    print(f"Successfully wrote gzipped TSV file to: {output_path}")