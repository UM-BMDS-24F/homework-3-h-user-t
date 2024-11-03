import os

# Add the BLAST+ bin directory to PATH
blast_bin_dir = r"C:\Program Files\NCBI\blast-2.16.0+\bin"
os.environ['PATH'] += os.pathsep + blast_bin_dir

# Now print the PATH to verify
print("Updated PATH environment variable in Python:")
print(os.environ['PATH'])
