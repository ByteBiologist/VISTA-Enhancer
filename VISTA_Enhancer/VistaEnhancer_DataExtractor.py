#!/usr/bin/env python

import requests
import re
import subprocess
import os
import glob

# Function to get the strand information using Samtools faidx
def get_strand(sequence):
    fa = "hg19.fa"
    output = subprocess.check_output(["samtools", "faidx", fa, sequence, "--mark-strand", "sign"]).decode("utf-8")
    strand = output.split('\n')[0].replace('>', '')
    # Extract the value inside the parentheses - chr16:78510608-78511944(+)
    value = re.search(r'\((.*?)\)', strand).group(1)
    return value

# Send a GET request to the webpage
url = 'https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page=1;show=1;form=ext_search;order=;search.gene=;search.result=yes;search.status=Positives;action=search;search.org=Human;page_size=100;search.sequence=1'
response = requests.get(url)
text_content = response.text

# Find all entries starting with ">Human"
entries = re.findall(r'>Human.*?(?=^>Human|\Z)', text_content, re.MULTILINE | re.DOTALL)

# Remove '</pre' from the last entry
last_entry = entries[-1]
entries[-1] = re.sub('</pre', '', last_entry)

# Define the output file name
output_file = 'VISTA_Human_enhancers_sequences.bed'

# Create the output file and write the header
with open(output_file, 'w') as file:
    # Write the header
    headers = ['#chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'annotation', 'sequence', 'expressionPattern']
    file.write('\t'.join(headers) + '\n')

    # Create a set to store the unique expression patterns
    tissue_categories = set()

    # Process and write each entry to the output file
    for entry in entries:
        # Remove '>' character and newlines, replace '|' with tab '\t' (removing surrounding spaces)
        entry = re.sub(r'[>\n]|(\s*\|\s*)', lambda m: '\t' if m.group(1) else '', entry)
        # Split the columns
        columns = entry.split('\t')
        # Rearrange the last column (ear[6/6]CTCCCCTgg)
        last_column_parts = re.split(r'(\])', columns[-1])
        # Combine the last_column_parts with the rest of the columns
        combined_columns = columns[:-1] + [last_column_parts[0] + last_column_parts[1], last_column_parts[2]]
        # Join the combined columns with '\t' separator
        combined_entry = '\t'.join(combined_columns)
        # Extract entries between 4th and last column
        column_4_to_last = combined_entry.split('\t')[4:-1]
        # Combine entries between 4th and last column
        joined_column_4_to_last = ';'.join(column_4_to_last)
        # Split the joined_column_4_to_last by ';'
        split_last_column = joined_column_4_to_last.split(';')
        # Get strand
        strand = get_strand(combined_entry.split('\t')[1])
        # Split sequence identifier into chromosome, start and end
        chromosome, position = combined_entry.split('\t')[1].split(':')
        start_position, end_position = position.split('-')

        # Iterate over expression patterns
        for item in split_last_column:
            # Parsing the expression patters ('trigeminal V (ganglion, cranial)[3/5]')
            expression_pattern = item.strip().split('[')[0].strip().replace(' ', '_').replace(',','').capitalize().replace('(','').replace(')','')
            rearranged_columns = [chromosome, start_position, end_position, combined_entry.split('\t')[2], strand, combined_entry.split('\t')[3], combined_entry.split('\t')[-1], item]
            # Write the entry with rearranged columns
            file.write('\t'.join(rearranged_columns) + '\n')
            # Add the expression pattern to the set
            tissue_categories.add(expression_pattern)

# Create the directory to store the expression pattern files
expression_pattern_files_directory = 'Tissue_specific_files'
os.makedirs(expression_pattern_files_directory, exist_ok=True)


# Split the entries based on tissue category
for tissue in tissue_categories:
    tissue_filename = f'{expression_pattern_files_directory}/{tissue}.{output_file}'
    with open(output_file, 'r') as input_file, open(tissue_filename, 'w') as tissue_file:
        # Write the header to the individual tissue files
        tissue_file.write('\t'.join(headers) + '\n')
        for line in input_file:
            # Compare the pattern to the expressionPattern column
            columns = line.split('\t')[-1].strip().split('[')[0].strip().replace(' ', '_').replace(',','').capitalize().replace('(','').replace(')','')
            # pattern_match = columns.strip().split('[')[0].strip().replace(' ', '_')
            if columns == tissue:
                tissue_file.write(line)

# Sort each tissue specific files individually
file_pattern = f"{expression_pattern_files_directory}/*.bed"
file_list = glob.glob(file_pattern)

for file_path in file_list:
#    print (file_path)
    sort_command = f"LC_ALL=C sort -k1,1 -k2,2n -k3,3n {file_path} -o {file_path}"
    subprocess.run(sort_command, shell=True)

# Compress and index the sorted files using bgzip and tabix
for file_path in file_list:
    bgzip_command = f"bgzip {file_path}"
    subprocess.run(bgzip_command, shell=True)

    tabix_command = f"tabix -p bed {file_path}.gz"
    subprocess.run(tabix_command, shell=True)

print("Output has been written to 'VISTA_Human_enhancers_sequences_06142023.bed' file.")
print(f"Tissue specific files have been created in the '{expression_pattern_files_directory}' directory.")

