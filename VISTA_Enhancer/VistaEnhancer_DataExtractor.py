#!/usr/bin/env python

import requests
import re
import subprocess
import os
import glob
from bs4 import BeautifulSoup

# Function to get the strand information using Samtools faidx
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
def get_strand(sequence):
    fa = "hg19.fa"
    output = subprocess.check_output(["samtools", "faidx", fa, sequence, "--mark-strand", "sign"]).decode("utf-8")
    lines = output.split('\n')
    strand = lines[0].replace('>', '')
    value = re.search(r'\((.*?)\)', strand).group(1)
    ref_sequence = ''.join(lines[1:])
    return value, ref_sequence

# Function to get gene information from experiment page
def get_experiment_info(experiment_id):
    # Construct the URL for the experiment ID
    url = f"https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?form=presentation&show=1&experiment_id={experiment_id}&organism_id=1"

    # Send a GET request to the URL
    response = requests.get(url)
    html_content = response.text

    # Parse the HTML content with BeautifulSoup
    soup = BeautifulSoup(html_content, 'html.parser')

    # Extract the position
    position_element = soup.find('b', text='Position:')
    if position_element:
        position = position_element.find_next_sibling(text=True).strip()
        position = position.replace(' (', '').replace(',', '')
    else:
        position = None

    # Extract the flanking genes
    flanking_genes_element = soup.find('b', text='Flanking genes:')
    if flanking_genes_element:
        flanking_genes = flanking_genes_element.find_next_siblings('a')
        flanking_genes = [gene.text.strip() for gene in flanking_genes]
        flanking_genes = '-'.join(flanking_genes)
    else:
        flanking_genes = None

    # Return the extracted information
    return position, flanking_genes

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
    headers = ['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'annotation', 'sequence', 'expressionPattern', 'Flanking_genes']
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
        combined_columns = columns[0:2] + ['_'.join(columns[2].split())] + columns[3:-1] + [last_column_parts[0] + last_column_parts[1], last_column_parts[2]]
        # Join the combined columns with '\t' separator
        combined_entry = '\t'.join(combined_columns)
        # Extract entries between 4th and last column
        column_4_to_last = combined_entry.split('\t')[4:-1]
        # Combine entries between 4th and last column
        joined_column_4_to_last = ';'.join(column_4_to_last)
        # Split the joined_column_4_to_last by ';'
        split_last_column = joined_column_4_to_last.split(';')
        # Get strand and reference sequence
        strand, ref_sequence = get_strand(combined_entry.split('\t')[1])
        # Split sequence identifier into chromosome, start and end
        chromosome, position = combined_entry.split('\t')[1].split(':')
        start_position, end_position = position.split('-')

        # Iterate over expression patterns
        for item in split_last_column:
            # Parsing the expression patterns ('trigeminal V (ganglion, cranial)[3/5]')
            expression_pattern = item.strip().split('[')[0].strip().replace(' ', '_').replace(',', '').capitalize().replace('(', '').replace(')', '')
            # Extract the score from within the square brackets and convert to decimal
            score_parts = re.search(r'\[(.*?)\]', item).group(1).split('/')
            score_decimal = int(score_parts[0]) / int(score_parts[1])
            # Format the score as a decimal with two decimal places
            score_formatted = "{:.2f}".format(score_decimal)
            # Add the expression pattern to the set
            tissue_categories.add(expression_pattern)
            # Get experiment ID
            element_id = combined_columns[2].split('_')[1]
            # Replace the prefix 'element_' with 'hg'
            element_name = combined_columns[2].replace('element_', 'hs')
            # Get Flanking genes infomation from the experiment page
            position, gene = get_experiment_info(element_id)
            # Handle None value for gene
            flanking_gene = gene if gene else ''
            # Write the entry with rearranged columns
            columns = [chromosome, start_position, end_position, element_name, score_formatted, strand, combined_columns[3], combined_columns[-1], item, ref_sequence, position]

            # Check if the enhancer_sequence match with reference_sequence
            if len(columns) >= 10:
                # Get column 8 and column 10 (ignoring capitalization)
                enh_seq = columns[7].strip()
                ref_seq = columns[9].strip()

                position1 = combined_entry.split('\t')[1]
                position2 = columns[10]

                if position1 == position2:
                    flanking_gene = gene
                else:
                   flanking_gene = "None"

                # Compare column 7 and column 9 (ignoring capitalization)
                if enh_seq.lower() == ref_seq.lower() and position1 == position2:
                    strand = "+"
                    rearranged_columns = [chromosome, start_position, end_position, element_name, score_formatted, strand, combined_columns[3], combined_columns[-1], item, flanking_gene]
                else:
                    strand = '-'
                    rearranged_columns = [chromosome, start_position, end_position, element_name, score_formatted, strand, combined_columns[3], combined_columns[-1], item, flanking_gene]

            # Write the entry with rearranged columns
            file.write('\t'.join(rearranged_columns) + '\n')

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
            expression_pattern = line.split('\t')[-2].strip().split('[')[0].strip().replace(' ', '_').replace(',','').capitalize().replace('(','').replace(')','')
            # pattern_match = columns.strip().split('[')[0].strip().replace(' ', '_')
            if expression_pattern == tissue:
                tissue_file.write(line)

# Sort each tissue specific files individually
file_pattern = f"{expression_pattern_files_directory}/*.bed"
file_list = glob.glob(file_pattern)

for file_path in file_list:
    #print (file_path)
    sort_command = f"LC_ALL=C sort -k1,1 -k2,2n -k3,3n {file_path} -o {file_path}"
    subprocess.run(sort_command, shell=True)

# Compress and index the sorted files using bgzip and tabix
for file_path in file_list:
    bgzip_command = f"bgzip -f {file_path}"
    subprocess.run(bgzip_command, shell=True)

    tabix_command = f"tabix -f -p bed {file_path}.gz"
    subprocess.run(tabix_command, shell=True)

print("Output has been written to 'VISTA_Human_enhancers_sequences.bed' file.")
print(f"Tissue specific files have been created in the '{expression_pattern_files_directory}' directory.")

## Processed Output
## Ouput is a tissue specific enhancer sequences with expression pattern in BED like format. (BED6+4)

# BED format
print("#Processed Output")

with open(output_file, 'r') as file:
    first_line = file.readline().strip()
    second_line = file.readline().strip()
# Print the first two lines
print(first_line)
print(second_line)

### Tissue categories:
### 21 uniqe tissue categories

# Tissue categories
print ("Tissue categories")
sorted_categories = sorted(tissue_categories, key=str.lower)  # Sort categories while ignoring case
unique_categories = []
for category in sorted_categories:
    capitalized_category = category.capitalize()
    formatted_category = capitalized_category.replace("_", " ")
    if category.lower() not in unique_categories and formatted_category not in unique_categories:
        unique_categories.append(formatted_category)

for category in unique_categories:
    print(category)


