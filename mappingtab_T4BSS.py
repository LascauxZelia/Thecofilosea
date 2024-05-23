import os

dir = "/home/zelia/Desktop/Nanopore/T4BSS/mafft/"

# Path of the output table file
output_file = "mapping_tab.txt"

# Function to extract lines starting with '>' and store them in an array
def extract_lines(file):
    lines = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                lines.append(line.strip())
    return lines

# Open the output file in write mode
with open(output_file, 'w') as f_out:
    # Write the table header
    f_out.write("Gene\tID\tOrganism\n")
    
    # Traverse through all files in the directory
    for filename in os.listdir(dir):
        # Check if the file is a .mafft file
        if filename.endswith(".mafft"):
            file_name = os.path.splitext(filename)[0]
            file_name = file_name.replace("_all", "")
            file_path = os.path.join(dir, filename)
            # Extract lines starting with '>'
            lines = extract_lines(file_path)
            # Write data into the output file
            for line in lines:
                # Split the line into columns using space as separator
                columns = line.split()
                # Write data into the output file
                f_out.write(f"{file_name}\t{columns[0][1:]}\t{' '.join(columns[1:])}\n")
