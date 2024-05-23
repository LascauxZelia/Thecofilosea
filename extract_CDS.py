# Reading file names from Fisci_genes_name.txt
with open('Fisci_genes_name.txt', 'r') as names_file:
    file_names = [line.strip() for line in names_file]

# Reading the .faa file and creating files for each protein
with open('Fisci_all_genes.faa', 'r') as input_file:
    proteins = input_file.read().split('>')[1:]  # Splitting the proteins
    
    # Iterating through each protein
    for index, protein in enumerate(proteins):
        header, sequence = protein.split('\n', 1)
        header = header.strip() #+ " Fisci"  # Adding "Fisci" at the end of the header
        sequence = sequence.replace('\n', '')

        # Getting the corresponding file name from Fisci_genes_name.txt
        file_name = file_names[index] + '_Fisci.faa'
        
        # Writing the protein into a separate .faa file
        with open(file_name, 'w') as output_file:
            output_file.write('>' + header + '\n')
            output_file.write(sequence + '\n')
