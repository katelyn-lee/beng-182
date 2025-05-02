#Database index creation. Given a FASTA file, create two index files.
  #(a)first file data_seq.txt is a concatenation of all sequences in the file, with @ in between each sequence
  #(b)The second file data_in.txt contains lines, each with 2 terms representing a sequence's ‘gi number’ and the offset in “data_seq.txt” where the sequence starts
#call the program with format: python3 database_index_creation.py datafile.txt

import sys

# My algorithm:
#     data_in is filled out with header line, then each subsequent sequence line is added to data_seq. 
#     I use offset to keep track of characters in data_seq for ease in adding offset info to data_in.txt

def find_gi(header) -> str:
    start_pos = header.find('|') + 1
    end_pos = header.find('|',start_pos)
    return header[start_pos:end_pos]

offset = 0 

try:
    infile_name = sys.argv[1]

    with open(infile_name,'rt') as file, open('data_seq.txt','w') as data_seq, open('data_in.txt','w') as data_in:
        sequence = ''
        for line in file:
            if line[0] == '>': 
                if offset != 0:
                    data_seq.write(sequence + '@')
                    offset += 1
                    sequence = ''

                # add seq gi num and seq start pos(offset) to data in:
                data_in.write(find_gi(line) + ' ')
                data_in.write(str(offset) + '\n')

            else:
                if line.strip() != '':
                    sequence = sequence + line.strip()
                    offset += len(line.strip()) 
        data_seq.write(sequence)

except FileNotFoundError:
    print(f"Error: Input file '{infile_name}' not found.")
except Exception as e:
    print(f"An error occurred: {e}")

