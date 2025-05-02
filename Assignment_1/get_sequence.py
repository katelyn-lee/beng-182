#Takes a short sequence as query, and, using the index files created in the previous problem, 
#returns the gi number of the database sequence containing the query string.
# Call this program in this format: python3 get_sequence.py [query string]

import sys

# My algorithm:
#     I grab all of data_seq.txt, and scan (left to right) for a matching substring that's same length as 
#     the query with a for loop. As I pass by an "@", I keep track of the sequence start position by updating 
#     my offset variable to it's position+1. Once I find a match, I look for corresponding gi number of
#     my matched_offset variable.

try:
    query = sys.argv[1]
    k = len(query)
    offset = 0
    match_offset = None
    with open('data_seq.txt','rt') as data_seq, open('data_in.txt','rt') as data_in:
        content = data_seq.readline()
        n = len(content)
        
        for i in range(n-k+1):
            if content[i] == '@':
                offset = i+1
            if '@' in content[i:i+k]:
                continue
            elif content[i:i+k] == query:
                match_offset = offset
                # print(match_offset)

        if match_offset != None:
            for line in data_in:
                if str(match_offset) in line:
                    end_pos = line.index(' ')
                    print(f"Match found: query in sequence with gi number {line[:end_pos]}")
        else:
            print(f"No match found: {query} is not a substring of any sequence in datafile.txt")

except Exception as e:
    print(f"An error occurred: {e}")

