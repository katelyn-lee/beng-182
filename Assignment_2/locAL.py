# Implemented a local alignment progam, using dynamic programing, t align 2 long DNA sequences with linear-gap penalty
# Must call program by using such commands: python3 locAL.py -q <query_file_name> -d <db_file_name> -m <match> -s <mismatch> -i <indel> -a
#   Input match, mismatch, and indel as positive integers. input query and db as fasta text files, each containing just one fasta sequence
#   Output for the best local alignment is: begin_query, end_query, begin_db, end_db, score, length
#      the -a flag is optional. without this flag, the alignment itself will not be outputted. otherwise, the alignment will be outputted on 
#       the command line and in a text file alignment1.txt

import sys
from typing import List, Dict, Tuple

# python3 locAL.py -q query -d db -m match -s mismatch -i indel -a
# python3 locAL.py -q query1.fasta -d db1.fasta -m 1 -s 10 -i 1 -a
# python3 locAL.py -q testq.fasta -d testdb.fasta -m 1 -s 2 -i 2 -a 0 24 888 912

# receive command line args and initialize variables
query_file, db_file = "", ""
query, db = "", ""
match, mismatch, indel = 0, 0, 0
output = False
forQ1 = False
forQ2 = False
forQ4 = False
s_query, e_query = 0, 0   #forQ4
s_db, e_db = 0, 0         #forQ4


if len(sys.argv) == 11:
    print("Proper commands inputted for locAL.py. Running program...")
    query_file = sys.argv[2]
    db_file = sys.argv[4]
    match = sys.argv[6]
    mismatch = sys.argv[8]
    indel = sys.argv[10]
    forQ1 = True
    print(f"These commands are received: query file name {query_file}, db file name {db_file}, match {match}, mismatch {mismatch}, indel {indel}.")
elif len(sys.argv) == 12:
    print("Proper commands inputted. Running program...") 
    query_file = sys.argv[2]
    db_file = sys.argv[4]
    match = sys.argv[6]
    mismatch = sys.argv[8]
    indel = sys.argv[10]
    forQ1 = True
    output = True   #output alignment on terminal and in out file
    print(f"These commands are received: query file name {query_file}, db file name {db_file}, match {match}, mismatch {mismatch}, indel {indel}.")
elif len(sys.argv) == 13:   #for internal use on Q2 locAL 
    query = sys.argv[2]
    db = sys.argv[4]
    match = int(sys.argv[6])
    mismatch = float(sys.argv[8])
    indel = float(sys.argv[10])
    output = False  #for q2, just need length of alignment, not actaully alignment itself
    forQ2 = True    #gives special not file reading variable assignment
elif len(sys.argv) == 16:   #for internal use on Q4 locAL 
    s_query, e_query = (sys.argv[12]), (sys.argv[13])
    s_db, e_db = (sys.argv[14]), (sys.argv[15])
    query = sys.argv[2]
    db = sys.argv[4]
    match = int(sys.argv[6])
    mismatch = float(sys.argv[8])
    indel = float(sys.argv[10])
    output = True  #for q4, yes must get actual alignment, just only put in file alignment4.txt
    forQ4 = True    #gives special not file reading variable assignment
else: 
    print("improper command input. try again using format: python locAL.py -q <query> -d <db> -m <match> -s <mismatch> -i <indel> -a")

if (forQ4 == False) and (forQ2 == False):
    with open(query_file,'rt') as qf, open(db_file,'rt') as dbf:
        next(qf)
        query = (qf.readline()).strip()
        next(dbf)
        db = (dbf.readline()).strip()
else:
    query = query   #just feels/looks nicer, doesnt do anything
    db = db

# run DP algorithm to get local alignment
begin_query, end_query, begin_db, end_db, score, length  = 0, 0, 0, 0, 0, 0
aligned_q, aligned_db = "", ""

def out_matrix(matrix: List[List[int]]):
    for row in matrix:
        print(row)
    
def make_matrix(q,db):
    m, n = len(q), len(db)
    matrix = []
    for i in range(m+1):
        row = []
        for j in range(n+1):
            row.append(0)
        matrix.append(row)
    return matrix

def local_alignment(match:int, mismatch:int, indel:int, q:str, db:str):
    global begin_query, end_query, begin_db, end_db, score, length
    global aligned_q, aligned_db

    matrix = make_matrix(q,db)
    # print(out_matrix(make_matrix(q,db)))
    m, n = len(q), len(db)

    highest_score = -1
    highest_pos = (-1,-1)

    # fill in DP matrix w local alignment
    for i in range(1,m+1):
        for j in range(1,n+1):
            match_val = -1
            if q[i-1] == db[j-1]:
                match_val = matrix[i-1][j-1] + int(match)
            else:
                match_val = int(matrix[i-1][j-1]) - int(mismatch)

            insert = matrix[i][j-1] - int(indel)
            delete = matrix[i-1][j] - int(indel)
            max_score = max(match_val,insert,delete)

            if max_score < 0:
                matrix[i][j] = 0
            else:
                matrix[i][j] = max_score
            if highest_score <= matrix[i][j]:
                highest_score = max_score   #tracking highest score in matrix as you fill it out
                highest_pos = (i,j)
    # out_matrix(matrix)

    # trace back to reconstruct alignment
    i, j = highest_pos[0],highest_pos[1]
    end_query, end_db = i-1, j-1
    score = highest_score

    if output == True:  #if going to output, then store aligned_sequences
        while i > 0 and j > 0 and matrix[i][j] > 0:
            current_score = matrix[i][j]
            if i > 0 and (matrix[i-1][j] == current_score + int(indel)):    #deletion first
                aligned_q = q[i-1] + aligned_q
                aligned_db = "-" + aligned_db
                i -= 1
                length+=1
            elif j > 0 and (matrix[i][j-1] == current_score + int(indel)):  #then insertion
                aligned_q = "-" + aligned_q
                aligned_db = db[j-1] + aligned_db
                j -= 1
                length+=1
            else:
                aligned_q = q[i-1] + aligned_q
                aligned_db = db[j-1] + aligned_db
                i -= 1
                j -= 1
                length+=1
        begin_query, begin_db = i, j
    else:   #not going to output, so dont store aligned_sequences, wastes memory, especially for q4
        while i > 0 and j > 0 and matrix[i][j] > 0:
            current_score = matrix[i][j]
            if i > 0 and (matrix[i-1][j] == current_score + int(indel)):    #deletion f irst
                i -= 1
                length+=1
            elif j > 0 and (matrix[i][j-1] == current_score + int(indel)):  #insertion
                j -= 1
                length+=1
            else:
                i -= 1
                j -= 1
                length+=1
        begin_query, begin_db = i, j

    return -1

local_alignment(match,mismatch,indel,query,db)

if forQ1 == True:
    print(f"{begin_query}, {end_query}, {begin_db}, {end_db}, {score}, {length}")
elif forQ2 == True:
    print(length)

#output alignment due to -a flag
if forQ1==True and output==True:
    with open("alignment1.txt","w") as file:
        between_full = ""
        last_output = -1

        for i in range(length):
            if aligned_q[i] == aligned_db[i]:
                between_full += "|"
            else:
                between_full += " "

        for i in range(0,length,150):
            file.write(aligned_q[i:i+150] + "\n")
            file.write(between_full[i:i+150] + "\n")
            file.write(aligned_db[i:i+150] + "\n" + "\n")
            print(aligned_q[i:i+150])
            print(between_full[i:i+150])
            print(aligned_db[i:i+150] + "\n")

# just for q4, output alignments only to file, not on command line bc using subprocess, and that would take too long
if forQ4==True and output==True:
    with open("4top15_alignmentsT15.txt","a") as file:
        # print("in top15_alignments.txt")
        between_full = ""
        last_output = -1

        file.write(f"Local alignment with score of {score} for ({s_query}, {e_query}, {s_db}, {e_db}, {score})" + "\n")

        for i in range(length):
            if aligned_q[i] == aligned_db[i]:
                between_full += "|"
            else:
                between_full += " "

        for i in range(0,length,150):
            file.write(aligned_q[i:i+150] + "\n")
            file.write(between_full[i:i+150] + "\n")
            file.write(aligned_db[i:i+150] + "\n" + "\n")
            # print(aligned_q[i:i+150])
            # print(between_full[i:i+150])
            # print(aligned_db[i:i+150] + "\n")

