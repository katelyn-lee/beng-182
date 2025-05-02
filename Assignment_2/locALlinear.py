# Mock up of iBLAST tool at NCBI. modified version of locAL.py still using dynamic programming matrices, but in linear space. 
# Given a score cut-off, finds all local alignments of query to database with a minimum score threshold of T.
#  Outputs a preliminary hits table with all qualifying local alignment hits sorted in decreasing order. Then outputs a final hits table that 
#    pruned away all alignments that overlapepd with other hits on the database. 
#    Each hit is outputted on its own line as: begin_query, end_query, begin_db, end_db, score 
#    Finally, the top 15 highest scoring hits have their actual alignments outputted in top15_alignments.txt


# Algorithm:
#   locAL.py modified to use linear space with 2 column method on DP matrix. uses subprocess to align top 15 to out file in locAL.py itself

import sys
from typing import List, Dict, Tuple
import subprocess

# python3 locALlinear.py -q <query> -d <db> -m <match> -s <mismatch> -i <indel> -T <threshold_value>

# initialize variables
query_file, db_file = "", ""
query, db = "", ""
match, mismatch, indel = 1, 2, 2
T = 0

try: 
  if len(sys.argv) == 14:
      query_file = sys.argv[2]
      db_file = sys.argv[4]
      match = sys.argv[6]
      mismatch = sys.argv[8]
      indel = sys.argv[10]
      T = sys.argv[12]
except Exception as e:
    print(f"An unexpected error occurred: {e}")
  
with open(query_file,'rt') as qf, open(db_file,'rt') as dbf:
    next(qf)
    query = (qf.readline()).strip()
    next(dbf)
    db = (dbf.readline()).strip()

# run DP algorithm to get local alignment
high_scoring_alignments = []

def locAL_linear_space(match:int, mismatch:int, indel:int, q:str, db:str):
    global high_scoring_alignments, T

    m, n = len(q), len(db)

    # 2 cols for linear space alignment
    prev_col = [0] * (m + 1)
    curr_col = [0] * (m + 1)

    # track positions of curr pos's origin
    prev_starts = [] 
    for i in range(m+1):
        prev_starts.append((i,0))
    curr_starts = [(0, 0)] * (m + 1)  

    for j in range(1,n+1):
        print(f"Doing local alignment, currently {j/n}% of matrix done")
        curr_col[0] = 0
        curr_starts[0] = (0, j)  # Alignment starts fresh at top of column

        for i in range(1,m+1):
            match_val = prev_col[i - 1] + (match if q[i - 1] == db[j - 1] else (-mismatch))
            insert = prev_col[i] - indel
            delete = curr_col[i - 1] - indel 
            score  = max(match_val,insert,delete,0)

            curr_col[i] = score
            if score == 0:
                curr_starts[i] = (i,j)
            elif score == match_val:
                curr_starts[i] = prev_starts[i-1]
            elif score == delete:
                curr_starts[i] = curr_starts[i-1]
            elif score == insert:
                curr_starts[i] = prev_starts[i]

            # record high-scoring local alignments in high_scoring_alignments w  {begin_query, {end_query}, {begin_db}, {end_db}, {score}
            if score >= T:
                end_pos = (i-1,j-1)
                start_pos = curr_starts[i]
                high_scoring_alignments.append((start_pos[0], end_pos[0], start_pos[1], end_pos[1], score))
                print(f"adding to high_scoring_alignments.")

        # reset for next col in DP  
        prev_col, curr_col = curr_col, [0] * (m + 1)
        prev_starts, curr_starts = curr_starts, [(0, 0)] * (m + 1)
    return high_scoring_alignments

locAL_linear_space(match,mismatch,indel,query,db)

if len(high_scoring_alignments) == 0:
    print(f"There are no alignments with at least a score of {T}")

b = len(high_scoring_alignments)

# preliminary hits table:
high_scoring_alignments.sort(key=lambda x: x[-1], reverse=True)
print(f"Now have prelim hits table with {b} entries. first 100 here: {high_scoring_alignments[:100]}" + "\n" + "\n")

with open("prelim_hits.txt","w") as prelim_out_file:
    prelim_out_file.write(f"Preliminary hits table has {b} entries in format of (begin_query, end_query, begin_db, end_db, score):" + "\n")
    for entry in high_scoring_alignments:
        prelim_out_file.write(str(entry) + "\n")

# prune table for only distinct alignments. 
def is_50overlap(a,b):  #check for overlap on db subsequences only
    start_a, start_b = a[2], b[2]
    end_a, end_b = a[3], b[3]

    # calc overlap range
    overlap_start = max(start_a, start_b)
    overlap_end = min(end_a, end_b)

    if overlap_end < overlap_start:
        return False #bc no overlap

    overlap_length = overlap_end - overlap_start + 1
    db_sub_length = end_a - start_a + 1

    # print(f"overlap_length: {overlap_length}")
    # print(f"db_sub_length: {db_sub_length}")
    if (overlap_length/db_sub_length) >= 0.5:
        return True
    else:
        return False

final_hits = []

# for prelim hits, if theres even one overlap with kept hits, can't keep it, so break 
for hit in high_scoring_alignments:
    keep = True
    for kept in final_hits:
        if is_50overlap(kept,hit)==True:
            keep = False
            break
    if keep:
        final_hits.append(hit)

x = len(final_hits)
print(f"final hits table done with {x} hits, first 100 here: {final_hits[:100]}")

with open("final_hits.txt","w") as final_out_file:
    final_out_file.write(f"Final hits table has {x} entries in format of (begin_query, end_query, begin_db, end_db, score):" + "\n")
    for hit in final_hits:
        final_out_file.write(str(hit) + "\n")

top_15 = final_hits[:15]
print(f"top_15: {top_15}")

# generate locAL for top 15 hits on final hits table
def run_locAL(hit):
    start_query, end_query = hit[0], hit[1]
    start_db, end_db = hit[2], hit[3]
    # score = hit[4]

    # run locAL on substrings of query and db to get local alignment
    command1 = ["python3","locAL.py", "-q", query[start_query:end_query+1], "-d", db[start_db:end_db+1], "-m", "1", "-s", "2", "-i", "2", "-a", str(start_query), str(end_query), str(start_db), str(end_db)]
    try:
        # subprocess.run(command1,capture_output=True,text=True,check=True)
        subprocess.run(command1)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
    return 

# top_15 = [(0, 24, 888, 912, 25), (0, 13, 3589, 3602, 14), (6, 19, 4684, 4697, 14), (12, 24, 1634, 1646, 13), (0, 11, 888, 899, 12), (14, 24, 3561, 3571, 11), (15, 24, 2137, 2146, 10), (16, 24, 3, 11, 9), (16, 24, 211, 219, 9), (0, 8, 858, 866, 9), (16, 24, 2762, 2770, 9), (14, 21, 1178, 1185, 8), (0, 7, 3554, 3561, 8), (7, 20, 158, 173, 7), (2, 8, 623, 629, 7)]
n = 1
for hit in top_15:
    output = run_locAL(hit)
    print(f"aligned {n} hit: {hit}")
    n+=1


# TEST CASE HAS 11 INSERTIONS of query into db, but somewhat modified. one exact insertion though should have length 200
