from Bio import SeqIO
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SearchIO

in_file =  'Tyr-FISH_sequences.fasta'
DB_file = 'Chromosome_seq_reconstr_351.fasta'
Similarity = 80
Evalue = 1e-15
no_hits_ind = 0
query_list = []
Coverage = 85
query_indexed_db  = SeqIO.index(in_file,'fasta')
print('Number of query sequences:' + str(len(query_indexed_db)))

file_name = '{0}_vs_{1}'.format(in_file,DB_file)
#with open(file_name, 'w') as file:
    #print('File opening')
#cmd = r'C:\Users\ikirov.CLO\Desktop\blast-2.2.31+\bin\makeblastdb.exe -in %s -dbtype nucl' % DB_file
#os.system(cmd)
#proga = r'C:\Users\ikirov.CLO\Desktop\blast-2.2.31+\bin\blastn.exe'
#blast = NcbiblastnCommandline(proga, query=in_file, db=DB_file, out=file_name, outfmt=5, word_size = 12, evalue = 0.001)#strand = 'plus'
#stdout, stderr = blast()

count=0

#  function from stack http://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
def merge_intervals(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    added_length = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
            added_length.append(0)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
                added_length[-1] = 0
            elif higher[0] - lower[1] <= 3000:
                upper_bound = max(lower[1], higher[1])
                added_length[-1] = abs(higher[0] - lower[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
            else:
                merged.append(higher)
                added_length.append(0)
    print(added_length)
    merged.append(added_length)
    return merged



with open('BED', 'w') as hit:
    file2 = open(file_name)
    s1 = SearchIO.parse(file2, 'blast-xml')
    count_hits = 0

    #  query
    for recor in s1:
        number_of_hits=0

        # hit sequence
        for HSP in recor:
            len_query = len(str(query_indexed_db[recor.id].seq))
            indicator = 0
            hs = HSP.hsps
            hit_ranges = []

            # HSP of the hit sequence
            for u in hs:
                if u.evalue <= Evalue and u.ident_num * 100 / len(u.query) >= Similarity:
                    e_value = u.evalue
                    indicator += 1
                    count_hits += 1

                    hit_ranges.append([u.hit_range[0],u.hit_range[1]])

            # Coverage calculation based on the merged intervals
            added_values = merge_intervals(hit_ranges)[-1]
            for i, int in enumerate(merge_intervals(hit_ranges)[0:-1]):
                if (abs(int[1] - int[0]) - added_values[i])*100/len_query >= Coverage and indicator >= 1:
                    hit.write(HSP.id + '\t' + recor.id  + '\t' + str(sorted(int)[0]) + '\t' + str(sorted(int)[1]) + '\n')
                    count+=1
print(count)