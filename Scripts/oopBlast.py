from Bio import SeqIO
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SearchIO


class BlastParser():
    def __init__(self, in_file, DB_file, Similarity=70,
                 Evalue=1e-5, Coverage=70, BLAST = True,
                 DB_build = True, outfmt=5, wordsize=12,
                 blastpath= r'C:\Users\ikirov.CLO\Desktop\blast-2.2.31+\bin\blastn.exe',blastbuild=r"C:\Users\ikirov.CLO\Desktop\blast-2.2.31+\bin\makeblastdb.exe"):
        self.in_file = in_file
        self.DB_file = DB_file
        self.Similarity = Similarity
        self.Evalue = Evalue
        self.Coverage = Coverage
        self.BLAST = BLAST
        self.DB_build = DB_build
        self.file_name = '{0}_vs_{1}'.format(in_file,DB_file)
        self.oufmt = outfmt
        self.wordsize = wordsize
        self.blastpath = blastpath
        self.blastbuild = blastbuild

    def runblast(self):
        if self.DB_build:
            cmd = self.blastbuild + ' -in %s -dbtype nucl' % self.DB_file
            os.system(cmd)
        proga = self.blastpath
        blast = NcbiblastnCommandline(proga, query=self.in_file, db=self.DB_file, out=self.file_name, outfmt=self.oufmt, word_size=self.wordsize,
                                      evalue=0.001)  # strand = 'plus'
        stdout, stderr = blast()

    #  function from stack http://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
    def merge_intervals(self,intervals):
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
        #print(added_length)
        merged.append(added_length)
        return merged

    def parseBlastXml(self):
        query_indexed_db = SeqIO.index(self.in_file,'fasta')
        count = 0
        with open('BED', 'w') as hit:
            file2 = open(self.file_name)
            s1 = SearchIO.parse(file2, 'blast-xml')
            count_hits = 0

            #  query
            for recor in s1:
                number_of_hits = 0

                # hit sequence
                for HSP in recor:
                    len_query = len(str(query_indexed_db[recor.id].seq))
                    hs = HSP.hsps
                    hit_ranges = []

                    # HSP of the hit sequence
                    for u in hs:
                        if u.evalue <= self.Evalue and u.ident_num * 100 / len(u.query) >= self.Similarity:
                            e_value = u.evalue
                            count_hits += 1

                            hit_ranges.append([u.hit_range[0], u.hit_range[1]])

                    if hit_ranges:
                        # Coverage calculation based on the merged intervals
                        added_values = self.merge_intervals(hit_ranges)[-1]
                        for i, intervals in enumerate(self.merge_intervals(hit_ranges)[0:-1]):
                            if (abs(intervals[1] - intervals[0]) - added_values[i]) * 100 / len_query >= self.Coverage:
                                hit.write(
                                    HSP.id + '\t' + recor.id + '\t' + str(sorted(intervals)[0]) + '\t' + str(sorted(intervals)[1]) + '\n')
                                count += 1
        print(count)

#Bl = BlastParser("Tyr-FISH_sequences.fasta","Chromosome_seq_reconstr_351.fasta", DB_build=True)

#Bl.runblast()
#Bl.parseBlastXml()