from Bio import SeqIO
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SearchIO


class BlastParser():
    def __init__(self, in_file, DB_file, Similarity=30,
                 Evalue=1e-5, Coverage=40, BLAST=True,
                 DB_build=True, outfmt=5, wordsize=12,
                 blastpath=r'blastn', blastbuild=r"makeblastdb", dbtype="nucl", threads = 4):
        self.in_file = in_file
        self.DB_file = DB_file
        self.Similarity = Similarity
        self.Evalue = Evalue
        self.Coverage = Coverage
        self.BLAST = BLAST
        self.DB_build = DB_build
        self.file_name = self.getOutFilename()
        self.oufmt = outfmt
        self.wordsize = wordsize
        self.blastpath = blastpath
        self.blastbuild = blastbuild
        self.dbtype = dbtype
        self.unique = []
        self.bedFileSbj = self.file_name + ".bed"
        self.threads = threads
        self.all_input_sequences = []
        self.sequences_with_hit = []

    def getOutFilename(self):
        ifile = ""
        dbfile = ""
        if "/" in self.in_file:
            ifile = self.in_file.split("/")[-1]
        if "/" in self.DB_file:
            dbfile = self.DB_file.split("/")[-1]

        return ifile + "_vs_" + dbfile
    def runblast(self):
        if self.DB_build:
            cmd = self.blastbuild + ' -in {0} -dbtype {1}'.format(self.DB_file, self.dbtype)
            os.system(cmd)
        proga = self.blastpath
        blast = NcbiblastnCommandline(proga, query=self.in_file, db=self.DB_file, out=self.file_name, outfmt=self.oufmt,
                                      word_size=self.wordsize,
                                      evalue=0.001,num_threads=self.threads)  # strand = 'plus'
        stdout, stderr = blast()

    def checkCompletenesQuery(self, query0, query1):
        return abs(query0 - query1)

    def merge_intervals(self, intervals):
            if intervals:
                sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
                sorted_by_lower_bound.append([0 for i in range(len(intervals))])
                print(sorted_by_lower_bound)
                return sorted_by_lower_bound
            else:
                return intervals


    def getFileUnique(self):
        with open(self.file_name + "_list.unique", "w") as out:
            out.write("\n".join(set(self.unique)))


    def parseBlastXml(self,file_exo=""):
        if file_exo != "":
            self.file_name = file_exo
            self.bedFileSbj = file_exo+".bed"
        self.all_input_sequences = [seq.id for seq in SeqIO.parse(self.in_file,'fasta')]

        query_indexed_db = SeqIO.index(self.in_file, 'fasta')
        count = 0
        with open(self.bedFileSbj, 'w') as hit:
            hit.write("Query\tHit\tIdentity\tFrom\tTo\tEvalue\tCoverage\n")
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

                    # HSP of the hit sequence
                    for u in hs:
                        iden = u.ident_num * 100 / len(u.query)
                        if u.evalue <= self.Evalue and iden >= self.Similarity: # and self.checkCompletenesQuery(u.hit_range[0], u.hit_range[1]) >= (self.Coverage/100)*len_query:
                            count_hits += 1

                            # Coverage calculation based on the merged intervals

                            self.sequences_with_hit.append(recor.id)
                            self.unique.append(recor.id)
                            hit.write(
                                recor.id + '\t' +
                                HSP.id + '\t' +
                                str(u.hit_range[0]) + '\t' +
                                str(u.hit_range[1]) + '\t' +
                                str(u.evalue) + "\t" +
                                str(u.hit.seq)+ "\t" +
                                str(u.query.seq) + "\t" +
                                str(self.checkCompletenesQuery(u.hit_range[0], u.hit_range[1]))+ "\t" +
                                str(u.hit_strand) + "\t" +
                                str(u.query_range[0]) + "\t" +
                                str(u.query_range[1]) + '\n' )
                            count += 1





Bl = BlastParser(r"/home/ilia/sORF/sORFfinder2/moss/FINAL_sORF_SELECTED.fa", r"/home/ilia/SOLID_moss/Ppatens_318_v3_index.fa", DB_build=False)

#Bl.runblast()
Bl.parseBlastXml("blast_translated_sORFs_vs_Ppatens_genome.xml")
