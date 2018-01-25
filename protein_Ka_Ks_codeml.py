
"""
Script run dn/ds test using protein queries and genome sequence;
It depends on:

1. perl, python, biopython
2. pal2nal perl script, available here https://github.com/HajkD/orthologr/tree/master/inst/pal2nal/pal2nal.v14.
 !!!set 'pal2nal' variable to run pal2nal script e.g. pal2nal=r"perl ./pal2nal.v14/pal2nal.pl" (default)

3.  blast. Please, if it is required, change 'blast' variable in class e.g.  blast=r"tblastn" (default)
4. codeml. .ctl file can be found here
5. clustalo.
6. codeml Parser module (see at codemlParser directory)

!! change variable query_seq_ind to the path to the file with sORF nucleotide sequences
e.g. (default)
query_seq_ind = SeqIO.index("/home/ilia/sORF/sORFfinder2/moss/FINAL_sORF_SELECTED.fa", "fasta")

7. oopBLASTv_forsORF.py script. Available in this repository
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq
from codemlParser.CodemlParser import CodemlParser
from oopBLASTv_forsORF import BlastParser
import argparse
import os.path



class runKaKsCalculatio():
    def __init__(self, fasta_query, fasta_hit, codeml="codeml", clustalo="clustalo", pal2nal=r"perl ./pal2nal.v14/pal2nal.pl",
                 blast=r"tblastn",blast_run = True, parse_blat=True, threads = 4, makedb=False):
        self.codeml = codeml
        self.clustal = clustalo
        self.pal2nal = pal2nal
        self.blast=blast
        self.fasta_query = fasta_query
        self.fasta_hit = fasta_hit
        self.blast_run = blast_run
        self.parse_blast = parse_blat
        self.threads = threads
        self.makedb = makedb

        self.Bl = BlastParser(self.fasta_query, self.fasta_hit, DB_build=self.makedb,blastpath=self.blast, BLAST=self.blast_run, threads=self.threads,
                              Similarity=30, Coverage = 30, wordsize=3)

        self.bedBlast = self.Bl.bedFileSbj
        self.ind = SeqIO.index(self.fasta_query, "fasta")
        #self.fall_cnt_bl2seq_eval = 0 # Number of sequences with too high evalue from bl2seq. They are not passed through
        self.fall_indentical_sequences = 0 # Number of identical sequences as revealed by bl2seq
        self.fall_short_length_cnt = 0
        self.pal2nal_outfile = "temp_for_codeml.phylip" # output file from PAL2NAL whic will be used for codeml
        self.outputFileTable = self.getOutFilename()
        self.finalTable = open(self.outputFileTable, "w")
        self.finalTable.write("Sequence Query\tSequence Hit\tt\tS\tN\tomega\tdN\tdS\tLRT\tSignificance\n")
        self.analysed_dnds_cnt = 0
        self.list_queries_for_codeml = []

    def runPal2Nal(self,query_id,sbj_id,bl2seq_output):
        os.system(self.pal2nal + " {0} {1} {2} -nogap -output paml > {3}".format("tmp_file_palnal.aln", "tmp_blast1.fasta", "tmp_blast3.fasta", self.pal2nal_outfile))
        self.runCodeml(query_id,sbj_id,bl2seq_output)

    def getOutFilename(self):
        ifile = ""
        dbfile = ""
        if "/" in self.fasta_query:
            ifile = self.fasta_query.split("/")[-1]
        if "/" in self.fasta_hit:
            dbfile = self.fasta_hit.split("/")[-1]

        return ifile + "_vs_" + dbfile + "_dnds.out"


    def getFilePal2Nal(self,query_id,query_seq,sbj_id,sbjct_seq,bl2seq_output):
        print("*** run PAL2NAL  ***")
        with open("tmp_file_palnal.aln", "w") as outfile:
            outfile.write(">query" + "\n" + query_seq + "\n")
            outfile.write(">hit" + "\n" + sbjct_seq + "\n")
        self.runPal2Nal(query_id,sbj_id,bl2seq_output)



    def parseBedfromBl(self):
        query_seq_ind = SeqIO.index("/home/ilia/sORF/sORFfinder2/moss/FINAL_sORF_SELECTED.fa", "fasta")
        hit_ind = SeqIO.index(self.fasta_hit, "fasta")

        with open(self.Bl.bedFileSbj) as sbj:
            outfile = ""
            for i,lines in enumerate(sbj):
                if i != 0:
                    print(i)
                    sp = lines.split("\t")
                    #sORF_seq = str(query_seq_ind[sp[0]].seq)
                    protein_sORFs = str(self.ind[sp[0]].seq)

                    if len(protein_sORFs) >= 25:

                        if "-" in sp[9]:
                            nucl_sbj = str(Seq(str(hit_ind[sp[1]].seq)[int(sp[3]):int(sp[4])]).reverse_complement())
                        else:
                            nucl_sbj = str(hit_ind[sp[1]].seq)[int(sp[3]):int(sp[4])]

                        nucl_query = str(query_seq_ind[sp[0]].seq)[int(sp[-2])*3:int(sp[-1])*3]

                        outfile = sp[1] + "_vs_" + sp[0]

                        with open("tmp_blast1.fasta", "w") as tmp1, \
                                open("tmp_blast2.fasta", "w") as tmp2,\
                                open("tmp_blast3.fasta", "w") as tmp3:



                            tmp1.write(">" + "query" + "\n" + nucl_query  + "\n")
                            tmp2.write(">" + "hit" + "\n" + protein_sORFs + "\n")
                            tmp3.write(">" + "hit" + "\n" + nucl_sbj  + "\n")
                        print(outfile)



                        if nucl_query != nucl_sbj:
                            self.list_queries_for_codeml.append(sp[0])
                            self.getFilePal2Nal(sp[0],sp[7],sp[1],sp[6],outfile)

                        else:
                            self.fall_indentical_sequences+=1

                        #try:
                            #print("*** run bl2seq  ***")
                            #os.system("bl2seq -i {0} -j {1} -p blastx -o {2}".format("tmp_blast1.fasta","tmp_blast2.fasta",outfile))
                            #self.parseblast2seqOut(outfile)
                        #except:
                            #print("ERROR: bl2seq")

                    else:
                        self.fall_short_length_cnt+=1


    def run_pipeline(self):
        print("BLAST is running...")
        self.runBlast()
        print("BLAST output parsing...")
        self.parseBedfromBl()

    def runCodeml(self,query_id,sbj_id,bl2seq_output):
        # generate tree
        self.analysed_dnds_cnt+=1
        print("*** run CODEML  ***")
        sORF_hit = bl2seq_output.split("_vs_")
        lis= [sORF_hit[1],sORF_hit[0]]
        codemlP = CodemlParser(query_id,sbj_id)
        lin_to_write= "\t".join(lis) + "\t" + codemlP.getLineforFinalTable()
        self.finalTable.write(lin_to_write)

    def runBlast(self):
        #os.system(self.blast + " {0} {1} -t=dna -q=dna -stepSize=1 -minScore=0 -minIdentity=0 {2}".format(self.fasta_hit,self.fasta_query,self.out_blat_file))
        if self.blast_run:
            self.Bl.runblast()
            self.Bl.parseBlastXml()
        elif self.parse_blast:
            self.Bl.parseBlastXml()

if __name__ == '__main__':
    makedb = False
    blastrun = True
    parser = argparse.ArgumentParser(description = 'Run dn/ds test using protein queries')
    parser.add_argument("infile", help='name of sORF fasta file')
    parser.add_argument("dbfile", help='name of target genome fasta file')
    parser.add_argument("--threads", help='number of threads',default=2)
    parser.add_argument("--makedb", help="Make database for blast or it already exists", default="F")
    parser.add_argument("--blast", help="Run blst or not", default="T")


    args = parser.parse_args()
    if args.makedb != "F":
        makedb = True

    if args.blast != "T":
        blastrun = False

    KKC = runKaKsCalculatio(args.infile, args.dbfile,blast_run=blastrun, threads=args.threads, makedb=makedb)
    KKC.run_pipeline()

    print("+++++ BLAST step: +++++++")
    print("Number of query sequences: " + str(len(KKC.Bl.all_input_sequences)))
    print("Number of uniquely named query sequences: " + str(len(set(KKC.Bl.all_input_sequences))))
    print("Number of sequences with similarity: " + str(len(set(KKC.Bl.sequences_with_hit))))
    nohit = [seq for seq in KKC.Bl.all_input_sequences if seq not in KKC.Bl.sequences_with_hit]
    print("Sequences with no similarity to the subjects: " + str(len(nohit)))

    print("+++++ CODEML step +++++ ")
    print("number of sORFs used for analysis: ", len(set(KKC.list_queries_for_codeml)))
    print("identical sequences ", KKC.fall_indentical_sequences)
    print("number of events annalysed for dn\ds ", KKC.analysed_dnds_cnt)
    print("number of too short sequence which were not used for dn\ds calculation ",KKC.fall_short_length_cnt)

    os.system("rm -f results*")

