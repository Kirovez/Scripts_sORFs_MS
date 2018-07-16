
from Bio import SeqIO
from Bio.Seq import Seq

class sORFCompletenessChecker():

    """
    class takes genome fasta file, fasta file protein translated sORFs and table from oopBLASTv_forsORF.py \
            return table with following column Query\tHit\tStart\tStop\tLength


    THIS VERSION TAKES INTO ACCOUNT WHETHER START CODON INSIGHT ALIGNMENT


    """
    def __init__(self, sORF_fasta, full_hit_fasta,bed_blast_out):

        self.sORF_fasta, self.full_hit_fasta, self.bed_blast_out = \
            sORF_fasta, full_hit_fasta,bed_blast_out

        self.outputFinalTab = open(self.sORF_fasta + "_vs_" + self.getfullhitfasta() + ".compcheckNEW", "w")
        self.bed_file_non_truncated = open(self.sORF_fasta + "_vs_" + self.getfullhitfasta() + ".BED_non_truncated", "w")

        self.HitIndex = SeqIO.index(self.full_hit_fasta, "fasta")
        self.QueryIndex = SeqIO.index(self.sORF_fasta, "fasta")

        self.runPip()

    def __doc__(self):
        print(" class takes genome fasta file, fasta file protein translated sORFs and table from oopBLASTv_forsORF.py \
            return table with following column Query\tHit\tStart\tStop\tLength")

    def getfullhitfasta(self):
        if "/" in self.full_hit_fasta:
            return self.full_hit_fasta.split("/")[-1]
        else:
            return self.full_hit_fasta

    def runPip(self):
        self.readBed()

    # read Bed file and return line one by one to function for getLine
    def readBed(self):
        with open(self.bed_blast_out) as bedfile:
            for line in bedfile:
                #print(line)
                self.getLine(line)

    def writeInBed(self, lineFromFinaltab):
        sp = lineFromFinaltab.rstrip().split("\t")
        if sp[-2] != "0":

            strand = ""
            if "-" in sp[-9]:
                strand = "-"
            else:
                strand = "+"

            if strand == "-":
                start = sp[-4]
                end = min([int(i) for i in [sp[-3], sp[-5]] if i != "-"])
            else:
                start = sp[-5]
                end = min([int(i) for i in [sp[-3], sp[-4]] if i != "-"])


            self.bed_file_non_truncated.write("\t".join([sp[1], str(start), str(end), sp[0], "0", strand]) + "\n")

    def getLine(self, line_from_bed):
        if "Query" not in line_from_bed: #not first line

            sp = line_from_bed.rstrip().split("\t")
            QueryId = sp[0]
            HitID = sp[1]

            StartHit, EndHit = int(sp[2]), int(sp[3])

            HitStrand = sp[-3]

            from_checking = self.checkCompletenesQuery(StartHit, EndHit, HitID, HitStrand)

            self.outputFinalTab.write(line_from_bed.rstrip() + "\t" + from_checking + "\t" + str(len(self.QueryIndex[QueryId])-1) + "\n")

            #self.writeInBed(line_from_bed.rstrip() + "\t" + from_checking + "\t" + str(len(self.QueryIndex[QueryId])-1) + "\n")



    def checkCompletenesQuery(self, coord1, coord2, hitseqID, strand):
        plusStart = "ATG"
        minusStart = "CAT"

        StartCodonCoord = "-"
        StopCodonCoord = "-"

        length = 0

        prematureStopCodon = "-"
        prematureStart = "-"

        plusStop = ["TAA", "TGA", "TAG"]
        minusStop = ["TTA", "TCA", "CTA"]

        seq = str(self.HitIndex[hitseqID].seq)


        startCur = plusStart
        stopCur = plusStop

        if "-" in strand :
            startCur = minusStop
            stopCur = minusStart



        # startCodonSearch and StopCodonCoord (if '-' strand) upstream
        ##chech start codon in alignment
        seq_tr = seq[coord1:coord2]
        if "-" in strand and minusStart in seq_tr:
            prematureStart = str(seq_tr.index(minusStart))
        elif "-" not in strand and plusStart in seq_tr:
            prematureStart = str(seq_tr.index(plusStart))


        init = coord1
        for stcodon in range((len(seq[0:coord1])//3) + 1):
            codon = seq[init:init+3]

            if  codon in startCur:

                if "-" in strand:
                    StopCodonCoord = str(init)
                else:
                    StartCodonCoord = str(init)
                break
            elif codon in stopCur and "-" not in strand:
                prematureStopCodon = str(init)
                break

            init -= 3


        # startCodonSearch (if '-' strand) and StopCodonCoord downstream
        initStop = coord2
        for stopcodon in range((len(seq[coord2:])//3) + 1):
            codon = seq[initStop-3:initStop]


            if codon in stopCur:

                if "-" in strand:
                    StartCodonCoord = str(initStop)
                else:
                    StopCodonCoord = str(initStop)
                break
            elif codon in stopCur and "-" in strand:
                prematureStopCodon = str(initStop)
                break

            initStop += 3

        if prematureStopCodon == "-":
            # check any stop codons in blast hit sequence part
            if "-" in strand:
                seq_tr = str(Seq(seq[coord1:coord2]).reverse_complement().translate())[:-1:-1]
            else:
                seq_tr = str(Seq(seq[coord1:coord2]).translate())[:-1]

            if "*" in seq_tr[:-1]:
                StopCodonCoord = str((seq_tr.index("*")*3) + coord1)


         ### calculate lengths ####
        length = 0
        # if stop codon before start length = 0
        if prematureStopCodon != "-":
            length = 0

        elif StartCodonCoord != "-":
            if StopCodonCoord != "-":
                length = abs(int(StopCodonCoord) - int(StartCodonCoord))

        # no start but premature start exists
        elif prematureStart != "-":
            if StopCodonCoord != "-":
                length = abs(int(StopCodonCoord) - int(prematureStart))



        toRet = [StartCodonCoord, StopCodonCoord, prematureStopCodon, str(int(length/3)-1)]

        return "\t".join(toRet)


#translated_sORF = r"/home/ilia/sORF/dnds/Ka_Ks_codeml/CDS_trans.fasta"
"""
path = r"/home/ilia/contigs/"
species_fasta = ["Zmays_284_Ensembl-18_2010-01-MaizeSequence.transcript.fa", "Vcarteri_317_v2.1.transcript.fa",
                 "Spolyrhiza_290_v2.transcript.fa",
                 "Smoellendorffii_91_v1.0.transcript.fa",
                 "Sfallax_310_v0.5.transcript.fa",
                 "Osativa_323_v7.0.transcript.fa",
                 "Mpolymorpha_320_v3.1.transcript.fa",
                 "Creinhardtii_281_v5.5.transcript.fa",
                 "Cpuppureus_trinity",
                 "Athaliana_167_TAIR10.transcript.fa"]


for i in species_fasta:
    query = "translated_FINAL_sORF_SELECTED.fa"
    species_file = path + i
"""

sORFCompletenessChecker("translated_FINAL_sORF_SELECTED.fa", "/home/ilia/RNA_seq_reads_moss/TopHat/Ppatens_318_v3.fa",
                            "blast_translated_sORFs_vs_Ppatens_genome.xml.bed")
