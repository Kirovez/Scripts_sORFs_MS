
from collections import defaultdict
from Bio import SeqIO
__version__ = 2.0

class Feature():
    """
    class represents a single gff line as an object
    """
    def sortedCoordinates(self, coord1, coord2):
        coord1 = int(coord1)
        coord2 = int(coord2)

        return [min(coord1,coord2), max(coord1,coord2)]

    def __init__(self, name, chromosome_name, strand, RNAname = ""):
        self.name = name
        self.chromosome_name = chromosome_name
        self.strand = strand
        self.exons = []
        self.features = defaultdict(list)  # feature name : [coordinates]
        self.introns = []
        self.RNAname = RNAname

    def addFeature(self, feature, coordinate1, coordinate2):
        if feature in self.features:
            self.features[feature].append(self.sortedCoordinates(coordinate1,coordinate2))
        else:
            self.features[feature] = [self.sortedCoordinates(coordinate1,coordinate2)]
        if feature == "exon":
            self.exons.append(self.sortedCoordinates(coordinate1,coordinate2))

    def determineIntrons(self):
        for i, exon in enumerate(self.getExons()):
            if i != 0: # start from second exon
                intron = [self.exons[i-1][1] + 1, exon[0] - 1]
                self.introns.append([min(intron),max(intron)])

    def getIntrons(self):
        self.determineIntrons()
        return sorted(self.introns, key=lambda x:x[0])

    def getExons(self):
        self.exons = [[min(i), max(i)] for i in self.exons]
        self.exons = sorted(self.exons, key=lambda x: x[0])
        return self.exons

    def getFeature(self, feature):
       return sorted(self.features[feature], key=lambda x:x[0])




class GffParser():

    filename = ""

    feature_dictionary = defaultdict(Feature)

    def __init__(self, filename):
        self.filename = filename
        self.readGff()
        print("Parser object has been created")

    def totalGenes(self):
        return len(self.feature_dictionary)

    def parse(self):
        feature_list = [self.feature_dictionary[i] for i in self.feature_dictionary]
        return feature_list

    def readGff(self):
        with open(self.filename) as gff_file:
            for lines in gff_file:
                if not lines.startswith("##"):
                    sp = lines.rstrip().split("\t")
                    # split gene name (e.g. ID=Pp3c1_10000;Name=Pp3c1_10000;ancestorIdentifier=Pp3c1_10000.v3.1).
                    # !!!Change if it is not necessary
                    feature_name = sp[8].split(";")[0].split("=")[1]
                    chromosome_name = sp[0]
                    strand = sp[6]
                    feature = sp[2]
                    coordinate1 = sp[3]
                    coordinate2 = sp[4]

                    if feature == "mRNA":
                        for_dic_name = feature_name + "_" + chromosome_name + "_" + strand
                        mRNAname = sp[8].split(";")[1].split("=")[1]
                        self.feature_dictionary[for_dic_name] = Feature(feature_name, chromosome_name, strand, RNAname=mRNAname)

                    # if you need to extract features per gene delete following lines and change mRNA to gene above
                    if feature != "gene":
                        self.feature_dictionary[for_dic_name].addFeature(feature,coordinate1,coordinate2)

    def writeIntrons(self, filenameout):
        """
        function to extract introns from gff3 file of physco
        :param filenameout: output filename
        :return: create file with intron coordinates where a line contains: chromosome name, start, end, strand
        """
        printed = {}
        with open(filenameout, "w") as out:
            for features in self.feature_dictionary:
                intro = self.feature_dictionary[features].getIntrons()
                print(self.feature_dictionary[features].introns)
                for introns in intro:
                    if "-".join([str(i) for i in introns]) not in printed:
                        out.write(self.feature_dictionary[features].chromosome_name + "\t"
                                  + str(introns[0]) + "\t" + str(introns[1]) + "\t" + self.feature_dictionary[features].strand + "\n")
                        printed["-".join([str(i) for i in introns])] = 0
    def checkUTRs(self):
        """
        :return: dictionary with names of RNAs possesing both UTRs
        """
        d = {}
        for feat in self.feature_dictionary:
            if "three_prime_UTR" in self.feature_dictionary[feat].features and \
                            "five_prime_UTR" in self.feature_dictionary[feat].features:
                d[self.feature_dictionary[feat].RNAname] = 0

        print(d)
        return d


    def writeFeature(self, feat, fileout):
        printed = {}
        with open (fileout, "w") as out:
            for genes in self.feature_dictionary:
                if self.feature_dictionary[genes].getFeature(feat):
                    for features in self.feature_dictionary[genes].getFeature(feat):
                        if "-".join([str(i) for i in features]) not in printed:
                            out.write(self.feature_dictionary[genes].chromosome_name + "\t"
                                      + str(features[0]) + "\t" + str(features[1]) + "\t" + self.feature_dictionary[genes].strand +"\n")
                            printed["-".join([str(i) for i in features])] = 0

    def getORFCoordinates(self, annotation=False):
        # find pacs with homology to Athaliana
        pacs_wih_annotation = {}
        if annotation:
            with open("Pp3_gene_annotation.txt") as ann:
                for lines in ann:
                    sp = lines.rstrip().split("\t")
                    if sp[-1] != "":
                        pacs_wih_annotation["PAC:" + sp[0]] = 0

        with open("ORF_Pp3.bed", "w") as out:
            for genes in self.feature_dictionary:
                #print(genes)
                if annotation and genes.split("_")[0] in pacs_wih_annotation:
                    if self.feature_dictionary[genes].getFeature("CDS"):
                        cds_coord = self.feature_dictionary[genes].getFeature("CDS")
                        out.write(genes + "\t" + self.feature_dictionary[genes].chromosome_name + "\t" +
                                  str(cds_coord[0][0]) + "\t" + str(cds_coord[-1][1]) + "\t" + "CDS" + "\t" +
                                  self.feature_dictionary[genes].strand + "\n")


"""
usage:

to extract introns:

GP = GffParser("Ppatens_318_v3.3.gene_exons.gff3")
GP.writeIntrons("Ppatens_318_v3.3_introns.bed")

"""

#print(len(GP.feature_dictionary))
#print(len(GP.checkUTRs()))
#GP.getORFCoordinates(annotation=True)
#

#GP.writeFeature("CDS", "Ppatens_318_v3.3_CDS.bed")

if __name__ == "main":
    GP = GffParser("Ppatens_318_v3.3.gene_exons.gff3")
    GP.writeIntrons("Ppatens_318_v3.3_introns.bed")
    GP.writeFeature("CDS", "Ppatens_318_v3.3_CDS.bed")


