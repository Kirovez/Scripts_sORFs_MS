from Bio.Phylo.PAML import codeml

"""
class takee two file after pal2nal and perform dn/ds analysis.
!!!! Change following variables:
# path to cnt file required for codeml
cntfile  = "codeml.ctl"

#working directory
cml.working_dir = "/home/ilia/sORF/dnds/Ka_Ks_codeml"

"""
class GKhit():
    def __init__(self, line):
        self.t, self.S, self.N, self.omega, self.dN, self.dS = "","","","","",""
        self.formatLine(line)


    def formatLine(self, line):
        s1 = line.rstrip().replace(" ","")
        st = s1.split("=")

        self.t = st[1][:-1]
        self.S = st[2][:-1]
        self.N = st[3][:-5]
        self.omega = st[4][:-2]
        self.dN = st[5][:-2]
        self.dS = st[6]

    def getForTable(self):
        return [self.t, self.S, self.N, self.omega, self.dN, self.dS]

class CodemlParser():
    def __init__(self,query_id, sbj_id):
        self.query_id, self.sbj_id = query_id, sbj_id
        self.fin_line = []
        self.runCodeml(self.query_id, self.sbj_id)


    def getLineforFinalTable(self):
        return "\t".join(self.fin_line) + "\n"

    def compareLnl(self,lns):
        tocheck = round(2 * (lns[0] - lns[1]),2)
        if tocheck > 6.635:
            return "{}\t**".format(tocheck)
        elif 2*(lns[0] - lns[1]) > 3.841:
            return "{}\t*".format(tocheck)
        else:
            return "{}\t-".format(tocheck)


    def readCodemlOutFiles(self,outfile):
        lnllist = []
        final_line = []
        final_line.append(self.query_id)
        final_line.append(self.sbj_id)
        ind = True
        with open(outfile) as outfile_1, \
                open(outfile.replace("XXXXX1","XXXXX0")) as outfile_0:
            for files in [outfile_0, outfile_1]:
                for lines in files:
                    if lines.startswith("lnL"):
                        lnl = float(lines.rstrip().replace(" ","").split("=")[-1])
                        lnllist.append(lnl)

                    elif lines.startswith("t=") and ind:
                        GK = GKhit(lines)
                        self.fin_line+=GK.getForTable()
                        ind = False # for both models the line with dnds is similar therefore it is not necessary to parse it again

            self.fin_line.append(self.compareLnl(lnllist))



    def runCodeml(self,query_id, sbj_id):
        # generate tree
        cntfile  = "codeml.ctl"

        cml = codeml.Codeml()
        cml.read_ctl_file(cntfile)
        #cml.alignment = "temp.storage"
        #cml.tree = "species.tree"

        cml.working_dir = "/home/ilia/sORF/dnds/Ka_Ks_codeml"

        with open(cml.working_dir + "/species.tree", "w") as outfile:
            outfile.write("(query,hit);")
        to_parse = []
        out_file = cml.working_dir + "/results_{0}_vs_{1}_omegaXXXXX.out".format(query_id, sbj_id)
        for i in [0,1]:
            if i == 0:
                out_file= out_file.replace("XXXXX", "XXXXX0")
                cml.out_file = out_file
                cml.set_options(fix_omega=0)
                cml.run(verbose = False)
            else:
                out_file = out_file.replace("XXXXX0", "XXXXX1")
                cml.out_file = out_file
                cml.set_options(fix_omega=1)
                cml.run(verbose=False)

        self.readCodemlOutFiles(out_file)

