import os
translated_sORF = r"translated_FINAL_sORF_SELECTED.fa"
path = r"/home/ilia/contigs/"
species_fasta = ["Zmays_284_Ensembl-18_2010-01-MaizeSequence.transcript.fa",
                 "Vcarteri_317_v2.1.transcript.fa",
                 "Spolyrhiza_290_v2.transcript.fa",
                 "Smoellendorffii_91_v1.0.transcript.fa",
                 "Sfallax_310_v0.5.transcript.fa",
                 "Osativa_323_v7.0.transcript.fa",
                 "Mpolymorpha_320_v3.1.transcript.fa",
                 "Creinhardtii_281_v5.5.transcript.fa",
                 "Cpuppureus_trinity",
                 "Athaliana_167_TAIR10.transcript.fa"]

for i in species_fasta:
    species_file = path + i
    cmd = "nohup python3 protein_Ka_Ks_codeml.py {0} {1} --blast T --threads 10 --makedb F".format(translated_sORF, species_file)
    os.system(cmd)

