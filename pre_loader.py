import pandas as pd
from Parser import nrps_parser
import os

class PreLoader:
	def __init__(self, path):
	################# required input files ###############
		print "#############"
		print "Root_Path: ", path
		print "###############"
		print "#############"
		print os.path.join(path,"files2Read/MIBiG_webapp_clean.csv")
		print "###############"
		self.MIBiG_DB1 = pd.read_csv(os.path.join(path,"files2Read/MIBiG_webapp_clean.csv"))
		self.NCBI_DB1 = pd.read_csv(os.path.join(path,"files2Read/NCBI_webapp_clean.csv"))
		self.MIBiG_DB2 = pd.read_csv(os.path.join(path, "files2Read/MIBiG_webapp_not_clean.csv"))
		self.NCBI_DB2 = pd.read_csv(os.path.join(path, "files2Read/NCBI_webapp_not_clean.csv"))
		self.NCBI_DB3 = pd.read_csv(os.path.join(path, "files2Read/NCBI_RefSeq_genomes_ftp.csv"))
        ############ required input instances ################
		self.nrps_parser = nrps_parser.NRPS_Parser()



