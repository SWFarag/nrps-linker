import pandas as pd
from Parser import nrps_parser
from Design import nrps_design
import os
import dropbox

class PreLoader:
	def __init__(self, path):
	################# required input files ###############
		print "#############"
		print "Root_Path: ", path
		print "###############"
		print "#############"
		print os.path.join(path,"files2Read/MIBiG_webapp_clean.csv")
		print "###############"
		self.MIBiG_DB1 = pd.read_csv(os.path.join(path,"files2Read/MIBiG_webapp_clean3.csv"))
		self.NCBI_DB1 = pd.read_csv(os.path.join(path,"files2Read/NCBI_webapp_clean2.csv"))

		self.MIBiG_DB2 = pd.read_csv(os.path.join(path, "files2Read/MIBiG_webapp_not_clean3.csv"))
		self.NCBI_DB2 = pd.read_csv(os.path.join(path, "files2Read/NCBI_webapp_not_clean2.csv"))

		self.NCBI_DB3 = pd.read_csv(os.path.join(path, "files2Read/NCBI_RefSeq_genomes_ftp2.csv"))

		self.NCBI_DB4 = pd.read_csv(os.path.join(path, "files2Read/NCBI_linkers_with_domain_pos2.csv"))
		self.MIBiG_DB5 = pd.read_csv(os.path.join(path, "files2Read/MIBiG_linkers_with_domain_pos3.csv"))

		self.NCBI_HTML = os.path.join(path, "files2Read/NCBI_HTML.html")

		#self.peptides_ids_xmlFileNames = pd.read_csv(os.path.join(path, "files2Read/Nourine_info/peptides_ids_xmlFileNames.csv"))
		self.peptides_ids_graphs = pd.read_csv(os.path.join(path, "files2Read/Nourine_info/peptides_ids_graphs.csv"))
		self.peptides_ids_xmlFileNames = pd.read_csv(os.path.join(path, "files2Read/Nourine_info/design_templates_live.csv"))

		self.templates_info= pd.read_csv(os.path.join(path, "files2Read/templates/templates_info.csv"))
        ############ required input instances ################
		self.nrps_parser = nrps_parser.NRPS_Parser()
		self.nrps_design = nrps_design.NRP_Design()
		self.dbx = dropbox.Dropbox('oXa01q2Xf-4AAAAAAAAFat8UzDzV7AkX-cBrN8coWKDcFK-ANG9HIF8P8g-MmdJy')

	def tohtml_design(self, df, tableId):
		cols = df.columns
		s = '''<table id="''' + tableId + '''" class="display table table-striped table-bordered" width="100%">'''
		s = s + '''
	        <col width="220">
	        <col width="220">
	        <col width="220">
	        <col width="220">
	        <col width="220">
	        '''
		s = s + '''
	        <thead>
	            <tr>'''
		temp = '''
	            <th>Index</th>
	        '''
		s = s + temp
		for col in cols:
			temp = '''
	                <th>''' + col + '''</th>
	            '''
			s = s + temp
		end_header = '''</tr>
	        </thead>'''
		s = s + end_header

		s = s + '''
	        <tfoot>
	            <tr>'''
		temp = '''
	            <th>Index</th>
	            '''
		s = s + temp
		for col in cols:
			temp = '''
	                <th>''' + col + '''</th>
	            '''
			s = s + temp

		end_foot = '''</tr>
	        </tfoot>'''
		s = s + end_foot

		s = s + '''
	        <tbody>'''
		for index, row in df.iterrows():
			s = s + '''
	            <tr>'''
			temp = '''
	                <td>''' + str(index) + '''</td>
	                '''
			s = s + temp
			for i in xrange(0, len(row)):
				line = str(row[i])
				if len(line) <= 50:
					temp = '''
	                    <td>''' + line + '''</td>
	                    '''
					s = s + temp
				else:
					temp = '''
	                    <td>''' + line[0:50] + '''...''' + '''</td>
	                    '''
					s = s + temp
			s = s + '''
	            </tr>'''
		end_table = '''
	        </tbody>
	    </table>'''
		s = s + end_table
		s = s.encode('ascii', 'xmlcharrefreplace')
		return s

	def tohtml_library_MIBiG(self, df, tableId):
		cols = df.columns
		s = '''<table id="''' + tableId + '''" class="display table table-striped table-bordered" width="100%">'''
		s = s + '''
	        <thead>
	            <tr>'''
		for col in cols:
			temp = '''
	                <th>''' + col + '''</th>
	            '''
			s = s + temp
		end_header = '''</tr>
	        </thead>'''
		s = s + end_header

		s = s + '''
	        <tfoot>
	            <tr>'''
		for col in cols:
			print col
			if col=="BGC Product":
				col="BGC"
				temp = '''
						<th>''' + col + '''</th>
					'''
				s = s + temp
				continue
			elif col=="MIBiG Accession":
				col="MIBiG ac."
				temp = '''
						<th>''' + col + '''</th>
					'''
				s = s + temp
				continue
			else:
				temp = '''
					<th>''' + col + '''</th>
			    	'''
				s = s + temp

		end_foot = '''</tr>
	        </tfoot>'''
		s = s + end_foot

		s = s + '''
	        <tbody>'''
		for index, row in df.iterrows():
			s = s + '''
	            <tr>'''
			for i in xrange(0, len(row)):
				line = str(row[i])
				if len(line) <= 50:
					temp = '''
	                    <td>''' + line + '''</td>
	                    '''
					s = s + temp
				else:
					temp = '''
	                    <td>''' + line[0:50] + '''...''' + '''</td>
	                    '''
					s = s + temp
			s = s + '''
	            </tr>'''
		end_table = '''
	        </tbody>
	    </table>'''
		s = s + end_table
		s = s.encode('ascii', 'xmlcharrefreplace')
		return s

	def tohtml_library_NCBI(self, df, tableId, df2):
		cols = df.columns
		s = '''<table id="''' + tableId + '''" class="display table table-striped table-bordered" width="100%">'''
		s = s + '''
	        <thead>
	            <tr>'''
		for col in cols:
			temp = '''
	                <th>''' + col + '''</th>
	            '''
			s = s + temp
		end_header = '''</tr>
	        </thead>'''
		s = s + end_header

		s = s + '''
	        <tfoot>
	            <tr>'''
		for col in cols:
			print col
			if col=="Genome_FTP_path":
				col="FTP"
				temp = '''
						<th>''' + col + '''</th>
					'''
				s = s + temp
				continue
			elif col=="Assembly_accession":
				col="Assembly ac."
				temp = '''
						<th>''' + col + '''</th>
					'''
				s = s + temp
				continue
			else:
				temp = '''
					<th>''' + col + '''</th>
			    	'''
				s = s + temp

		end_foot = '''</tr>
	        </tfoot>'''
		s = s + end_foot

		s = s + '''
	        <tbody>'''
		for index, row in df.iterrows():
			s = s + '''
	            <tr>'''
			for i in xrange(0, len(row)):
				line = str(row[i])
				if len(line) <= 50:
					temp = '''
	                    <td>''' + line + '''</td>
	                    '''
					s = s + temp
				elif line.startswith("ftp://"):
						#line = "https"+line[3:]
						df1=df2[df2['Genome_FTP_path']==line]
						gene_id = df1.iloc[0, 4]
						temp = '''
							<td><a target="_blank" href="'''+line+'''">''' + gene_id + '''</a></td>
							'''
						s = s + temp
				elif " " in line:
					temp = '''
							<td>''' + line + '''</td>
							'''
					s = s + temp
				else:
					temp = '''
						<td>''' + line[0:50] + '''...''' + '''</td>
				    	'''
					s = s + temp
			s = s + '''
	            </tr>'''
		end_table = '''
	        </tbody>
	    </table>'''
		s = s + end_table
		s = s.encode('ascii', 'xmlcharrefreplace')
		return s

	def tohtml_library_parser(self, df, tableId):
		cols = df.columns
		s = '''<table id="''' + tableId + '''" class="display table table-striped table-bordered" width="100%">'''
		s = s + '''
	        <thead>
	            <tr>'''
		for col in cols:
			temp = '''
	                <th>''' + col + '''</th>
	            '''
			s = s + temp
		end_header = '''</tr>
	        </thead>'''
		s = s + end_header

		s = s + '''
	        <tfoot>
	            <tr>'''
		for col in cols:
			temp = '''
				<th>''' + col + '''</th>
				'''
			s = s + temp

		end_foot = '''</tr>
	        </tfoot>'''
		s = s + end_foot

		s = s + '''
	        <tbody>'''
		for index, row in df.iterrows():
			s = s + '''
	            <tr>'''
			for i in xrange(0, len(row)):
				line = str(row[i])
				if len(line) <= 50:
					temp = '''
	                    <td>''' + line + '''</td>
	                    '''
					s = s + temp
				else:
					temp = '''
	                    <td>''' + line[0:50] + '''...''' + '''</td>
	                    '''
					s = s + temp
			s = s + '''
	            </tr>'''
		end_table = '''
	        </tbody>
	    </table>'''
		s = s + end_table
		s = s.encode('ascii', 'xmlcharrefreplace')
		return s
