from Bio import SeqIO
import re
import csv
from os import listdir
from os.path import isfile, join

class NrpParser:

    def __init__(self):
        self.path_in = None
        self.results = None
        self.num_filesNoATCA = None
        self.A_domains = {}
        self.A_domains2 = {}
        self.All_domains = {}
        self.sec_met=[]

    def sec_Extractor(self, sec_met_list):
        self.sec_met.append(sec_met_list)
        domains_list = []
        subDomains = ["NRPS/PKS Domain: Condensation_Starter", "NRPS/PKS Domain: AMP-binding", "NRPS/PKS Domain: PCP","NRPS/PKS Domain: ACP",
                      "NRPS/PKS Domain: Condensation_", "NRPS/PKS Domain: Epimerization", "NRPS/PKS Domain: Cglyc", "NRPS/PKS Domain: PKS_AT","NRPS/PKS Domain: Heterocyclization"]
        for domain in sec_met_list:
            if (any(sub_domain in domain for sub_domain in subDomains)):
                domains_list.append(domain)
        return domains_list


    def purfiyA_domain(self, A_domain):
        searchObj = re.search('.{4}\(SANDPUMA ensemble\)', A_domain, re.M | re.I)
        aa = re.search('\w{3}', searchObj.group(), re.M | re.I)
        res= aa.group()
        if res=="all": res="nrp"
        return res

    def purfiyT_domain(self, T_domain):
        searchObj = re.search('-\d+', T_domain, re.M | re.I)
        t_site = re.search('\d+', searchObj.group(), re.M | re.I)
        return t_site.group()

    def purfiyC_domain(self, C_domain):
        c_site = re.search('\d+', C_domain, re.M | re.I)
        return c_site.group()

    def purfiyE_domain(self, E_domain):
        searchObj = re.search('-\d+', E_domain, re.M | re.I)
        E_site = re.search('\d+', searchObj.group(), re.M | re.I)
        return E_site.group()

    def purfiy_Translation(self, Trans_domain):
        trans = re.search('[A-Z]+', Trans_domain, re.M | re.I)
        return trans.group()

    def purfiyPos(self, domain):
        aa1 = re.search('\d+', domain, re.M | re.I)
        aa2 = re.search('-\d+', domain, re.M | re.I)
        return aa1.group() + aa2.group()

    def purfiyPos2(self, domain):
        aa1 = re.search('\d+', domain, re.M | re.I)
        aa2 = re.search(':\d+', domain, re.M | re.I)
        return aa1.group() + aa2.group()

    def purification1(self, domains):
        chunk_list = []
        for domain in domains:
            if ("AMP-binding" in domain):
                aa = self.purfiyA_domain(domain)
                pos_a = self.purfiyPos(domain)
                tu = (aa,pos_a)
                chunk_list.append(tu)
                continue

            if ("PKS_AT" in domain):
                aa = self.purfiyA_domain(domain)
                pos_as = self.purfiyPos(domain)
                tu = (aa, pos_as)
                chunk_list.append(tu)
                continue

            if ("PCP" in domain):
                t_site = self.purfiyT_domain(domain)
                pos_t = self.purfiyPos(domain)
                tu = (t_site, pos_t)
                chunk_list.append(tu)
                continue

            if ("ACP" in domain):
                t_site = self.purfiyT_domain(domain)
                pos_ts = self.purfiyPos(domain)
                tu = (t_site, pos_ts)
                chunk_list.append(tu)
                continue

            if ("Condensation" in domain):
                c_site = self.purfiyC_domain(domain)
                pos_c = self.purfiyPos(domain)
                tu = (c_site, pos_c)
                chunk_list.append(tu)
                continue

            if ("Cglyc" in domain):
                c_site = self.purfiyC_domain(domain)
                pos_cg = self.purfiyPos(domain)
                tu = (c_site, pos_cg)
                chunk_list.append(tu)
                continue

            if ("Heterocyclization" in domain):
                c_site = self.purfiyC_domain(domain)
                pos_cg = self.purfiyPos(domain)
                tu = (c_site, pos_cg)
                chunk_list.append(tu)
                continue

            if ("Epimerization" in domain):
                e_site = self.purfiyE_domain(domain)
                pos_e = self.purfiyPos(domain)
                tu = (e_site, pos_e)
                chunk_list.append(tu)
                continue

        return chunk_list

    def purification1A(self, domains):
        chunk_list = []
        for domain in domains:
            if ("AMP-binding" in domain):
                aa = self.purfiyA_domain(domain)
                chunk_list.append(aa)
                continue

            if ("PKS_AT" in domain):
                aa = self.purfiyA_domain(domain)
                chunk_list.append(aa)
                continue
        return chunk_list

    def purification2(self, domains):
        chunk_list = []
        for domain in domains:
            if ("AMP-binding" in domain):
                chunk_list.append('A')
                continue

            if ("PKS_AT" in domain):
                chunk_list.append('As')
                continue

            if ("PCP" in domain):
                chunk_list.append('T')
                continue

            if ("ACP" in domain):
                chunk_list.append('T')
                continue

            if ("Cglyc" in domain):
                chunk_list.append('C')
                continue

            if ("Heterocyclization" in domain):
                chunk_list.append('C')
                continue

            if ("X" in domain):
                chunk_list.append('X')
                continue

            if ("Condensation" in domain):
                chunk_list.append('C')
                continue

            if ("Epimerization" in domain):
                chunk_list.append('E')
                continue

        return chunk_list

    def purification2A(self, domains):
        chunk_list = []
        for domain in domains:
            if ("AMP-binding" in domain):
                chunk_list.append('A')
                continue
            if ("PKS_AT" in domain):
                chunk_list.append('As')
                continue
        return chunk_list

    def flatten(self, myListList):
        flatList = []
        for mylist in myListList:
            for e in mylist:
                flatList.append(e)
        return flatList


    def ATCA_builder(self, NRPS_list1, NRPS_list2):
        ATCA_block = []
        ATCA = ""
        ATCA_block_list = []
        for i in xrange(0, len(NRPS_list2)):
            d = NRPS_list2[i]
            if (d == "A"):
                if (ATCA == ""):
                    ATCA += d
                    ATCA_block.append(NRPS_list1[i])
                    continue
                elif ((ATCA == "ATC") or (ATCA == "AEC") ):
                    ATCA += d
                    ATCA_block.append(NRPS_list1[i])
                    ATCA_block_list.append(ATCA_block)
                    ATCA_block=[]
                    ATCA = "A"
                    ATCA_block.append(NRPS_list1[i])
                    continue
                else:
                    ATCA_block=[]
                    ATCA = "A"
                    ATCA_block.append(NRPS_list1[i])
                    continue

            if (d == "T"):
                if (ATCA == "A"):
                    ATCA += d
                    ATCA_block.append(NRPS_list1[i])
                    continue
                else:
                    ATCA = ""
                    ATCA_block = []
                    continue

            if (d == "C"):
                if (ATCA == "AT"):
                    ATCA += d
                    ATCA_block.append(NRPS_list1[i])
                    continue
                elif (ATCA == "AE"):
                    ATCA += d
                    ATCA_block.append(NRPS_list1[i])
                    continue
                else:
                    ATCA = ""
                    ATCA_block = []
                    continue

            if (d == "E"):
                if (ATCA == "AT"):
                    ATCA = "AE"
                    ATCA_block[1]= NRPS_list1[i]
                    continue
                else:
                    ATCA = ""
                    ATCA_block = []
                    continue

        return ATCA_block_list

    def ATCA_builder2(self, NRPS_list1, NRPS_list2):
        ATCA_block = []
        ATCA_block_list=[]
        ATCA = ""

        for i in range(len(NRPS_list2)):
            d = NRPS_list2[i]
            if (d == "A"):
                if (ATCA == "" and (i <= (len(NRPS_list2)) - 4)):
                    ATCA += d
                    ATCA_block.append(NRPS_list1[i])
                    continue
                if (ATCA == "ATC"):
                    ATCA += d
                    ATCA_block.append(NRPS_list1[i])
                    ATCA_block_list.append(ATCA_block)
                    ATCA_block=[]
                    if(i <= (len(NRPS_list2)) - 4):
                        ATCA = "A"
                        ATCA_block.append(NRPS_list1[i])
                        continue
                elif (i <= (len(NRPS_list2) - 4)):
                        ATCA_block = []
                        ATCA = "A"
                        ATCA_block.append(NRPS_list1[i])
                        continue

            if (d == "T"):
                if (ATCA == "A" and (i <= (len(NRPS_list2)) - 3)):
                    ATCA += d
                    ATCA_block.append(NRPS_list1[i])
                    continue
                else:
                    ATCA = ""
                    ATCA_block=[]
                    continue

            if (d == "C"):
                if (ATCA == "AT" and (i <= (len(NRPS_list2)) - 2)):
                    ATCA += d
                    ATCA_block.append(NRPS_list1[i])
                    continue
                else:
                    ATCA = ""
                    ATCA_block=[]
                    continue
        return ATCA_block_list

    def readFiles(self, path):
        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        return onlyfiles

    def NRPS_extractor(self, path):
        cluster_NRPS_blocks = {}
        cluster_All_domains = {}
        cluster_A_domains = {}
        cluster_A_domains2 = {}
        counter = 0
        files = self.readFiles(path)
        for cluster in files:
            cluster = cluster.strip()
            genome = SeqIO.read(path + cluster, "genbank")
        #n = len(filenames_list)
        #print "len_filenames:", n
        #for i in xrange(0, n):
            #filename = filenames_list[i]
            #f = preloader.dbx.files_download(path='/uploadedFiles/' + filename)
            #data = f[1].content
            #cluster = StringIO(data)
            #cluster = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            print "cluster:", cluster
            #genome = SeqIO.read(cluster, "genbank")  # you MUST tell SeqIO what format is being read
            print "Genome_ID:", genome.id
            NRPS_blocks = []
            NRPS_list3 = []

            ##### For A domains ####
            NRPS_list1A = []
            NRPS_list2A = []

            for feat in genome.features:
                NRPS_block = []
                ATCA_blocks = []
                qual = feat.qualifiers
                if (qual.has_key('sec_met')):
                    sec = qual.get('sec_met')
                    domains = self.sec_Extractor(sec)
                    if domains:
                        trans = qual.get('translation')
                        trans = trans[0]
                        NRPS_list1 = []
                        NRPS_list2 = []

                        NRPS_list1.append(self.purification1(domains))
                        NRPS_list2.append(self.purification2(domains))
                        #print NRPS_list1
                        #print NRPS_list2

                        NRPS_list3.append(self.purification2(domains))
                        NRPS_list3.append("-")

                        flattened1 = self.flatten(NRPS_list1)
                        flattened2 = self.flatten(NRPS_list2)
                        flattened3 = self.flatten(NRPS_list3)
                        #print flattened1
                        #print flattened2

                        ######## for A domains ########

                        NRPS_list1A.append(self.purification1A(domains))
                        NRPS_list2A.append(self.purification2A(domains))
                        NRPS_list1A.append("-")
                        NRPS_list2A.append("-")

                        flattened1A = self.flatten(NRPS_list1A)
                        flattened2A = self.flatten(NRPS_list2A)

                        ATCA_blocks.append(self.ATCA_builder(flattened1, flattened2))
                        ATCA_blocks = self.flatten(ATCA_blocks)
                        NRPS_block.append(ATCA_blocks)
                        NRPS_block.append(trans)
                        NRPS_block.append(feat.location)
                        NRPS_block.append(feat.strand)
                        NRPS_blocks.append(NRPS_block)
            counter += 1
            description = genome.description
            NRPS_blocks.append(description)
            NRPS_blocks.append(genome.id)
            cluster_NRPS_blocks[cluster] = NRPS_blocks
            cluster_A_domains[cluster] = flattened1A
            cluster_A_domains2[cluster] = flattened2A
            cluster_All_domains[cluster] = flattened3
        self.A_domains = cluster_A_domains
        self.A_domains2 = cluster_A_domains2
        self.All_domains = cluster_All_domains

        print "FINAL:", "Number of processed clusters:", counter, " Number of NRPS_blocks:", len(cluster_NRPS_blocks)
        return cluster_NRPS_blocks

    def linker_extractor(self, cluster_NRPS_blocks):
        cluster_linkers_dic = {}
        all_linkers = []
        for f in cluster_NRPS_blocks.keys():
            files_linkers = []
            NRPS_blocks = cluster_NRPS_blocks.get(f)
            for NRPS_block in NRPS_blocks:
                if not NRPS_block[0]:
                    continue
                if(isinstance(NRPS_block, str)):
                   continue
                trans = NRPS_block[1]
                gene_loc = NRPS_block[2]
                gene_loc = str(gene_loc.start.position) + "-" + str(gene_loc.end.position)
                print gene_loc
                gene_strand = NRPS_block[3]
                ATCA_blocks = NRPS_block[0]
                for ATCA_block in ATCA_blocks:
                    aa1 = ATCA_block[0][0]
                    aa1_pos = ATCA_block[0][1]

                    t_index = ATCA_block[1][0]
                    t_pos = ATCA_block[1][1]

                    c_index = ATCA_block[2][0]
                    c_pos = ATCA_block[2][1]

                    aa2 = ATCA_block[3][0]
                    aa2_pos = ATCA_block[3][1]

                    linker = trans[int(t_index):(int(c_index)-1)]
                    new_ATCA = [aa1, linker, aa2, NRPS_blocks[-2], NRPS_blocks[-1], len(linker), f, aa1_pos, aa2_pos, t_pos, c_pos, trans, gene_loc,  gene_strand]
                    all_linkers.append(new_ATCA)
                    files_linkers.append(new_ATCA)
            cluster_linkers_dic[f] = files_linkers
        return cluster_linkers_dic

    def linker_extractor2(self,genome):
        cluster_linkers_dic = {}
        all_linkers = []
        for f in genome.keys():
            files_linkers = []
            NRPS_blocks = genome.get(f)
            for NRPS_block in NRPS_blocks:
                if not NRPS_block[0]:
                    continue
                if (isinstance(NRPS_block, str)):
                    continue
                trans = NRPS_block[1]
                gene_loc = NRPS_block[2]
                gene_loc = str(gene_loc.start.position) + "-" + str(gene_loc.end.position)
                print gene_loc
                gene_strand = NRPS_block[3]
                ATCA_blocks = NRPS_block[0]
                for ATCA_block in ATCA_blocks:
                    aa1 = ATCA_block[0][0]
                    aa1_pos = ATCA_block[0][1]

                    t_index = ATCA_block[1][0]
                    t_pos = ATCA_block[1][1]

                    c_index = ATCA_block[2][0]
                    c_pos = ATCA_block[2][1]

                    aa2 = ATCA_block[3][0]
                    aa2_pos = ATCA_block[3][1]

                    linker = trans[int(t_index):(int(c_index) - 1)]
                    new_ATCA = [aa1, linker, aa2, NRPS_blocks[-2], NRPS_blocks[-1], len(linker), f, aa1_pos, aa2_pos,t_pos, c_pos, trans, gene_loc, gene_strand]
                    all_linkers.append(new_ATCA)
                    files_linkers.append(new_ATCA)
            cluster_linkers_dic[f] = files_linkers
        return all_linkers

    def filesNoNRPS(self, file_linkers_dic):
        counter = 0
        files_no_linkers = []
        for f in file_linkers_dic.keys():
            if not file_linkers_dic.get(f):
                counter += 1
                files_no_linkers.append(f)
        result = [counter, files_no_linkers]
        return result

    def ATCA_extractor(self, path):
        cluster_NRPS_blocks = self.NRPS_extractor(path)
        print
        print "#######################"
        print "Linker Extracttor"
        print "#######################"
        print
        file_linkers_dic = self.linker_extractor(cluster_NRPS_blocks)
        all_linkers = self.linker_extractor2(cluster_NRPS_blocks)
        self.results = [file_linkers_dic, all_linkers]

    def linkers_writer(self,results,path_out):
        all_linkers = results[1]
        writer = csv.writer(path_out)
        writer.writerows(all_linkers)


    def startExecution(self, path):
        self.ATCA_extractor(path)

    def writeFinalRes(self,results,path_out):
       temp=self.filesNoNRPS(results[0])
       self.num_filesNoATCA = temp[0]
       self.linkers_writer(results,path_out)
       print "Total number of files processed:", len(results[0]),"\n", "Total number of linkers retrieved:", len(results[1]),"\n", "Total number of files with no linkers: ",self.num_filesNoATCA
