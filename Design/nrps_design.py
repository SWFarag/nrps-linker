import xml.etree.ElementTree as ET
import networkx as nx
from matplotlib import pyplot as plt
from flask import render_template, request
import pandas as pd
import os
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random, string
import pickle

class NRP_Design:

    def __init__(self):
        self.name = None

    def setNourinebasedDesignPage(self, preloader):
        name_list = preloader.peptides_ids_xmlFileNames.ix[:, 1]
        peptide_a = sorted(name_list)
        return render_template('Design/peptideBased.html', peptide_a=peptide_a)

    def getNourineBasedPeptide(self, npeptide, preloader, path):
        template_info = preloader.templates_info
        nrp_info = template_info[template_info['TemplateName'] == npeptide]
        nrp_comp = nrp_info.iloc[0, 1]
        nrps_counter = nrp_info.iloc[0, 2]

        npeptideInfo = self.getPeptideInfo_xml(npeptide=npeptide, preloader=preloader, path=path)
        record2 = preloader.peptides_ids_graphs[preloader.peptides_ids_graphs['Peptide'] == npeptide.upper()]
        g = record2.iloc[0, 3]
        c = g.split("@")
        c1 = c[0]
        nodes = c1.split(",")
        graph = "/static/nourineGraphs/" + npeptide + ".png"
        size = len(nodes)
        link = "http://bioinfo.lifl.fr/norine/result.jsp?ID=" + npeptideInfo[3]
        return render_template('Design/novoPeptide1.html', id=npeptideInfo[3], npeptide=npeptide, category=npeptideInfo[0],
                               activities=npeptideInfo[1], organism=npeptideInfo[2], composition=nodes, graph=graph,
                               size=size, link=link, nrps_counter=nrps_counter, nrp_comp=nrp_comp)

    def getNourineBasedPeptideResult(self, preloader, path, params):
        originalPep = request.form['np2']
        nrp_comp = request.form['np3']
        nrps_counter = request.form['np4']
        seq = request.form['peptideSeq']
        print
        print "%%%%%%%%%%%%%%%%%%%%%%%%%   Chosen Template Info  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print 'Template: ', originalPep
        print "Original sequence: ", nrp_comp
        print "Number of NRPSs: ", nrps_counter
        print "Modified sequence:", seq
        print
        record = preloader.peptides_ids_graphs[preloader.peptides_ids_graphs['Peptide'] == originalPep.upper()]
        npeptideInfo = self.getPeptideInfo_xml(originalPep, preloader=preloader, path=path)

        results = self.getEdge_nodes(record=record, seq=seq)
        link = "http://bioinfo.lifl.fr/norine/result.jsp?ID=" + npeptideInfo[3]
        graph1 = path + "static/nourineGraphs/" + originalPep + ".png"
        nrp_comp_list = nrp_comp.split("--")
        seq_list = seq.split("--")
        if len(nrp_comp_list) != len(seq_list):
            error1 = "Number of NRPSs genes do not match. It should be exactly " + str(nrps_counter) + " but " + str(len(seq_list)) + " is given."
            return render_template('Design/novoPeptide1.html', id=npeptideInfo[3], npeptide=originalPep,
                                   category=npeptideInfo[0], activities=npeptideInfo[1], organism=npeptideInfo[2],
                                   composition=results[0], graph=graph1, size=len(results[0]),
                                   link=link, error1=error1, nrps_counter=nrps_counter, nrp_comp=nrp_comp)

        elif not self.checklengthpergene(nrp_comp_list, seq_list):
            error2 = "Number of monomers within NRPSs genes do not match."
            return render_template('Design/novoPeptide1.html', id=npeptideInfo[3], npeptide=originalPep,
                                   category=npeptideInfo[0], activities=npeptideInfo[1], organism=npeptideInfo[2],
                                   composition=results[0], graph=graph1, size=len(results[0]),
                                   link=link, error2=error2, nrps_counter=nrps_counter, nrp_comp=nrp_comp)
        else:
            print "######## Compute Differences #############"
            diffs, num_mod = self.computeDiff(nrp_comp_list, seq_list)
            nodiff = self.checkdiffs(diffs)
            print
            if not nodiff:
                isGood = self.isSingleGood(nrp_comp_list, seq_list)
                if not isGood:
                    error4 = "A change within a single monomer containing gene has been detected." \
                           " Please modify only NRPSs genes with more than one monomer."
                    return render_template('Design/novoPeptide1.html', id=npeptideInfo[3], npeptide=originalPep,
                                           category=npeptideInfo[0], activities=npeptideInfo[1],
                                           organism=npeptideInfo[2],
                                           composition=results[0], graph=graph1, size=len(results[0]),
                                           link=link, error4=error4, nrps_counter=nrps_counter, nrp_comp=nrp_comp)
                print "####### Retrieve required Linkers ############"
                print
                t_new, t_origin = self.computeNewLinkers(diffs, nrp_comp_list, seq_list, preloader, originalPep)
                if not os.path.exists(os.path.join(path, "original_linkers")): os.makedirs(os.path.join(path, "original_linkers"))
                t_origin[0].to_csv(path+"/original_linkers/"+originalPep+".csv", index=False)
                ####### new data ##########
                new_linkers = t_new
                imls_new = new_linkers[0]
                if imls_new is None:
                    pairs = new_linkers[1]
                    no_linkers = new_linkers[2]
                    return render_template('Design/novoPeptide1.html', id=npeptideInfo[3], npeptide=originalPep,
                                           category=npeptideInfo[0], activities=npeptideInfo[1], organism=npeptideInfo[2],
                                           composition=results[0], graph=graph1, size=len(results[0]), link=link,
                                           nrps_counter=nrps_counter, nrp_comp=nrp_comp,
                                           pairs=pairs, no_linkers=no_linkers)
                else:
                    imls1_new = imls_new[['A1', 'Linker', 'A2', 'Length', 'Description', 'Cluster']]
                    data_new = preloader.tohtml_design(imls1_new, "table_new")

                ####### old_data ##########
                old_linkers = t_origin
                imls_old = old_linkers[0]
                imls1_old = imls_old[['A1', 'Linker', 'A2', 'Antibiotic', 'MIBiD_ID', 'Length', 'Cluster']]
                if imls1_old is None:
                    data_old = imls1_old
                else:
                    data_old = preloader.tohtml_design(imls1_old, "table_old")

                pairs = new_linkers[1]
                no_linkers = new_linkers[2]
                if no_linkers:
                    return render_template('Design/novoPeptide1.html', id=npeptideInfo[3], npeptide=originalPep,
                                           category=npeptideInfo[0], activities=npeptideInfo[1],
                                           organism=npeptideInfo[2],
                                           composition=results[0], graph=graph1, size=len(results[0]), link=link,
                                           nrps_counter=nrps_counter, nrp_comp=nrp_comp,
                                           pairs=pairs, no_linkers=no_linkers)
                else:
                    self.buildPeptideGraph(originalPep=originalPep, new_seq=seq_list, results=results, path=path)

                    graph1 = "/static/nourineGraphs/" + originalPep + '.png'
                    graph2 = "/static/designedGraphs/" + originalPep + '.png'
                    super_params_x = [graph1, id, originalPep, graph2, data_old,
                                    data_new, 'New Linkers Table:', "Old Linkers Table:", nrp_comp,
                                    seq, pairs, no_linkers, num_mod, t_origin[3]]
                    temp_name="super_params_"+str(super_params_x[-4])+".txt"
                    the_filename = params+temp_name
                    print "the the_filename", the_filename
                    with open(the_filename, 'wb') as f:
                        pickle.dump(super_params_x, f)
                    return render_template('Design/novoPeptide11.html', graph1=graph1, id=id, orginalPep=originalPep,
                                           graph2=graph2, table_old=data_old, table_new=data_new,title_new='New Linkers Table:', title_old='Old Linkers Table:',
                                           old_seq=nrp_comp, new_seq=seq, pairs=pairs, no_linkers=no_linkers,
                                           mod=num_mod, changes=t_origin[3], numberChanges=len(t_origin[3]), the_filename=temp_name)
            else:
                error3 = "No peptide modifications observed, please substitute at least a single residue!"
                return render_template('Design/novoPeptide1.html', id=npeptideInfo[3], npeptide=originalPep,
                                       category=npeptideInfo[0], activities=npeptideInfo[1], organism=npeptideInfo[2],
                                       composition=results[0], graph=graph1, size=len(results[0]),
                                       link=link, error3=error3, nrps_counter=nrps_counter, nrp_comp=nrp_comp)

############# Auxilary methods ####################

    def computeDiff(self, seq1, seq2):
        diffs = []
        for j in xrange(0, len(seq1)):
            comp1 = seq1[j].split(",")
            comp2 = seq2[j].split(",")
            diff = [((i - 1), i, (i + 1), comp2[i].lower()) for i in xrange(0, len(comp2)) if comp2[i].lower() not in comp1[i].lower()]
            diffs.append(diff)
            print ("Differences", len(diff), diff)
        number_of_modification=0
        for diff in diffs:
            number_of_modification += len(diff)
        return diffs, number_of_modification

    def getEdge_nodes(self, record, seq):
        id = record.iloc[0, 0]
        g = record.iloc[0, 3]
        c = g.split("@")
        edges = c[1:len(c)]
        c1 = c[0]
        old_seq = c1.split(",")
        #print "Original: ", old_seq
        new_seq = seq.split(",")
        #print "New Seq: ", new_seq
        results = [old_seq, new_seq, id, g, edges]
        return results

    def computeNewLinkers(self, diffs, old_seq, new_seq, preloader, originalPep):
        n = 0
        changes_total = []
        pairs_new = set()
        pairs_origin = set()
        for i in xrange(0, len(diffs)):
            diff = diffs[i]
            changes_local = []
            if diff:
                for j in xrange(0, len(diff)):
                    t = diff[j]
                    old_seq1 = old_seq[i].split(",")
                    #print "old:", old_seq1, len(old_seq1)
                    new_seq1 = new_seq[i].split(",")
                    #print "new:", new_seq1
                    ##### Right most terminal monomer ###
                    if t[1] == (len(old_seq1) - 1):
                        pair_new = (new_seq1[t[0]].lower(), new_seq1[t[1]].lower())
                        pairs_new.add(pair_new)

                        pair_origin = (old_seq1[t[0]].lower(), old_seq1[t[1]].lower())
                        pairs_origin.add(pair_origin)

                        if pair_new[0] == pair_origin[0]:
                            pair_new = (new_seq1[t[0]].lower(), new_seq1[t[1]].lower(), 2)
                        elif pair_new[1] == pair_origin[1]:
                            pair_new = (new_seq1[t[0]].lower(), new_seq1[t[1]].lower(), 1)
                        else:
                            pair_new = (new_seq1[t[0]].lower(), new_seq1[t[1]].lower(), 3)

                        changes_local.append((pair_origin, pair_new, n, i))
                        n +=1
                    ##### left most terminal monomer #####
                    elif t[1] == 0:
                        pair_new = (new_seq1[t[1]].lower(), new_seq1[t[2]].lower())
                        pairs_new.add(pair_new)

                        pair_origin = (old_seq1[t[1]].lower(), old_seq1[t[2]].lower())
                        pairs_origin.add(pair_origin)

                        #### modified monomer pos , either first or second pos)
                        if pair_new[0] == pair_origin[0]:
                            pair_new = (new_seq1[t[1]].lower(), new_seq1[t[2]].lower(), 2)
                        elif pair_new[1] == pair_origin[1]:
                            pair_new = (new_seq1[t[1]].lower(), new_seq1[t[2]].lower(), 1)
                        else: pair_new = (new_seq1[t[1]].lower(), new_seq1[t[2]].lower(), 3)

                        changes_local.append((pair_origin, pair_new, n, i))
                        n += 1

                    ##### Middle monomer ###
                    else:
                        pair1_new = (new_seq1[t[0]].lower(), new_seq1[t[1]].lower())
                        pair2_new = (new_seq1[t[1]].lower(), new_seq1[t[2]].lower())
                        pairs_new.add(pair1_new)
                        pairs_new.add(pair2_new)

                        pair1_origin = (old_seq1[t[0]].lower(), old_seq1[t[1]].lower())
                        pair2_origin = (old_seq1[t[1]].lower(), old_seq1[t[2]].lower())
                        pairs_origin.add(pair1_origin)
                        pairs_origin.add(pair2_origin)

                        if pair1_new[0] == pair1_origin[0]:
                            pair1_new = (new_seq1[t[0]].lower(), new_seq1[t[1]].lower(), 2)
                        elif pair1_new[1] == pair1_origin[1]:
                            pair1_new = (new_seq1[t[0]].lower(), new_seq1[t[1]].lower(), 1)
                        else:
                            pair1_new = (new_seq1[t[0]].lower(), new_seq1[t[1]].lower(), 3)

                        if pair2_new[0] == pair2_origin[0]:
                            pair2_new = (new_seq1[t[1]].lower(), new_seq1[t[2]].lower(), 2)
                        elif pair2_new[1] == pair2_origin[1]:
                            pair2_new = (new_seq1[t[1]].lower(), new_seq1[t[2]].lower(), 1)
                        else:
                            pair2_new = (new_seq1[t[1]].lower(), new_seq1[t[2]].lower(), 3)

                        changes_local.append((pair1_origin, pair1_new, n, i))
                        n += 1
                        changes_local.append((pair2_origin, pair2_new, n, i))
                        n += 1

            changes_l = self.unifychanges(changes_local)
            for c in changes_l:
                changes_total.append(c)
        print
        print "new pair: ", pairs_new
        print "original pair: ", pairs_origin
        print "changes: ", changes_total
        print
        print "####### Original Linkers ########"
        t_origin = self.buildTablesOrigin(pairs_origin, preloader, originalPep, changes_total)
        print
        print "####### New Linkers ########"
        t_new = self.buildTablesNew(pairs_new, preloader, changes_total)
        print
        return (t_new, t_origin)


    def buildTablesNew(self, pairs, preloader, changes):
        imls_db = preloader.NCBI_DB4
        frames = []
        emptyData = []
        for p in pairs:
            if '-' in p[0]:
                a1 = p[0].split('-')[1]
            else: a1 = p[0]

            if '-' in p[1]:
                a2 = p[1].split('-')[1]
            else: a2 = p[1]
            data1 = imls_db[(imls_db["A1"] == a1) & (imls_db["A2"] == a2)]
            dim = data1.shape
            print "dim", dim
            if dim[0] == 0:
                emptyData.append((a1, a2))
            else:
                frames.append(data1)
        if frames:
            result_new = pd.concat(frames)
            return (result_new, pairs, emptyData, changes)
        else:
            result = None
            return (result, pairs, emptyData, changes)

    def buildTablesOrigin(self, pairs, preloader, originalPep, changes):
        originalPep2 = originalPep.split(" ")[0]
        imls_db = preloader.MIBiG_DB5
        frames = []
        emptyData = []
        locs = []
        if originalPep2=="Ile-polymyxin":
            data1 = imls_db[imls_db["Antibiotic"] == originalPep2 + " biosynthetic gene cluster"]
            data1 = data1["Gene_location"]
            set_d = {}
            counter = 0
            for d in set(data1):
                set_d[counter] = d
                counter += 1
            for p in changes:
                if '-' in p[0][0]:
                    a1 = p[0][0].split('-')[1]
                    locs.append(p[3])
                else:
                    a1 = p[0][0]
                    locs.append(p[3])

                if '-' in p[0][1]:
                    a2 = p[0][1].split('-')[1]
                    locs.append(p[3])
                else:
                    a2 = p[0][1]
                    locs.append(p[3])

                gene_loc = set_d.get(p[3])
                data1 = imls_db[(imls_db["A1"] == a1) & (imls_db["A2"] == a2) & (imls_db["Antibiotic"] == originalPep2+" biosynthetic gene cluster") &
                                (imls_db["Gene_location"] == gene_loc )]

                dim = data1.shape
                print "dim", dim
                if dim[0] == 0:
                    emptyData.append((a1, a2))
                else:
                    frames.append(data1)
        else:
            for p in changes:
                if '-' in p[0][0]:
                    a1 = p[0][0].split('-')[1]
                    locs.append(p[3])
                else:
                    a1 = p[0][0]
                    locs.append(p[3])

                if '-' in p[0][1]:
                    a2 = p[0][1].split('-')[1]
                    locs.append(p[3])
                else:
                    a2 = p[0][1]
                    locs.append(p[3])

                data1 = imls_db[(imls_db["A1"] == a1) & (imls_db["A2"] == a2) & (imls_db["Antibiotic"] == originalPep2+" biosynthetic gene cluster")]

                dim = data1.shape
                print "dim", dim
                if dim[0] == 0:
                    emptyData.append((a1, a2))
                else:
                    frames.append(data1)

        if frames:
            result_new = pd.concat(frames)
            result_new = result_new.drop_duplicates()
            return (result_new, pairs, emptyData, changes)
        else:
            result = None
            return (result, pairs, emptyData, changes)


    def buildPeptideGraph(self, originalPep, new_seq, results, path):
        old_seq = results[0]
        edges = results[4]

        seq_edi = []
        for s in new_seq:
            for i in s.split(","):
                seq_edi.append(i)

        new_seq = seq_edi
        if len(new_seq) != len(old_seq):
            monomer = [old_seq[0]]
            new_seq = monomer+new_seq

        G = nx.Graph(title=originalPep)
        edge_list = set()
        for j in xrange(0, len(edges)):
            outdegrees = edges[j]
            e = outdegrees.split(",")
            for n in e:
                t1 = (j, int(n))
                t2 = (int(n), j)
                if (t1 in edge_list) or (t2 in edge_list):
                    continue
                else:
                    G.add_edge(j, int(n))
                    edge_list.add(t1)
                    edge_list.add(t2)

        for n in xrange(0, len(old_seq)):
            G.add_node(n, label=new_seq[n])

        values = range(0, len(old_seq))
        mapL = {}
        for v in values:
            mapL[v] = new_seq[v]
        nx.draw_shell(G, cmap=plt.get_cmap('prism'), node_color=values, node_size=2500, with_labels=True, labels=mapL)
        plt.suptitle("Modified Peptide")
        saveImgPath = os.path.join(path, os.path.join("static/designedGraphs", originalPep + '.png'))
        plt.savefig(saveImgPath)  # save as png
        plt.clf()

    def getPeptideInfo_xml(self, npeptide, preloader, path):
        df_nourine = preloader.peptides_ids_xmlFileNames
        record1 = df_nourine[df_nourine['Antibiotic'] == npeptide]
        xml_file = record1.iloc[0, 2]
        id = record1.iloc[0, 0]
        tree = ET.parse(os.path.join(path,"files2Read/Nourine_info/Nourine_Xml/" + xml_file + ".xml"))
        root = tree.getroot()

        category = ''
        activities = []
        organism = ''

        for peptide in root.findall('.//category'):
            category = peptide.text

        for peptide in root.findall('.//activities/activity'):
            activities.append(peptide.text)

        for peptide in root.findall('.//organism/name'):
            organism = peptide.text

        link = "http://bioinfo.lifl.fr/norine/result.jsp?ID=" + id

        return (category,activities,organism,id, link)

    def checkdiffs(self, diffs):
        nodiff = True
        for diff in diffs:
            if not diff:
                continue
            else: nodiff = False
        return nodiff

    def isSingleGood(self, old_seq, new_seq):
        isGood = False
        ones = [(old_seq[i], i)for i in xrange(0, len(old_seq)) if len(old_seq[i].split(",")) == 1]
        if ones:
            for one in ones:
                if one[0] == new_seq[one[1]]:
                    isGood = True
                else:
                    isGood = False
                    break
        else:
            isGood = True
        return isGood

    def checklengthpergene(self, ori_list, new_list):
        subLength = False
        for i in xrange(0, len(ori_list)):
            sub_l1 = ori_list[i].split(",")
            sub_l2 = new_list[i].split(",")
            if len(sub_l1) == len(sub_l2):
                subLength = True
            else:
                subLength = False
                break
        return subLength


    ################################### Design New Cluster ################################

    def check_compatability1(self, record, start, end, translation):
        print "Determining Gene Directionality:"
        f = record.seq
        r = record.seq.reverse_complement()
        print "start:", start, "end:", end
        gene = f[start:end]
        trans = translation[0:4]

        try:
            print gene.translate(table="Bacterial", cds=True)[0:4] == trans
            print "Gene is ff"
            return gene, "f"
        except Exception as e:
            s = str(e)
            if (s.startswith("First codon")) and (s.endswith("is not a start codon")):
                gene_r = gene.reverse_complement()
                try:
                    print gene_r.translate(table="Bacterial", cds=True)[0:4] == trans
                    print "fr..good stop condon"
                    return gene_r, "rc"
                except Exception as e:
                    s = str(e)
                    if (s.startswith("Final codon")) and (s.endswith("stop codon")):
                        print gene_r.translate(table="Bacterial")[0:4] == trans
                        print "fr..bad stop condon"
                        return gene_r, "rc"
                    else:
                        print gene_r.translate(table="Bacterial")[0:4] == trans
                        print "fr..", s
                        return gene_r, "rc"
            else:
                print gene.translate(table="Bacterial")[0:4] == trans
                print "ff..bad stop condon"
                return gene, "f"

    def getLinker(self, gene, directionality1, end_t, start_c, directionality2):
        linker = gene[end_t * 3:(start_c) * 3]
        print "original", linker
        print("LINKER: ", len(linker))
        if directionality1 == 'f' and directionality2 == 'f':
            return linker
        elif directionality1 == 'f' and directionality2 == 'rc':
            linker = linker.reverse_complement()
            return linker
        else:
            return linker

    def getADomain(self, gene, directionality1, start_a, end_a, directionality2):
        a_domain = gene[start_a * 3:end_a * 3]
        print "original_a_domain", a_domain
        print("A_domain: ", len(a_domain))
        if directionality1 == 'f' and directionality2 == 'f':
            return a_domain
        elif directionality1 == 'f' and directionality2 == 'rc':
            a_domain = a_domain.reverse_complement()
            return a_domain
        else:
            return a_domain

    def dnaTranslator(self, dna):
        return dna.translate(table="Bacterial")

    def setNewGenesCreation(self, action_dic, toEdit_cluster, preloader, dropbox_u):
        new_gene_dic = {}
        for key in action_dic.iterkeys():
            print
            print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    Gene_location: " + str(key) + "    $$$$$$$$$$$$$$$$$$$$$$$$$$$$"
            print
            actions = action_dic.get(key)
            new_gene = self.createNewGene(actions, toEdit_cluster, preloader, dropbox_u)
            new_gene_dic[key] = new_gene
        return new_gene_dic

    def createNewGene(self, actions, toEdit_cluster, preloader, dropbox_u):
        if len(actions) > 1:
            print "************************************"
            print "Multiple modifications per NRPS gene:", str(len(actions)) + " modifications"
            print "************************************"
        else:
            print "*********************************"
            print "Single modification per NRPS gene:", str(len(actions)) + " modification(s)"
            print "*********************************"
        print
        print "........................................"
        print "Order A_domains by location ascendingly: "
        print "........................................"
        print
        act_loc_dic = {}
        actions_sorted_in_order = []
        for act in actions:
            act_ori = act[0]
            act_ori_pos = act_ori[-1]
            if act_ori_pos != 3:
                loc_a_domain = (act[0][5], act[0][6])
                if act_loc_dic.has_key(loc_a_domain):
                    act_loc_dic.get(loc_a_domain).append(act)
                else:
                    act_loc_dic[loc_a_domain] = [act]
            else:
                    loc_a_domain = (act[0][5], act[0][6], act[0][7], act[0][8])
                    if act_loc_dic.has_key(loc_a_domain):
                        act_loc_dic.get(loc_a_domain).append(act)
                    else:
                        act_loc_dic[loc_a_domain] = [act]

        sorted_a_domain_locs = act_loc_dic.keys()
        sorted_a_domain_locs.sort()
        print "Sorted_a_domains within that Gene:", sorted_a_domain_locs
        print
        for loc in sorted_a_domain_locs:
            for ac in act_loc_dic.get(loc):
                actions_sorted_in_order.append(ac)

        linkers = []
        a_domains = []
        genes_dic = []

        end_t_domains = []
        start_c_domains = []
        start_a_domains = []
        end_a_domains = []

        a_domains_ori_t = []
        a_domains_mod_t = []

        for i in xrange(0, len(actions_sorted_in_order)):
            action = actions_sorted_in_order[i]
            print "@@@@@@@@@@@@@@@@@@@@@@@@@"
            print "Action_"+str(i+1) + " in Gene_location: "+str(action[0][0]) + "_" +str(action[0][1])
            print "@@@@@@@@@@@@@@@@@@@@@@@@@"
            row_ori = action[0]
            row_mod = action[1]
            cluster1 = toEdit_cluster
            cluster2 = action[2]

            pos_ori = row_ori[-1]
            if pos_ori != 3:
                gene_start1 = row_ori[0]
                gene_end1 = row_ori[1]
                translation1 = row_ori[2]
                end_t_domain1 = row_ori[3]
                start_c_domain1 = row_ori[4]
                start_a_domain1 = row_ori[5]
                end_a_domain1 = row_ori[6]

                end_t_domains.append((end_t_domain1, "end_t_domain"))
                start_c_domains.append((start_c_domain1, "start_c_domain"))
                start_a_domains.append((start_a_domain1, "start_a_domain"))
                end_a_domains.append((end_a_domain1, "end_a_domain"))

            else:
                gene_start1 = row_ori[0]
                gene_end1 = row_ori[1]
                translation1 = row_ori[2]
                end_t_domain1 = row_ori[3]
                start_c_domain1 = row_ori[4]
                start_a_domain11 = row_ori[5]
                end_a_domain11 = row_ori[6]
                start_a_domain12 = row_ori[7]
                end_a_domain12 = row_ori[8]

                end_t_domains.append((end_t_domain1, "end_t_domain"))
                start_c_domains.append((start_c_domain1, "start_c_domain"))
                start_a_domains.append((start_a_domain11, "start_a_domain"))
                end_a_domains.append((end_a_domain11, "end_a_domain"))
                start_a_domains.append((start_a_domain12, "start_a_domain"))
                end_a_domains.append((end_a_domain12, "end_a_domain"))
            a1_ori = row_ori[-3]
            a2_ori = row_ori[-2]



            if pos_ori == 1:
                a_domains_ori_t.append((a1_ori, start_a_domain1, end_a_domain1))
            elif pos_ori == 2:
                a_domains_ori_t.append((a2_ori, start_a_domain1, end_a_domain1))
            else:
                a_domains_ori_t.append((a1_ori, start_a_domain11, end_a_domain11))
                a_domains_ori_t.append((a2_ori, start_a_domain12, end_a_domain12))


            print
            print "^^^^^^^^^^^^^^^^^^^^^^^"
            print "Cluster to be edited"
            print "^^^^^^^^^^^^^^^^^^^^^^^"
            print gene_start1
            print gene_end1
            print translation1[0:5]
            print end_t_domain1
            print start_c_domain1
            if pos_ori != 3:
                print start_a_domain1
                print end_a_domain1
            else:
                print start_a_domain11
                print end_a_domain11
                print
                print start_a_domain12
                print end_a_domain12

            pos_mod = row_mod[-1]
            if pos_mod != 3:
                gene_start2 = row_mod[0]
                gene_end2 = row_mod[1]
                translation2 = row_mod[2]
                end_t_domain2 = row_mod[3]
                start_c_domain2 = row_mod[4]
                start_a_domain2 = row_mod[5]
                end_a_domain2 = row_mod[6]
            else:
                gene_start2 = row_mod[0]
                gene_end2 = row_mod[1]
                translation2 = row_mod[2]
                end_t_domain2 = row_mod[3]
                start_c_domain2 = row_mod[4]
                start_a_domain21 = row_mod[5]
                end_a_domain21 = row_mod[6]
                start_a_domain22 = row_mod[7]
                end_a_domain22 = row_mod[8]

            a1_mod = row_mod[-3]
            a2_mod = row_mod[-2]


            if pos_mod == 1:
                a_domains_mod_t.append((a1_mod, start_a_domain2, end_a_domain2))
            elif pos_mod == 2:
                a_domains_mod_t.append((a2_mod, start_a_domain2, end_a_domain2))
            else:
                a_domains_mod_t.append((a1_mod, start_a_domain21, end_a_domain21))
                a_domains_mod_t.append((a2_mod, start_a_domain22, end_a_domain22))

            print
            print "^^^^^^^^^^^^^^^^^^^^^^^"
            print "Cluster needed for edit"
            print "^^^^^^^^^^^^^^^^^^^^^^^"
            print gene_start2
            print gene_end2
            print translation2[0:5]
            print end_t_domain2
            print start_c_domain2
            if pos_mod != 3:
                print start_a_domain2
                print end_a_domain2
            else:
                print start_a_domain21
                print end_a_domain21
                print
                print start_a_domain22
                print end_a_domain22
            print

            cluster_name = cluster2.split("/")[-1]
            print cluster_name
            print dropbox_u + cluster_name
            print cluster2
            print cluster1
            preloader.dbx.files_download_to_file(dropbox_u+cluster_name, cluster2)
            cluster2 = dropbox_u+cluster_name
            record1 = SeqIO.read(cluster1, "genbank")
            record2 = SeqIO.read(cluster2, "genbank")

            print "Retrieve Gene1:"
            print gene_start1, gene_end1
            gene1, directionality1 = self.check_compatability1(record1, gene_start1, gene_end1, translation1)
            genes_dic.append((gene1, directionality1))
            print
            print "Retrieve Gene2:"
            gene2, directionality2 = self.check_compatability1(record2, gene_start2, gene_end2, translation2)
            print
            print "###### linker #####"
            linker2 = self.getLinker(gene2, directionality1, end_t_domain2, start_c_domain2, directionality2)
            linkers.append((end_t_domain1, start_c_domain1, linker2))
            print self.dnaTranslator(linker2)
            print"Linker_2b_added_dna: ", linker2
            print
            print "###### A_domain #####"
            if pos_ori != 3:
                a_domain2 = self.getADomain(gene2, directionality1, start_a_domain2, end_a_domain2, directionality2)
                a_domains.append((start_a_domain1, end_a_domain1, a_domain2))
                print self.dnaTranslator(a_domain2)
                print "A_domain_AA length:", len(a_domain2) / 3, " A_domain_dna length:", len(a_domain2)
            else:
                a_domain21 = self.getADomain(gene2, directionality1, start_a_domain21, end_a_domain21, directionality2)
                a_domains.append((start_a_domain11, end_a_domain11, a_domain21))
                print self.dnaTranslator(a_domain21)
                print "A_domain_AA length:", len(a_domain21) / 3, " A_domain_dna length:", len(a_domain21)
                print
                a_domain22 = self.getADomain(gene2, directionality1, start_a_domain22, end_a_domain22, directionality2)
                a_domains.append((start_a_domain12, end_a_domain12, a_domain22))
                print self.dnaTranslator(a_domain22)
                print "A_domain_AA length:", len(a_domain22) / 3, " A_domain_dna length:", len(a_domain22)
            print

        final_a_domains = self.clean_a_domains(a_domains_ori_t, a_domains_mod_t, a_domains)
        print "----------------------------------------------------------------------"
        print "Clean and Unify needed A_domains based on the unique Orginal A_domains: Done"
        print "----------------------------------------------------------------------"
        print

        print
        print "###################"
        print "All_a_domains length", len(a_domains)
        print "Final_a_domains", len(final_a_domains)
        all_domains = self.makeSorteList([end_t_domains, start_c_domains, start_a_domains, end_a_domains])
        genes_dic = list(set(genes_dic))
        gene = genes_dic[0][0]
        dic = genes_dic[0][1]
        new_gene = self.get_right_case3(linkers, gene, final_a_domains, all_domains)

        if dic == "rc":
            new_gene = new_gene.reverse_complement()
            print new_gene[0:3]
        new_gene = self.checkforexception(new_gene)
        return new_gene

    def creat2bEdited(self, listOfEdits, toEdit_cluster, toEdit_cluster_info, changes, bigTable, preloader, dropbox_u, dropbox_d):

        genes_changes = self.set_cluster2b_edited(toEdit_cluster_info, changes)
        print
        print "######## Number of NRPS genes to be modified ###########"
        print
        print "Number of NRPS genes to be modified:", len(genes_changes)
        print
        nrps_action = self.find_right_edit_set1(genes_changes, changes, listOfEdits, bigTable, dropbox_d)
        print "**************************"
        #print "nrps_action", nrps_action
        print "genes_changes values", len(genes_changes.values())
        #print "list of edits:", listOfEdits
        print "changes:", changes
        print "**************************"
        print
        new_gene_list = []

        print "######### RE-engineer NRPS genes ##############"
        print
        action_dic = {}
        for act in nrps_action:
            loc1 = (act[0][0],act[0][1])
            if action_dic.has_key(loc1):
                action_dic.get(loc1).append(act)
            else:
                action_dic[loc1] = [act]

        new_gene_dic = self.setNewGenesCreation(action_dic, toEdit_cluster, preloader, dropbox_u)

        print
        sorted_locs = new_gene_dic.keys()
        print "Unsorted Gene loctions:",  sorted_locs
        sorted_locs.sort()
        print "Sorted Gene locations:", sorted_locs
        for loc in sorted_locs:
            new_gene_list.append((loc,new_gene_dic.get(loc)))

        new_gene_list.append(toEdit_cluster)
        print "new_gene_list_1:", len(new_gene_list)-1
        return new_gene_list


    def checkforexception(self, new_gene):
        try:
            new_gene.translate(table=11, cds=True)
            if new_gene[0:3] != "ATG":
                new_gene = "ATG" + new_gene[3:]
        except Exception as e:
            s = str(e)
            print s
            new_gene.translate(table=11)
            if new_gene[0:3] != "ATG":
                new_gene = "ATG" + new_gene[3:]
        return new_gene

    def find_right_edit_set1(self, genes_changes, changes, listOfEdits, bigTable, dropbox_d):
        nrps_action = []
        for nrps_gene in genes_changes.iterkeys():
            for r in genes_changes.get(nrps_gene):
                row_ori = r
                a1_ori = row_ori[-3]
                a2_ori = row_ori[-2]
                for key in changes:
                    a1_temp = key[0][0][-3:]
                    a2_temp = key[0][1][-3:]
                    if a1_ori == a1_temp and a2_ori == a2_temp:
                        a1_mod = key[1][0][-3:]
                        a2_mod = key[1][1][-3:]
                        pos_mod = key[1][2]
                        row_mod, cluster2 = self.find_right_edit_set2(a1_mod=a1_mod, a2_mod=a2_mod, listOfEdits=listOfEdits, pos_mod=pos_mod, bigTable=bigTable, dropbox_d=dropbox_d)
                        nrps_action.append([row_ori, row_mod, cluster2])
        return nrps_action

    def find_right_edit_set2(self, a1_mod, a2_mod, listOfEdits, pos_mod, bigTable, dropbox_d):
        for edit in listOfEdits:
            if edit[1] == a1_mod and edit[3] == a2_mod:
                index = int(edit[0])
                row = bigTable.iloc[[index], :]
                cluster_name = row.iloc[0, 6]
                cluster_edit = os.path.join(dropbox_d, cluster_name)
                gene_start2 = int(str(row.iloc[0, -2]).split("-")[0])
                gene_end2 = int(str(row.iloc[0, -2]).split("-")[1])
                translation2 = str(row.iloc[0, -3])
                end_t_domain2 = int(str(row.iloc[0, -6]).split("-")[1])
                start_c_domain2 = int(str(row.iloc[0, -5]).split("-")[0])

                if pos_mod == 1:
                    start_a_domain2 = int(str(row.iloc[0, -7]).split("-")[0])
                    end_a_domain2 = int(str(row.iloc[0, -7]).split("-")[1])
                elif pos_mod == 2:
                    start_a_domain2 = int(str(row.iloc[0, -4]).split("-")[0])
                    end_a_domain2 = int(str(row.iloc[0, -4]).split("-")[1])
                else:
                    start_a_domain21 = int(str(row.iloc[0, -7]).split("-")[0])
                    end_a_domain21 = int(str(row.iloc[0, -7]).split("-")[1])
                    start_a_domain22 = int(str(row.iloc[0, -4]).split("-")[0])
                    end_a_domain22 = int(str(row.iloc[0, -4]).split("-")[1])
                if pos_mod != 3:
                    strand2 = int(str(row.iloc[0, -1]))
                    row_mod = [gene_start2, gene_end2, translation2, end_t_domain2, start_c_domain2, start_a_domain2,
                                end_a_domain2, strand2, a1_mod, a2_mod, pos_mod]
                else:
                    strand2 = int(str(row.iloc[0, -1]))
                    row_mod = [gene_start2, gene_end2, translation2, end_t_domain2, start_c_domain2, start_a_domain21,
                               end_a_domain21, start_a_domain22, end_a_domain22, strand2, a1_mod, a2_mod, pos_mod]
                return row_mod, cluster_edit

    def set_cluster2b_edited(self, toEdit_cluster_info, changes):
        genes_changes = {}
        for index, row in toEdit_cluster_info.iterrows():
            a1 = str(row[0])
            a2 = str(row[2])
            gene_start1 = int(str(row[12]).split("-")[0])
            gene_end1 = int(str(row[12]).split("-")[1])
            translation1 = str(row[11])
            end_t_domain1 = int(str(row[9]).split("-")[1])
            start_c_domain1 = int(str(row[10]).split("-")[0])
            pos = self.find_the_right_a_domain(a1, a2, changes, index)
            if pos == 1:
                start_a_domain1 = int(str(row[7]).split("-")[0])
                end_a_domain1 = int(str(row[7]).split("-")[1])
            elif pos == 2:
                start_a_domain1 = int(str(row[8]).split("-")[0])
                end_a_domain1 = int(str(row[8]).split("-")[1])
            else:
                start_a_domain11 = int(str(row[7]).split("-")[0])
                end_a_domain11 = int(str(row[7]).split("-")[1])
                start_a_domain12 = int(str(row[8]).split("-")[0])
                end_a_domain12 = int(str(row[8]).split("-")[1])
            strand1 = int(row[13])
            if pos != 3:
                row_mods = [gene_start1, gene_end1, translation1, end_t_domain1, start_c_domain1, start_a_domain1,
                            end_a_domain1, strand1, a1, a2, pos]
                if genes_changes.has_key(row[-2]):
                    l = genes_changes.get(row[-2])
                    l.append(row_mods)
                else:
                    genes_changes[row[-2]] = [row_mods]
            else:
                row_mods = [gene_start1, gene_end1, translation1, end_t_domain1, start_c_domain1, start_a_domain11,
                            end_a_domain11, start_a_domain12, end_a_domain12, strand1, a1, a2, pos]
                if genes_changes.has_key(row[-2]):
                    l = genes_changes.get(row[-2])
                    l.append(row_mods)
                else:
                    genes_changes[row[-2]] = [row_mods]
        return genes_changes

    def find_the_right_a_domain(self, a1, a2, changes, index):
        for key in changes:
            if key[0][0][-3:] == a1 and key[0][1][-3:] == a2 and key[2] >= index:
                pos = key[1][2]
                return pos

    def clean_a_domains(self, a_domains_ori_t, a_domains_mod_t, a_domains):
        needed_indices=[]
        temp_ori = []
        for i in xrange(0, len(a_domains_ori_t)):
            ad = a_domains_ori_t[i]
            if ad in temp_ori: continue
            else:
                temp_ori.append(ad)
                needed_indices.append(i)
        final_a_domains_mod_t = [a_domains_mod_t[i] for i in needed_indices]
        print "curcial", final_a_domains_mod_t
        final_a_domains_seq = [a_domains[i] for i in needed_indices]
        return final_a_domains_seq

    def domainExist(self, final_a_domains, ad):
        domainExist = False
        for i in final_a_domains:
            if i[0] == ad[0]:
                domainExist = True
                break
        return domainExist

    def get_right_case3(self, linkers, gene, a_domains, all_domains):   
        linkers_copy = linkers[:]
        a_domains_copy = a_domains[:]

        everything = []
        for e in linkers:
            t = (e[0],"z",e[2])
            everything.append(t)
        for e in a_domains:
            t=(e[0],"z",e[2])
            everything.append(t)
        for e in all_domains:
            everything.append(e)
        print everything
        everything.sort()
        print everything

        print "Linkers:",  linkers
        print "domains:", a_domains

        dom_text = []
        s = ''
        print
        print "######################"
        print "Re-engineering Starts", s
        print "######################"
        first_element = everything[0]
        s += gene[0:first_element[0] * 3]
        print
        print "Adding first element", len(s)
        every_edit = everything[1:-1]
        for i in xrange(0, len(every_edit)):
            ele = every_edit[i]
            if ele[1] == "start_a_domain":
                a = dom_text[-1]
                s += gene[a[0] * 3:ele[0] * 3]
                dom_text.append(ele)

            elif ele[1] == "end_a_domain":
                dom_text.append(ele)
                continue

            elif ele[1] == "end_t_domain":
                a = dom_text[-1]
                s += gene[a[0] * 3:ele[0] * 3]
                dom_text.append(ele)

            elif ele[1] == "start_c_domain":
                dom_text.append(ele)
                continue

            else:
                s += ele[2]
                #print "edit", len(ele[2]), "total:", len(s), "I just added a linker or a domain"
        print
        print "Adding last domain to the Gene"
        last_element = everything[-1]

        if last_element[1] == "end_a_domain":
            s += gene[(last_element[0]) * 3:]
        else:s += gene[last_element[0] * 3:]

        print
        print "================================="
        print "TEST_EDITS WITHIN A SINGLE GENE:"
        print "================================="
        print

        print "length_old_gene:", len(gene)
        print "length_new_gene:", len(s)
        print
        gene_diff = abs(len(gene) - len(s))
        print "difference in genes:", gene_diff
        sum_links_mod = 0
        sum_links_ori = 0
        for link in linkers_copy:
            sum_links_mod += len(link[2])
            print "end_t_domain:", link[0], ":", "start_c_domain:", link[1]
            sum_links_ori += (abs(int(link[0])-int(link[1]))) * 3

        print "sum_links_mod", sum_links_mod
        print "sum_links_ori", sum_links_ori
        print
        print
        sum_dom_mod = 0
        sum_dom_ori = 0
        for do in a_domains_copy:
            print "domain_length", len(do[2])
            sum_dom_mod += len(do[2])
            print "star_a_domain:", do[0], ":", "end_a_domain:", do[1]
            sum_dom_ori += (abs(int(do[0]) - int(do[1]))) * 3

        print "sum_dom_mod", sum_dom_mod
        print "sum_dom_ori", sum_dom_ori

        a = sum_dom_mod + sum_links_mod
        b = sum_dom_ori + sum_links_ori
        print "Sum of All domains and Linkers within mod", a
        print "Sum of All domains and Linkers within ori", b
        mods_diff = abs(a - b)
        print "difference in mods:",  mods_diff
        print
        print "###############################"
        print "Are both differences Conserved:", gene_diff == mods_diff
        print "###############################"
        print
        return s

    def makeSorteList(self, l):
        all_list = []
        for e1 in l:
            for e2 in e1:
                all_list.append(e2)
        all_list = list(set(all_list))
        all_list.sort()
        return all_list

    def setDomains(self, domains):
        print "number of domains even:", len(domains) % 2 == 0
        d_t = []
        for i in xrange(0, len(domains), 2):
            d = domains[i]
            d_1 = domains[i+1]
            d_t.append((d[0], d_1[0], d_1[1]))
        return d_t

    def createNewCluster(self, new_gene_list, record):
        locs = []
        genes = []
        new_cluster = ""
        cluster = record.seq
        print "length of original cluster", len(cluster)
        for ng in new_gene_list:
            new_gene = ng[1]
            gene_start = ng[0][0]
            gene_end = ng[0][1]
            print "Gene_start_end", gene_start, gene_end
            genes.append(new_gene)
            locs.append(gene_start)
            locs.append(gene_end)
        genes.reverse()
        locs_edited = locs[1:-1]
        mid_locs = []
        if len(locs_edited) != 0:
            for i in xrange(0, len(locs_edited),2):
                mid_locs.append((locs_edited[i], locs_edited[i+1]))

            if locs[0] == 1:
                new_cluster += genes.pop()
                for loc in mid_locs:
                    new_cluster +=cluster[loc[0]:loc[1]] + genes.pop()
                new_cluster += cluster[locs[-1]:]
            else:
                new_cluster += cluster[0:locs[0]] + genes.pop()
                for loc in mid_locs:
                    new_cluster += cluster[loc[0]:loc[1]] + genes.pop()
                new_cluster += cluster[locs[-1]:]
        else:
            if locs[0] == 1:
                new_cluster += genes.pop()
                new_cluster += cluster[locs[-1]:]
            else:
                new_cluster += cluster[0:locs[0]] + genes.pop()
                new_cluster += cluster[locs[-1]:]

        print "length of new cluster", len(new_cluster)
        # if new_cluster[0:3]!= "ATG":
        #     print "Changing to the Correct Start Codon"
        #     new_cluster = "ATG" + new_cluster[3:]
        return new_cluster

    def createFastas(self, app, new_gene_list, record):
        for ng in new_gene_list:
            new_gene = ng[1]
            gene_start = ng[0][0]
            gene_end = ng[0][1]
            record_seq = SeqRecord(new_gene, id=record.id, description=record.description + " " + str(gene_start) + "-" +str(gene_end))
            filePath = os.path.join(app.config['CLUSTER_FOLDER'], str(record.id) + "_" + str(gene_start) +"_" + str(gene_end))
            SeqIO.write(record_seq, filePath, "fasta")

    def word(self, length):
        return ''.join(random.choice(string.lowercase) for i in range(length))

    def check_list_of_edits(self, list_of_edits, params, tempName):
        the_filename = params + tempName
        print "edilo", the_filename
        with open(the_filename, 'rb') as f:
            super_params = pickle.load(f)
        changes = super_params[-1]
        changes_pure = [(change[1][0][-3:], change[1][1][-3:]) for change in changes]
        print "changes_pure", changes_pure
        n = len(list_of_edits)

        rows_pass = False
        if n == len(changes):
            for i in xrange(0,n):
                edit = list_of_edits[i]
                row_chosen = (edit[1], edit[3])
                #print row_chosen
                if row_chosen in changes_pure:
                    rows_pass = True
                    changes_pure.remove(row_chosen)
                else:
                    rows_pass = False
                    break
        else:
            rows_pass = False
        return rows_pass, super_params

    def randomword(self, length):
        return ''.join(random.choice(string.lowercase) for i in range(length))

    def settingNewCluster(self, preloader, new_gene_list, app):
        ng_list = [i[1] for i in new_gene_list[0:-1]]
        loc_list = [i[0] for i in new_gene_list[0:-1]]
        gene_dic = {}
        for i in xrange(0, len(ng_list)):
            loc = loc_list[i]
            gene_dic[loc] = ng_list[i]

        record = SeqIO.read(new_gene_list[-1], "genbank")

        ### New_gene ####
        print
        preloader.nrps_design.createFastas(app, new_gene_list[0:-1], record)
        print "###########################################"
        print "Creating Fasta file for every edited genes: DONE"
        print "###########################################"


        ### New_Cluster ####
        print
        print
        print "######### Creating New Cluster ############"
        new_cluster = preloader.nrps_design.createNewCluster(new_gene_list[0:-1], record)
        record_seq = SeqRecord(new_cluster, id=record.id, description=record.description)

        filePath = os.path.join(app.config['CLUSTER_FOLDER'], record.id)
        SeqIO.write(record_seq, filePath, "fasta")

        rand_words = [preloader.nrps_design.randomword(6) for i in xrange(0, len(new_gene_list))]
        indices_pos = [i for i in xrange(0, len(new_gene_list[0:-1]))]

        return render_template('Design/newCluster.html', new_gene=new_gene_list[0:-1], new_cluster=new_cluster,
                               record_id=record.id, rand_words=rand_words, indices_pos=indices_pos)

    def unifychanges(self, changes):
        changes_clean = []
        bad_indices = []
        for i in xrange(0,len(changes)):
            change = changes[i]
            pair_old = change[0]
            pair_new = change[1]
            t_edit = (pair_old, pair_new)
            if t_edit not in changes_clean:
                changes_clean.append(t_edit)
            else: bad_indices.append(i)

        for j in bad_indices:
            del changes[j]
        return changes


