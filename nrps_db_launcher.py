import pre_loader
from flask import Flask, render_template, request, make_response, send_file
import StringIO
import os
import pandas as pd
import ast
############### App ##################3
app = Flask(__name__)
app.secret_key = 'some_secret'
###### Pre-loading files #####
path = app.root_path
preloader = pre_loader.PreLoader(path)
listOfEdits = []
new_gene_list = []
#app.config['UPLOAD_FOLDER'] = os.path.join(path, "uploaded/")
#app.config['EXTRACTION_FOLDER'] = os.path.join(path, "extracted/")
#app.config['CLUSTER_FOLDER'] = os.path.join(path, "new_clusters/")
#app.config['DROPBOX_FOLDER_U'] = os.path.join(path, "dropbox_folder_u/")
#app.config['DROPBOX_FOLDER_D'] = '/NRPS_clusters_NCBI_66K_2/'
#app.config['PARAMS_FOLDER'] = os.path.join(path, "params/")

app.config['UPLOAD_FOLDER'] = os.path.join("/opt/app-root/src/uploads/", "uploaded/")
app.config['EXTRACTION_FOLDER'] = os.path.join("/opt/app-root/src/uploads/", "extracted/")
app.config['CLUSTER_FOLDER']= os.path.join("/opt/app-root/src/uploads/", "new_clusters/")
app.config['DROPBOX_FOLDER_U']= os.path.join("/opt/app-root/src/uploads/", "dropbox_folder_u/")
app.config['DROPBOX_FOLDER_D'] = '/NRPS_clusters_NCBI_66K_2/'
app.config['PARAMS_FOLDER'] = os.path.join("/opt/app-root/src/uploads/", "params/")

app.config['ALLOWED_EXTENSIONS'] = set(['gbk'])
if not os.path.exists(app.config['UPLOAD_FOLDER']):os.makedirs(app.config['UPLOAD_FOLDER'])
if not os.path.exists(app.config['EXTRACTION_FOLDER']):os.makedirs(app.config['EXTRACTION_FOLDER'])
if not os.path.exists(app.config['CLUSTER_FOLDER']):os.makedirs(app.config['CLUSTER_FOLDER'])
if not os.path.exists(app.config['DROPBOX_FOLDER_U']):os.makedirs(app.config['DROPBOX_FOLDER_U'])
if not os.path.exists(app.config['PARAMS_FOLDER']):os.makedirs(app.config['PARAMS_FOLDER'])

############################################################# HOME ###########################################################################
@app.route('/')
def main():
    return render_template('Home/home.html')

############################################################  PARSER ##########################################################################

@app.route('/Parser')
def doParse():
    preloader.nrps_parser.deleteUploadedFiles(app.config['UPLOAD_FOLDER'])
    preloader.nrps_parser.deleteUploadedFiles(app.config['EXTRACTION_FOLDER'])
    return render_template('Parser/parser.html')

@app.route('/Parser/Upload', methods=['POST'])
def upload():
    uploaded_files = request.files.getlist("file[]")
    filenames_dic, bad_filenames_dic = preloader.nrps_parser.uploadFiles(uploaded_files=uploaded_files, app=app)
    return render_template('Parser/upload.html', filenames_dic=filenames_dic, bad_filenames_dic=bad_filenames_dic)

######## ATCA Extractor ################

@app.route('/Parser/Upload/NewLinkers')
def getNewLinker():
    results = preloader.nrps_parser.extract_atca_linkers(path=app.config['UPLOAD_FOLDER'])
    extracted_linkers = results[0]
    extracted_linkers.to_csv(os.path.join(app.config['EXTRACTION_FOLDER'], ".csv"), index=False)
    extracted_linkers = extracted_linkers[['A1', 'Linker', 'A2', 'Length', 'Description', 'Cluster']]
    extracted_linkers = preloader.tohtml_library_parser(extracted_linkers, "extracted_linkers")
    return render_template('Parser/newTable.html', table=extracted_linkers,
                           title='Newly extracted linkers', num_filesNoATCA=results[1], num_files=results[2], num_linkers=results[3])

@app.route('/Parser/Upload/NewLinkers/Download')
def download_tab():
    extracted_linkers = pd.read_csv(os.path.join(app.config['EXTRACTION_FOLDER'], ".csv"))
    db = extracted_linkers
    s = StringIO.StringIO()
    db.to_csv(s, index=False)
    csv = s.getvalue()
    response = make_response(csv)
    cd = 'attachment; filename=Extracted_linkers.csv'
    response.headers['Content-Disposition'] = cd
    response.mimetype ='text/csv'
    return response

############################################################# Libraries ##############################################################################

@app.route('/Libraries/MIBiG_db')
def getMIBiGLinkers():
    MIBiG_db = preloader.MIBiG_DB1
    MIBiG_db_short = MIBiG_db[['A1', 'Linker', 'A2', 'Length', 'MIBiG_ID',"Product", 'Organism']]
    MIBiG_db = preloader.tohtml_library_MIBiG(MIBiG_db_short, "mibig")
    return render_template('libraries/MIBiG.html', table=MIBiG_db, title='MIBiG Database IMLs')

@app.route('/Libraries/NCBI_db')
def getNCBILinkers():
    NCBI_db = preloader.NCBI_HTML
    with open(NCBI_db, 'r') as myfile_html:
        data_html = myfile_html.read().replace('\n', '')
    return render_template('libraries/NCBI.html', table=data_html, title="NCBI RefSeq Database IMLs")


################################################### DeNovo Design #####################################################################

@app.route('/DenovoDesign/peptide_based')
def setPeptidebasedDesignPage():
    return preloader.nrps_design.setNourinebasedDesignPage(preloader=preloader)

@app.route('/DenovoDesign/peptide_based/newPeptide', methods=['POST'])
def getPeptideInfo():
    npeptide = request.form['Npeptide']
    return preloader.nrps_design.getNourineBasedPeptide(npeptide=npeptide, preloader=preloader, path=path)

@app.route('/DenovoDesign/peptide_based/newPeptide/result', methods=['POST'])
def getNewPeptide():
    return preloader.nrps_design.getNourineBasedPeptideResult(preloader=preloader, path=path, params=app.config['PARAMS_FOLDER'])

############################################################# DOWNLOAD ###########################################################################
@app.route('/Downloads/')
def downloadLib():
    return render_template('Download/download.html')

@app.route('/Downloads/download_MIBiG1/')
def download_MIBiG1():
    db = preloader.MIBiG_DB1
    s = StringIO.StringIO()
    db.to_csv(s,index=False)
    csv = s.getvalue()
    response = make_response(csv)
    cd = 'attachment; filename=MIBiG_linkers1.csv'
    response.headers['Content-Disposition'] = cd
    response.mimetype ='text/csv'
    return response

@app.route('/Downloads/download_MIBiG2/')
def download_MIBiG2():
    db = preloader.MIBiG_DB2
    s = StringIO.StringIO()
    db.to_csv(s,index=False)
    csv = s.getvalue()
    response = make_response(csv)
    cd = 'attachment; filename=MIBiG_linkers2.csv'
    response.headers['Content-Disposition'] = cd
    response.mimetype ='text/csv'
    return response

@app.route('/Downloads/download_NCBI1/')
def download_NCBI1():
    db = preloader.NCBI_DB1
    s = StringIO.StringIO()
    db.to_csv(s, index=False)
    csv = s.getvalue()
    response = make_response(csv)
    cd = 'attachment; filename=NCBI_linkers1.csv'
    response.headers['Content-Disposition'] = cd
    response.mimetype='text/csv'
    return response

@app.route('/Downloads/download_NCBI2/')
def download_NCBI2():
    db = preloader.NCBI_DB2
    s = StringIO.StringIO()
    db.to_csv(s, index=False)
    csv = s.getvalue()
    response = make_response(csv)
    cd = 'attachment; filename=NCBI_linkers2.csv'
    response.headers['Content-Disposition'] = cd
    response.mimetype='text/csv'
    return response

@app.route('/Downloads/download_NCBI3/')
def download_NCBI3():
    db = preloader.NCBI_DB3
    s = StringIO.StringIO()
    db.to_csv(s, index=False)
    csv = s.getvalue()
    response = make_response(csv)
    cd = 'attachment; filename=NCBI_RefSeq_genomes_list.csv'
    response.headers['Content-Disposition'] = cd
    response.mimetype='text/csv'
    return response

@app.route('/Downloads/download_zip/')
def download_zip():
    filename = "clusters_gbks.zip"
    path_zip = os.path.join(path,"files2Read/BGC_clusters.zip")
    return send_file(filename_or_fp=path_zip, attachment_filename=filename, as_attachment=True, mimetype='application/zip')


############################################################# DATA_SUMMARY ###########################################################################

@app.route('/Data_Summary/')
def getDataSummary():
    return render_template('Data_Summary/data_summary.html')

############################################################# Analysis ###########################################################################

@app.route('/Analysis/')
def getAnalysis():
    return render_template('Analysis/analysis2.html')


################### Nav-Bar ###########
@app.route('/Contact')
def getinfo():
    return render_template("Navbar/personal.html")

@app.route('/Uploads_Examples')
def getuploadExamples():
    return render_template("Navbar/upload_example.html")

@app.route('/Uploads_Examples/Daptomycin_BGC')
def getDaptomycin():
    return send_file(os.path.join(path,'files2Read/templates/daptomycin.gbk'),  attachment_filename="Daptomycin_BGC.gbk", as_attachment=True)


@app.route('/Uploads_Examples/Vancomycin_BGC')
def getVancomycin():
    return send_file(os.path.join(path,'files2Read/templates/vancomycin.gbk'), attachment_filename="Vancomycin_BGC.gbk", as_attachment=True)


@app.route('/Tutorial')
def getTutorial():
    return render_template("Navbar/tutorial2.html")

@app.route('/postmethod_3', methods = ['POST', 'GET'])
def get_post_javascript_data3():
    global listOfEdits
    listOfEdits = request.form['javascript_data']
    listOfEdits = ast.literal_eval(listOfEdits)
    return "Done"

@app.route('/postmethod', methods = ['POST', 'GET'])
def get_post_javascript_data():
    global listOfEdits
    global new_gene_list

    postValues = request.form['np2']
    print "length", len(postValues.split("/"))
    print postValues
    peptideName = postValues.split("/")[0]


    changes = postValues.split("/")[1]

    the_filename = postValues.split("/")[2]
    print "1", the_filename
    dropbox_u = app.config['DROPBOX_FOLDER_U']
    dropbox_d = app.config['DROPBOX_FOLDER_D']
    changes = ast.literal_eval(changes)
    #the_filename = ast.literal_eval(the_filename)
    #print "2",  the_filename

    rows_pass, super_params =preloader.nrps_design.check_list_of_edits(listOfEdits, app.config['PARAMS_FOLDER'], the_filename)
    if not rows_pass:
        error = "Please make sure to select the correct number of linkers" + " [" + str(len(super_params[13])) + " linkers] " +\
                "and the right pairs of linkers, which are as follows:"
        return render_template('Design/novoPeptide11.html', graph1=super_params[0], id=super_params[1],
                               orginalPep=super_params[2], graph2=super_params[3], table_old=super_params[4],
                               table_new=super_params[5], title_new=super_params[6],
                               title_old=super_params[7], old_seq=super_params[8],
                               new_seq=super_params[9], pairs=super_params[10], no_linkers=super_params[11],
                               mod=super_params[12], changes=super_params[13], the_filename=the_filename, error=error)
    else:
        print "@@@@@@@@@@@@@@@@@@@      Rows chosen by the User     @@@@@@@@@@@@@@@@@@@"
        print "Peptidname according to client request:", peptideName
        print "changes according to client request:", changes
        print "______________________________________________________________________"
        bigTable = preloader.NCBI_DB4

        toEdit_cluster_info = pd.read_csv(path + "/original_linkers/" + peptideName + ".csv")
        toEdit_cluster = os.path.join(path, "files2Read/templates/"+peptideName+".gbk")

        new_gene_list = preloader.nrps_design.creat2bEdited(listOfEdits, toEdit_cluster, toEdit_cluster_info,changes, bigTable, preloader, dropbox_u, dropbox_d)
        return preloader.nrps_design.settingNewCluster(preloader, new_gene_list, app)


@app.route('/newGene/fasta_download', methods = ['POST', 'GET'])
def download_new_gene_fasta2():
    id = request.form['np2']
    print app.config['CLUSTER_FOLDER']
    filePath = os.path.join(app.config['CLUSTER_FOLDER'], id)
    filename = "new_cluster_"+id+".fasta"
    return send_file(filename_or_fp=filePath, attachment_filename=filename, as_attachment=True, mimetype='application/fasta')

if __name__ == '__main__':
    app.run(threaded=True, host='0.0.0.0')

