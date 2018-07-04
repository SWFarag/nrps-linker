import pre_loader
from flask import Flask, render_template, request, make_response, send_file
import StringIO
import os
import pandas as pd
import ast
############### App ##################3
app = Flask(__name__)
###### Pre-loading files #####
path = app.root_path
preloader = pre_loader.PreLoader(path)
listOfEdits = []
new_gene_list = []
app.config['UPLOAD_FOLDER'] = os.path.join(path, "uploaded/")
app.config['EXTRACTION_FOLDER'] = os.path.join(path, "extracted/")
app.config['DROPBOX_FOLDER_D'] = '/NRPS_clusters_NCBI_66K_2/'

# app.config['UPLOAD_FOLDER'] = os.path.join("/opt/app-root/src/uploads/", "uploaded/")
# app.config['EXTRACTION_FOLDER'] = os.path.join("/opt/app-root/src/uploads/", "extracted/")
# app.config['CLUSTER_FOLDER']= os.path.join("/opt/app-root/src/uploads/", "new_clusters/")
# app.config['DROPBOX_FOLDER_U']= os.path.join("/opt/app-root/src/uploads/", "dropbox_folder_u/")
# app.config['DROPBOX_FOLDER_D'] = '/NRPS_clusters_NCBI_66K_2/'
# app.config['PARAMS_FOLDER'] = os.path.join("/opt/app-root/src/uploads/", "params/")

app.config['ALLOWED_EXTENSIONS'] = set(['gbk'])
if not os.path.exists(app.config['UPLOAD_FOLDER']):os.makedirs(app.config['UPLOAD_FOLDER'])
if not os.path.exists(app.config['EXTRACTION_FOLDER']):os.makedirs(app.config['EXTRACTION_FOLDER'])
############################################################# HOME ###########################################################################
@app.route('/')
def main():
    return render_template('Home/home.html')

############################################################  PARSER ##########################################################################

@app.route('/Parser3/<int:wrongParser>/')
def doParseThree(wrongParser):
    preloader.nrps_parser.deleteUploadedFiles(app.config['UPLOAD_FOLDER'])
    preloader.nrps_parser.deleteUploadedFiles(app.config['EXTRACTION_FOLDER'])
    return render_template('Parser/parser.html', wrongParser=wrongParser)

@app.route('/Parser4/<int:wrongParser>/')
def doParseFour(wrongParser):
    preloader.nrps_parser.deleteUploadedFiles(app.config['UPLOAD_FOLDER'])
    preloader.nrps_parser.deleteUploadedFiles(app.config['EXTRACTION_FOLDER'])
    return render_template('Parser/parserFour.html', wrongParser=wrongParser)

@app.route('/Parser/Upload3', methods=['POST'])
def uploadThree():
    uploaded_files = request.files.getlist("file[]")
    filenames_dic, bad_filenames_dic = preloader.nrps_parser.uploadFiles(uploaded_files=uploaded_files, app=app)
    return render_template('Parser/upload.html', filenames_dic=filenames_dic, bad_filenames_dic=bad_filenames_dic)

@app.route('/Parser/Upload4', methods=['POST'])
def uploadFour():
    uploaded_files = request.files.getlist("file[]")
    filenames_dic, bad_filenames_dic = preloader.nrps_parser.uploadFiles(uploaded_files=uploaded_files, app=app)
    return render_template('Parser/uploadFour.html', filenames_dic=filenames_dic, bad_filenames_dic=bad_filenames_dic)

######## ATCA Extractor ################

@app.route('/Parser/Upload/NewLinkersAntiSMASH_3')
def getNewLinkerThree():
    wrongParser=False
    try:
        results = preloader.nrps_parser.extract_atca_linkers_Three(path=app.config['UPLOAD_FOLDER'])
    except Exception as e:
        wrongParser=True
        return doParseThree(wrongParser)
    ##results = preloader.nrps_parser.extract_atca_linkers_Three(path=app.config['UPLOAD_FOLDER'])
    extracted_linkers = results[0]
    extracted_linkers.to_csv(os.path.join(app.config['EXTRACTION_FOLDER'], "extracted.csv"), index=False)
    extracted_linkers = extracted_linkers[['A1', 'Linker', 'A2', 'Length', 'Description', 'Cluster']]
    extracted_linkers = preloader.tohtml_library_parser(extracted_linkers, "extracted_linkers")
    return render_template('Parser/newTable.html', table=extracted_linkers,
                           title='Newly extracted linkers', num_filesNoATCA=results[1], num_files=results[2], num_linkers=results[3])

@app.route('/Parser/Upload/NewLinkersAntiSMASH_4')
def getNewLinkerFour():
    wrongParser=False
    try:
        results = preloader.nrps_parser.extract_atca_linkers_Four(path=app.config['UPLOAD_FOLDER'])
    except Exception as e:
        wrongParser=True
        return doParseFour(wrongParser)
    #results = preloader.nrps_parser.extract_atca_linkers_Four(path=app.config['UPLOAD_FOLDER'])
    extracted_linkers = results[0]
    extracted_linkers.to_csv(os.path.join(app.config['EXTRACTION_FOLDER'], "extracted.csv"), index=False)
    extracted_linkers = extracted_linkers[['A1', 'Linker', 'A2', 'Length', 'Description', 'Cluster']]
    extracted_linkers = preloader.tohtml_library_parser(extracted_linkers, "extracted_linkers")
    return render_template('Parser/newTable.html', table=extracted_linkers,
                           title='Newly extracted linkers', num_filesNoATCA=results[1], num_files=results[2], num_linkers=results[3])

@app.route('/Parser/Upload/NewLinkers/Download')
def download_tab():
    extracted_linkers = pd.read_csv(os.path.join(app.config['EXTRACTION_FOLDER'], "extracted.csv"))
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
    MIBiG_db = preloader.tohtml_library_MIBiG2(MIBiG_db_short, "mibig")
    return render_template('libraries/MIBiG.html', table=MIBiG_db, title='MIBiG IMLs')

@app.route('/Libraries/NCBI_db')
def getNCBILinkers():
    NCBI_db = preloader.NCBI_HTML
    with open(NCBI_db, 'r') as myfile_html:
        #data_html = myfile_html.read().replace('\n', '')
        data_html = myfile_html.read()
    return render_template('libraries/NCBI.html', table=data_html, title="Putative NRPS IMLs")

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
    return render_template('Analysis/analysis.html')


################### Nav-Bar ###########
@app.route('/Contact')
def getinfo():
    return render_template("Navbar/personal.html")

@app.route('/Uploads_Examples')
def getuploadExamples():
    return render_template("Navbar/upload_example.html")

@app.route('/Uploads_Examples/AntiSMASH3.0_files/')
def getAntiSMASH3_files():
    filename = "antiSMASH3.0_BGC_cluster_files.zip"
    path_zip = os.path.join(path,"files2Read/templates/clusters_files_3.zip")
    return send_file(filename_or_fp=path_zip, attachment_filename=filename, as_attachment=True, mimetype='application/zip')

@app.route('/Uploads_Examples/AntiSMASH4.0_files/')
def getAntiSMASH4_files():
    filename = "antiSMASH4.0_BGC_cluster_files_4.zip"
    path_zip = os.path.join(path,"files2Read/templates/clusters_files_4.zip")
    return send_file(filename_or_fp=path_zip, attachment_filename=filename, as_attachment=True, mimetype='application/zip')


@app.route('/Tutorial')
def getTutorial():
    return render_template("Navbar/tutorial.html")

if __name__ == '__main__':
    app.run(threaded=True, host='0.0.0.0')
