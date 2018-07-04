from werkzeug.utils import secure_filename
import atca_parser_latest
import atca_parser_latest_anti4
import pandas as pd
import os

class NRPS_Parser:

    def __init__(self):
        self.atca_parser = atca_parser_latest.NrpParser()
        self.atca_parser_Four = atca_parser_latest_anti4.NrpParser()

    # For a given file, return whether it's an allowed type or not
    def allowed_file(self, filename, app):
        return '.' in filename and \
               filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']


    def uploadFiles(self, uploaded_files, app):
        # Get the name of the uploaded files
        filenames_dic = {}
        bad_filenames_dic = {}
        for file in uploaded_files:
            if file:
               # Check if the file is one of the allowed types/extensions
               if self.allowed_file(file.filename, app):
                   # Make the filename safe, remove unsupported chars
                   filename = secure_filename(file.filename)
                   file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                   # Save the filename into a list, we'll use it later
                   filenames_dic[filename]=file
                   # Redirect the user to the uploaded_file route, which
                   # will basicaly show on the browser the uploaded file
                   # Load an html page with a link to each uploaded file
               else:
                   filename = secure_filename(file.filename)
                   bad_filenames_dic[filename] = "wrong-extension"
            else:filenames_dic
        return filenames_dic, bad_filenames_dic

    def extract_atca_linkers_Three(self, path):
        self.atca_parser.startExecution(path)
        r = self.atca_parser.results
        cols = ["A1", "Linker", "A2", "Description", "Accession", "Length", "Cluster", "A1_pos", "A2_pos", "t_pos", "c_pos", "trans", "gene_loc", "gene_strand"]
        extracted_data = pd.DataFrame(r[1], columns=cols)
        temp = self.atca_parser.filesNoNRPS(r[0])
        num_filesNoATCA = temp[0]
        num_files = len(r[0])
        num_linkers = len(r[1])
        results = [extracted_data, num_filesNoATCA, num_files, num_linkers]
        self.deleteUploadedFiles(path)
        return results

    def extract_atca_linkers_Four(self, path):
        self.atca_parser_Four.startExecution(path)
        r = self.atca_parser_Four.results
        cols = ["A1", "Linker", "A2", "Description", "Accession", "Length", "Cluster", "A1_pos", "A2_pos", "t_pos", "c_pos", "trans", "gene_loc", "gene_strand"]
        extracted_data = pd.DataFrame(r[1], columns=cols)
        temp = self.atca_parser_Four.filesNoNRPS(r[0])
        num_filesNoATCA = temp[0]
        num_files = len(r[0])
        num_linkers = len(r[1])
        results = [extracted_data, num_filesNoATCA, num_files, num_linkers]
        self.deleteUploadedFiles(path)
        return results

    def deleteUploadedFiles(self, path):
        folder = path
        for the_file in os.listdir(folder):
            file_path = os.path.join(folder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)
