#!/Users/Stephen/Library/Enthought/Canopy_64bit/User/bin/python
import os
import subprocess

import MySQLdb as mdb
from flask import current_app

from get_effective_pathways import PathwaySummaryBasic, \
    GenericPathwayFileProcessor, PathwayDetailedFileWriter


class PathwaySummaryParsed(PathwaySummaryBasic):
    """Holds pathway information, pulled from anaysis output."""
    def __init__(self, pathway_number):
        PathwaySummaryBasic.__init__(self, pathway_number)
        # above adds path_id, n_actual, n_effective, p_value
        self.name = None
        self.root_url = "http://www.broadinstitute.org/gsea/msigdb/cards/"
        self.url = self.fetch_url()
        self.gene_set = set()

    def as_string(self):
        """Return string for javascript."""
        name_str = self.name.replace('_', ' ')
        name_str = name_str.replace('BIOCARTA ', '')
        name_str = name_str.replace('REACTOME ', '')
        name_str = name_str.replace('KEGG ', '')
        name_str = name_str.replace('TEL PATHWAY', 'TELOMERASE PATHWAY')
        name_str = name_str.replace('RNA PATHWAY', 'PKR SIGNALING PATHWAY')
        outstr = "{ id:" + self.path_id + ", " + "name:'" + name_str \
                 + "', pval:'" + self.p_value + "', size:" + str(self.n_actual) \
                 + ", effective:" + str(
            self.n_effective) + ", url:'" + self.url + "', geneSet: "
        if self.gene_set:
            gene_set_str = "','".join(self.gene_set)
            outstr = outstr + "['" + gene_set_str + "']}"
        else:
            outstr = outstr + "[]}"
        return outstr

    def fetch_url(self):
        """Look up database (refs.pathways) to get info_url for path_id."""
        url_row = None
        cmd = """SELECT info_url FROM refs.pathways WHERE path_id = {};""".format(
            self.path_id)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd)
            rowCount = cur.rowcount
            if not rowCount == 1:
                raise Exception("Result contains %g rows Ids for pathway %s."
                                % (rowCount, self.path_id))
            url_row = cur.fetchone()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        url = url_row[0]
        url = url.split(self.root_url)[1]
        return url


def make_js_file(results_path, out_path):
    # BUILD JS FILE
    allPathways = list()  # will hold pathway objects
    ignoreList = ['CANCER', 'GLIOMA', 'MELANOMA', 'LEUKEMIA', 'CARCINOMA']
    # project_id = raw_input("Project id? ")  # <TODO:sgg> change to input in python3
    # 'pathways_pvalues_{}_pretty.txt'.format(project_id)
    with open(results_path, 'r') as file:
        for line in file:
            vals = line.strip("\n").split('\t')
            pway = PathwaySummaryParsed(vals[0])
            # skip line if part of name is in ignoreList
            skipLine = False
            name = vals[1]
            for ignore in ignoreList:
                if ignore in name:
                    skipLine = True
                    break
            if (skipLine):
                continue
            pway.name = name
            pway.n_actual = int(vals[2])
            pway.n_effective = int(vals[3])
            pway.p_value = vals[4]
            gene_set = set.union(set(eval(vals[6])), set(eval(vals[7])))
            pway.gene_set = gene_set
            allPathways.append(pway)
    with open(out_path, 'w') as out:  # 'pathways_pvalues_{}.js'
        out.write("root_url = '{}';\n".format(pway.root_url))
        out.write("pwayList = [\n")
        for ind, pway in enumerate(allPathways):
            if pway.gene_set:
                out.write(pway.as_string())
                if ind < len(allPathways) - 1:
                    out.write(",")
                out.write("\n")
        out.write("];\n")


def create_svgs(txt_path):
    """Run matlab script that builds svgs in directory containing
    pathways txt."""
    # ORIG cmd = """matlab -nosplash -nodesktop -r "plot_pway_targets('{txtpath}');" < /dev/null >{root_dir}tempstdout.txt 2>{root_dir}tempstderr.txt &"""
    cmd = 'matlab -nosplash -nodesktop -r \"plot_pway_targets({txtpath!r},' \
          '\'--svg\',\'--skipfew\');\"'.format(txtpath=txt_path)
    with open(os.devnull, "r") as fnullin:
        with open(os.devnull, "w") as fnullout:
            subprocess.check_call(cmd, stdin=fnullin, stdout=fnullout,
                                  stderr=fnullout, shell=True)


def create_matrix_svgs(txt_path):
    """Run matlab script that builds matrix svgs. SVGs stored in matrix_
    txt/matrix_svg."""
    # ORIG cmd = """matlab -nosplash -nodesktop -r "plot_pway_targets('{txtpath}');" < /dev/null >{root_dir}tempstdout.txt 2>{root_dir}tempstderr.txt &"""
    cmd = 'matlab -nosplash -nodesktop -r \"plot_patient_genes(' \
          '{txtpath!r});\"'.format(txtpath=txt_path)
    with open(os.devnull, "r") as fnullin:
        with open(os.devnull, "w") as fnullout:
            subprocess.check_call(cmd, stdin=fnullin, stdout=fnullout,
                                  stderr=fnullout, shell=True)


def split_pretty_path(txt_path):
    """Splits pretty.txt path
    (e.g. ORIGDIR/pathways_pvalues_PROJSTR_pretty.txt)."""
    name_prefix = GenericPathwayFileProcessor.base_str  # 'pathways_pvalues_'
    name_postfix = PathwayDetailedFileWriter.name_postfix  # '_pretty.txt'
    path_prefix = txt_path.split(name_postfix)[0]
    orig_dir = path_prefix.split(name_prefix)[0]  # ends with /
    project_str = path_prefix.split(name_prefix)[1]
    js_name = name_prefix + project_str + '.js'
    return orig_dir, project_str, js_name


def generate_files():

    skipfew = True
    # (orig_dir, project_str, js_name) = split_pretty_path(args.path)
    out_dir = current_app.config['UPLOAD_FOLDER']
    project_str = str(UPLOAD.file_id)
    descriptive_name = UPLOAD.filename
    # js_name = <UPLOAD.file_id>.js

    make_js_file(detail_path, js_out_path)
    # html_name = create_html(detail_path, out_dir, project_str, descriptive_name,
    #                         skipfew)
    create_svgs(detail_path)
    create_matrix_svgs(detail_path)




if __name__ == '__main__':
    dbvars = {'host': 'localhost', 'db': 'pway_test',
              'read_default_file': "~/.my.cnf.pway"}
    main()
