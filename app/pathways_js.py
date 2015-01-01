#!/Users/Stephen/Library/Enthought/Canopy_64bit/User/bin/python
import os
import subprocess


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
    print cmd
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
    print cmd
    with open(os.devnull, "r") as fnullin:
        with open(os.devnull, "w") as fnullout:
            subprocess.check_call(cmd, stdin=fnullin, stdout=fnullout,
                                  stderr=fnullout, shell=True)


def generate_files(out_dir, detail_path, user_upload):

    # skipfew = True
    descriptive_name = str(user_upload.file_id) + user_upload.filename
    js_out_path = os.path.join(out_dir, descriptive_name + ".js")

    make_js_file(detail_path, js_out_path)
    # html_name = create_html(detail_path, out_dir, project_str, descriptive_name,
    #                         skipfew)
    create_svgs(detail_path)
    create_matrix_svgs(detail_path)




if __name__ == '__main__':
    generate_files()
