#!/usr/bin/env python
###
#    possibile citation:
#    ```bibtxt
#    @article{sun2020whole,
#    title={Whole-genome sequencing and bioinformatics analysis of apiotrichum mycotoxinivorans: predicting putative zearalenone-degradation enzymes},
#    author={Sun, Jinyuan and Xia, Yan and Ming, Dengming},
#    journal={Frontiers in Microbiology},
#    pages={1866},
#    year={2020},
#    publisher={Frontiers}
#    }
#    ```
###

from signal import raise_signal
from pywebio import session, config, start_server
from pywebio.output import *
from pywebio.pin import *
from pywebio.input import *

from time import time
import distutils.dir_util
import hashlib
from CMPP_plot_web import *

DATABASE_PATH = '/home/jsun/functional/'
TEST_MODE = 0

def check_file_format(input_file_path):
    name_list = []
    seq_list = []
    with open(input_file_path, 'r') as ifile:
        for line in ifile:
            if line.startswith('>'):
                name_list.append(line.strip())
                seq_list.append('')
            else:
                try:
                    seq_list[-1] += line.strip()
                except IndexError:
                    return 0
    all_aa = ''.join(seq_list)
    if len(all_aa) < 1:
        put_text("No amino acid found! Check your input file").style('color: red; font-size: 40px')
        return 0
    else:
        seq_num = len(name_list)
        put_text(f'There are {seq_num} protein sequences in your file.').style('color: grey; font-size: 40px')
        return 1


def timer_func(func):
    # This function shows the execution time of 
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        func(*args, **kwargs)
        t2 = time()
        put_text(f'Sequences search done in {(t2-t1):.4f}s').style('color: black; font-size: 30px')
        return 0
    return wrap_func

@timer_func
def run_system(cmd:str, db_name:str="db"):
    if TEST_MODE:
        pass
    else:
        os.system(cmd)
        print(cmd)
    return cmd

class CMPP_web():
    def __init__(self, input_file:str, database_path:str, output_path:str = 'CMPP_out/') -> None:
        self.token = input_file.replace(input_file.split("/")[-1],',')[:-1]
        self.output_path = input_file.replace(input_file.split("/")[-1],output_path)
        distutils.dir_util.mkpath(self.output_path)
        self.input_file_path = input_file
        self.input_file_name = input_file.split("/")[-1]
        self.database_path = database_path
        self.fixed_blast_param = '-evalue 1e-10 -outfmt 6 -max_target_seqs 5 -num_threads 4 -out'
        
    def run_search(self):
        # TODO: parsing file with python code.
        job_list = [
            f"blastp -query {self.input_file_path} -db {self.database_path}/merops/merops_scan.lib {self.fixed_blast_param} {self.output_path}/{self.input_file_name}.merops.btb",
            f"blastp -query {self.input_file_path} -db {self.database_path}/PHI/phi-base_current.fas {self.fixed_blast_param} {self.output_path}/{self.input_file_name}.phi.btb",
            f"hmmscan -o P450.out --tblout {self.output_path}/{self.input_file_name}.p450.htb --noali --cpu 2 -E 1e-5 {self.database_path}/P450/P450.hmm.txt {self.input_file_path}",
            f"hmmscan -o CAZy.out --tblout {self.output_path}/{self.input_file_name}.cazy.htb --noali --cpu 2 -E 1e-5 {self.database_path}/CAZy/dbCAN-HMMdb-V9.txt {self.input_file_path}"
        ]
        # t0 = time.time()
        put_text('Running Blast search against merops database...').style('color: black; font-size: 30px')
        run_system(job_list[0])
        put_text('Running Blast search against PHI database...').style('color: black; font-size: 30px')
        run_system(job_list[1])
        put_text('Running hmmscan search against P450 database...').style('color: black; font-size: 30px')
        run_system(job_list[2])
        put_text('Running hmmscan search against CAZy database...').style('color: black; font-size: 30px')
        run_system(job_list[3])

        distutils.dir_util.mkpath(self.output_path)
        os.system("grep -v \"#\" %s.p450.htb|awk -F \" \" '{print$3}'|sort|uniq > %sp450.glist"%(self.output_path+self.input_file_name, self.output_path))
        os.system("grep -v \"#\" %s.cazy.htb|awk -F \" \" '{print$3}'|sort|uniq > %scazy.glist"%(self.output_path+self.input_file_name, self.output_path))
        os.system("awk '{print$1}' %s.merops.btb |sort|uniq > %smerops.glist"%(self.output_path+self.input_file_name, self.output_path))
        os.system("awk '{print$1}' %s.phi.btb |sort|uniq > %sphi.glist"%(self.output_path+self.input_file_name, self.output_path))

    
    def run_plot(self):
        path = self.output_path
        dataset_dict = read_glist(path=self.output_path)
        venn_plot(dataset_dict,path=path)

        phi_dict = make_phi_dict(btbfile = self.output_path+self.input_file_name+'.phi.btb')
        phi_plot(phi_dict, short_phi_names = True,  dark = False, path=path)

        merops_map_dict = read_merops_map(merops_path = "common/merops.txt")
        merops_dict = mk_merops_dict(merops_map_dict, btbfile = self.output_path+self.input_file_name+'.merops.btb')
        merops_plot(merops_dict, short_merops_names = False, path=path)

        cazy_dict = mk_cazy_dict(htbfile = self.output_path+self.input_file_name+'.cazy.htb')
        cazy_plot(cazy_dict, short_cazy_names=True, path=path)
        os.system(f'touch {self.token}/DONE')


def input_seq():
    put_text("Protein sequence:")
    res = textarea('Protein sequence', type=TEXT)
    md5 = hashlib.md5()
    md5.update(res.encode('utf-8'))
    hash_code = md5.hexdigest()
    # TODO: if hash_code in database, return results 
    with open(f'{hash_code}.fa', 'w+') as ifile:
        ifile.write(res)
    put_text(f"your jobs {hash_code} was sumbmitted.")
    return hash_code+".fa"

def creat_token(text):
    md5 = hashlib.md5()
    md5.update(text)
    token = md5.hexdigest()
    return token

def get_genome():
    f = file_upload("Upload a fungal genome for annotation", 
                    max_size='100M',
                    help_text="File show be in fasta file format.")
    print("Loading file")
    token = creat_token(f['content'])
    print("Creating token")
    if token in os.listdir():
        put_text('The genome is uploaded!\nYour job id is '+token).style('color: black; font-size: 40px')
        if 'DONE' in os.listdir(token):
            put_text('get your results from '+token)
        else:
            put_text('A job with the same input is running! Please be patient!')
        return token, f['filename']
    else:
        distutils.dir_util.mkpath(token)
        open(token+"/"+f['filename'], 'wb').write(f['content'])
        check = check_file_format(token+"/"+f['filename'])
        while not check:
            os.system(f'rm -rf ./{token}')
            put_text(f'The uploaded file is invalid! Upload again please!').style('color: red; font-size: 40px')
            f = file_upload("Upload a fungal genome for annotation", 
                    max_size='100M',
                    help_text="File show be in fasta file format.")

            token = creat_token(f['content'])
            check = check_file_format(token+"/"+f['filename'])

        put_text(f'The genome is uploaded!\nYour job id is {token}').style('color: black; font-size: 40px')
        return token+"/"+f['filename']


def run_app(input_file, token, DATABASE_PATH):
    cmpp = CMPP_web(input_file, DATABASE_PATH)
    cmpp.run_search()
    cmpp.run_plot()
    os.system(f'rm {input_file}')
    os.system(f"tar vczf {token}.tar.gz {token}")


def load_results(token):
    
    put_text("Here are your results:").style('text-align:center;font-size:40px')
    for file in os.listdir(f'{token}/CMPP_out/'):
        if file.endswith("png"):
            img = open(f'{token}/CMPP_out/{file}', 'rb').read()  
            name = file.replace("_",' ').replace(".png","")
            put_text(f'{name} plot').style('font-size:30px')
            put_image(img).style('display: block;margin-left: auto;margin-right: auto;width: 60%;')
    content = open(f'./{token}.tar.gz', 'rb').read()
    put_file(f'{token}.tar.gz', content, 'Download all results').style('text-align:center, font-size:40px')

def run_job():
    input_file = get_genome()
    print("Genome loaded!")
    if "/" in input_file:
        # 2. run the app
        token = input_file.split("/")[0]
        run_app(input_file, token, DATABASE_PATH)
        load_results(token)
    else:
        token, file_name = input_file
        if 'DONE' in os.listdir(token):
            load_results(token)
        else:
            pass

def load_job():
    token = input("Your token:",help_text="Demo token: ce428d17e648815ce5cd305c77751a91")
    load_results(token)

def main():
    put_markdown('''
    # CMPP: Annotate CAZy, MEROPS, PHI, P450 for fungi
    To annotate carbohydrate-active enzymes, proteases and inhibitors, pathogen-host interaction proteins, and Cytochrome P450 in a fungi proteome.
    [GitHub](https://github.com/JinyuanSun/CMPP)
    [Colab version](https://colab.research.google.com/github/JinyuanSun/CMPP/blob/main/ColabCMPP.ipynb)
    ''')
    img = open('./logo_CMPP.png', 'rb').read()
    put_image(img).style('display: block;margin-left: auto;margin-right: auto')
    put_buttons(['Load a job', 'Start a new job'], onclick=[load_job, run_job])
    

if __name__ == '__main__':
    # load_results('.')
    # main()
    start_server(main, session_expire_seconds=600, port=3000, debug=True, cdn=False, max_payload_size='100M')

