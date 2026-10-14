import os
import subprocess

def run():
    workdir = os.path.abspath(os.getcwd())
    os.chdir("..")

    fastadir = os.path.abspath(os.getcwd())

    filenames = []

    for file in os.listdir(fastadir):
        if file.endswith(".fasta"):
            filenames.append(fastadir + "/" + file)

    filenames.sort()

    os.chdir(workdir)

    for file in filenames:
        dir_out = file.split("subset_")[1].split(".fasta")[0] + "_fimo"
        meme_folder = "../meme/" + file.split("subset_")[1].split(".fasta")[0] + "_rem_meme"
        meme = meme_folder + "/meme.txt"
        proc = subprocess.Popen("grep 'Description' " + meme, shell=True, stdout=subprocess.PIPE, text=True)
        res = proc.communicate()[0].split('\n')
        res.remove("")
        for motif in res:
            motif_id = motif.split(" ")[1]
            motif_no = motif.split(" ")[2]
            os.system(f'fimo -o {dir_out}_{motif_id}_{motif_no} --motif {motif_id} {meme_folder}/meme.txt {file}')


#        os.system(f"")

if __name__ == "__main__":
    run()
