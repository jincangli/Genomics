import pandas as pd
import argparse
import subprocess
import os
#import logging
#LOG = logging.getLogger(__name__)

#__author__ = ("Jincang Li",)
#__email__ = ""

parser = argparse.ArgumentParser(description='Stats')
parser.add_argument("--inputChr", "-i", type=str, help="Chr fasta",required=True)
parser.add_argument("--changelist", "-b", type=str, help="",default="",required=True)
parser.add_argument("--outputfa", "-o", type=str, help="output file name eg. output.fasta",default="rename.fasta")
args = vars(parser.parse_args())

chr_file=args["inputChr"]
change_id=args["changelist"]
outputfa=args["outputfa"]
os.system(f'rm {outputfa}')

# with open(chr_file) as chr_fa:
with subprocess.Popen(f'seqkit seq -w 0 {chr_file}', shell = True, stdout=subprocess.PIPE, text=True) as chr_fa:
    chr_fa_dict = {}
    for line in chr_fa.stdout:
        line = line.replace('\n','')
        if line.startswith('>'):
            seq_name = line.split(' ')[0][1:] # 1:去掉'>',2:去掉多余的空格
            chr_fa_dict[seq_name] = ''
        else:
            chr_fa_dict[seq_name] += line.replace('\n','')

endIDlist=list(chr_fa_dict.keys())

change_id_df=pd.read_csv(change_id, sep='\t', header=None)
if change_id_df.shape[1] == 3:
    change_id_df.columns = ['oldId','newId','chain']
elif change_id_df.shape[1] == 2:
    change_id_df.columns = ['oldId','newId']

trantab = str.maketrans('ACGTacgt', 'TGCAtgca')

for i in endIDlist:
    new_chr_fa_dict = {}
    if i not in change_id_df['oldId'].values:
        continue
    changelist_df = change_id_df[change_id_df['oldId']==i]

    index = changelist_df.index.tolist()[0]
    oldId = change_id_df.iloc[index].at['oldId']
    newId = change_id_df.iloc[index].at['newId']
    chain = change_id_df.iloc[index].at['chain'] if len(change_id_df) == 3 else 'none'
    # print(index, oldId, newId,chain, '√')
    if chain == "+":
        seq = str(chr_fa_dict.get(oldId))
        with open(outputfa, "a") as file:
            file.write('>{}\n{}\n'.format(newId,seq))
    elif chain == "-":
        seq = str(chr_fa_dict.get(oldId)).translate(trantab)[::-1]
        with open(outputfa, "a") as file:
            file.write('>{}\n{}\n'.format(newId,seq))
    else:
        seq = str(chr_fa_dict.get(oldId))
        with open(outputfa, "a") as file:
            file.write('>{}\n{}\n'.format(newId,seq))
 
