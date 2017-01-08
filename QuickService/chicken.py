#!/usr/bin/env python3.5
"""
Created by JRLi on 2017/01/08 for chicken contigs
"""
import sys

sequence_file = sys.argv[1]
target_name = sys.argv[2]
start_site = int(sys.argv[3])
end_site = int(sys.argv[4])
id, seq = None, []
id2seq_dict = {}
with open('./'+sequence_file, 'r') as input_file, open('./'+sys.argv[2]+'_'+sys.argv[3]+'_'+sys.argv[4]+'.fa', 'w') as out:
    for line in input_file:
        if line.startswith('>'):
            if id is not None:
                sequence = ''.join(seq)
                id2seq_dict[id] = sequence
                seq = []
            id = line.replace('\r\n', '').replace('\n', '')[1:]
        else:
            seq.append(line.replace('\r\n', '').replace('\n', ''))
    sequence = ''.join(seq)
    id2seq_dict[id] = sequence
    out.write('>'+target_name+'\n'+id2seq_dict[target_name][start_site:end_site + 1]+'\n')
