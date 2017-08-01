#!/usr/bin/env python3.6
import argparse
import os
import sys
from collections import defaultdict
check = 'annot_0.0'
blast_d = 'Blast_tmp'
miRNA_d = 'forMiRNA'
out_d = 'blast_out_result_after_check'
unm_d = 'unmatch_AfterBlastmiRNA_formirplant'
exp_d = 'exp_files'


def args_parse():
    parser = argparse.ArgumentParser(description='For parse fa to miRNA or degradome')
    parser.add_argument('-q', '--query', type=float, default=1.0, help='Query coverage threshold, default is 0.9')
    parser.add_argument('-t', '--target', type=float, default=1.0, help='Target coverage threshold, default is 0.9')
    parser.add_argument('-i', '--identity', type=float, default=1.0, help='Identity threshold, default is 0.9')
    args = parser.parse_args()
    return args


def get_osa_dict(path):
    mirna_seq_dict = {}
    seq_mirna_dict = defaultdict(list)
    with open(path) as osa_fa:
        for line in osa_fa:
            if line.startswith('>'):
                mirna = line[1:].split(' ')[0]
                seq = next(osa_fa).rstrip()
                mirna_seq_dict[mirna] = seq
                seq_mirna_dict[seq].append(mirna)
    return mirna_seq_dict, seq_mirna_dict


def unmatched_get(read_list, rp):
    r_set = set(read_list)
    all_c, out_c = 0, 0
    with open(os.path.join(miRNA_d, rp + '.collapse_m.fa')) as in_fa, \
            open(os.path.join(unm_d, rp + '.unmatched_collapse_m.fa'), 'w') as out_fa:
        for line in in_fa:
            if line.startswith('>'):
                all_c += 1
                read_id = line[1:].rstrip()
                seq = next(in_fa)
                if read_id not in r_set:
                    out_c += 1
                    out_fa.write(line + seq)
        print('{}\tall: {}\tout: {}\tblast: {}\tcheck: {}'.format(rp, all_c, out_c, len(r_set), all_c - out_c))


def blast_get(b_file, qc, tc, idt, rp):
    read_list = []
    mirna_read_list_dict = defaultdict(list)
    with open(os.path.join(blast_d, b_file)) as in_a, open(os.path.join(out_d, rp + '.osa_miRNA_{}_{}_[]'.format(
            qc, tc, idt)), 'w') as out_a:
        out_a.write('QueryID\tQ_length\tTargetID\tT_length\tQ_st\tQ_et\tT_st\tT_et\tE-value\tBit-score\tGap_open\t'
                    'Match\tAlign\tQ_coverage\tT_coverage\tIdentity\n')
        for line in in_a:
            lf = line.rstrip().split('\t')
            query_cov = int(lf[12]) / int(lf[1])
            target_cov = int(lf[12]) / int(lf[3])
            identity = int(lf[11]) / int(lf[12])
            if (query_cov >= qc) and (target_cov >= tc) and (identity >= idt):
                target = lf[2].split(' ')[0]
                read_list.append(lf[0])
                mirna_read_list_dict[target].append(lf[0])
                o_tmp = lf[0:13] + list(map(str, [query_cov, target_cov, identity]))    # map is a object, need list()
                out_a.write('\t'.join(o_tmp) + '\n')
        return read_list, dict(sorted(mirna_read_list_dict.items()))


def get_exp(br_list, m2rl_dict, m2s_dict, rp):
    with open(os.path.join(exp_d, rp + '_exp_with_mean.csv'), 'w') as exp_wm, \
            open(os.path.join(exp_d, rp + '_exp_no_mean.csv'), 'w') as exp_nm:
        mir_count = 0
        for miRNA, r_list in m2rl_dict.items():
            mir_count += 1
            seq = m2s_dict.get(miRNA)
            count_wm, count_nm = 0, 0
            for r_ID in r_list:
                r_count = int(r_ID.split('-')[1])
                count_nm += r_count
                count_wm += (r_count / br_list.count(r_ID))
            exp_nm.write('{},{},{}\n'.format(miRNA, count_nm, seq))
            exp_wm.write('{},{},{}\n'.format(miRNA, count_wm, seq))
            if mir_count <=5:
                print('[check]{}\t{}\t{}\t{}\t{}\t{}'.format(rp, miRNA, count_nm, count_wm, r_list, seq))


def main(argv=None):
    if argv is None:
        argv = args_parse()
        mirna2seq_dict, seq2mirna_dict = get_osa_dict(blast_d + '/osa_mature.fa')   # for expression output file
        print('miRNA2Seq_dict: {}\tSeq2miRNA: {}'.format(len(mirna2seq_dict), len(seq2mirna_dict)))
        with open('/'.join([blast_d, 'osa_seq2miRNA']), 'w') as seq2m:
            seq2m.write('miRNA_sequence\tmiRNA_id\n')
            for k, v in seq2mirna_dict.items():
                seq2m.write(k + '\t' + ','.join(v) + '\n')
        # Start to process blast results
        file_list = os.listdir(blast_d)
        file_count = 0
        for file in file_list:
            if file.endswith(check):
                file_count += 1
                rice_p = file.replace('.osa_mature.fa.w8e-3.annot_0.0', '')
                b_read_list, mr2read_dict = blast_get(file, argv.query, argv.target, argv.identity, rice_p)
                print('{}\tmir2ReadList_dict: {}\tb_reads_list: {}'.format(rice_p, len(mr2read_dict), len(b_read_list)))
                unmatched_get(b_read_list, rice_p)  # get unmatched fa
                get_exp(b_read_list, mr2read_dict, mirna2seq_dict, rice_p)  # get exp file

        # open(os.path.join(out_d, rp + '.osa_miRNA_exp_{}_{}_[]'.format(qc, tc, idt)), 'w') as out_e
if __name__ == '__main__':
    sys.exit(main())
