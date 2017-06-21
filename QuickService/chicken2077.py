def reverse_complement(line, seq = False):
    complement_dict = {'a':'t','t':'a','c':'g','g':'c','A':'T','T':'A','C':'G','G':'C'}
    r_line = line[::-1]
    if seq:
        rc_line = "".join(complement_dict.get(bp, bp) for bp in r_line)
        return rc_line
    else:
        return r_line


def insert_end(line, index):
    return "".join([line[j:j+index]+"\n" for j in range(0, len(line), index)])


def get_2077():
    with open('Annotation_2077.gff','r') as input, open('Annotation_c2077.gff','w') as output:
        for line in input:
            if line.startswith('2077'):
                output.write(line)



def rev_fa():
    seq_list = []
    with open('E://vb1/contig2077.fa') as input, open('E://vb1/contig2077_r.fa', 'w') as output:
        for line in input:
            line = line.replace('\r\n','').replace('\n', '').replace('\"', '').strip()
            if line.startswith('>'):
                output.write(line + '\n')
            else:
                seq_list.append(line)
        seq = ''.join(seq_list)
        print(len(seq))
        seqr = reverse_complement(seq, True)
        print(len(seqr))
        seqrb = insert_end(seqr, 60)
        output.write(seqrb)


def gtfr_2077():
    s_dict = {'+':'-','-':'+'}
    with open('E://vb1/2077.gff') as input, open('E://vb1/2077_r.gff', 'w') as output:
        for line in input:
            line = line.replace('\r\n', '').replace('\n', '')
            lf = line.split('\t')
            lf[6] = s_dict.get(lf[6], lf[6])
            line_c = '\t'.join(lf)
            output.write(line_c + '\n')


def rv_gtf():
    seqlen = 361966
    with open('E://vb1/2077_r.gff') as input, open('E://vb1/2077_rv.gff', 'w') as output:
        for line in input:
            line = line.replace('\r\n', '').replace('\n', '')
            lf = line.split('\t')
            lf[3] = str(seqlen - int(lf[3]) + 1)
            lf[4] = str(seqlen - int(lf[4]) + 1)
            lf[3], lf[4] = lf[4], lf[3]
            line_c = '\t'.join(lf)
            output.write(line_c + '\n')


def main():
    rv_gtf()

if __name__ == '__main__':
    main()