#! /usr/bin/env python2.7

from __future__ import division
import os,sys,re,argparse,gzip,commands


def get_analysis_option():
    my_cmdparser = argparse.ArgumentParser(description='transfrom fq reads to tag with numbers\n')
    my_cmdparser.add_argument('-o',action='store',dest='outfile',type=str,help='outfile')
    my_cmdparser.add_argument('-i',action='store',dest='infile',help='in fq files')
    my_cmdparser.add_argument('-c',action='store',dest='cut_off',default = 1, type = int, help='cut off for reads count')
    options=my_cmdparser.parse_args()
    configure={i:getattr(options,i) for i in dir(options) if (i[0] != '_' and i not in ('ensure_value', 'read_file','read_module'))}
    return configure

def fq2tag(configure):
    fq_len_frequency_dict = {}
    tag_dict = {}
    reads_n = 0
    out_reads_n = 0
    file_type = 'txt'
    fq_list = configure['infile'].strip().split(',')
    len_stat_lst = [0,0,0,0]
    for i in fq_list:
        file_t = commands.getoutput('file %s' % (i))
        if re.search(r'gzip compressed data',file_t):
            file_type='gz'
        elif re.search(r'ASCII text',file_t):
            file_type='txt'
        if file_type =='gz':
		    f=gzip.open(i,'rb')
        elif file_type == 'txt':
	        f=open(i,'r')
        for n,i in enumerate(f.readlines()):
            line_number=n+1
            if line_number % 4 !=2:continue
            fq_seq = i.strip()
            try:
	            tag_dict[fq_seq]+= 1
            except KeyError:
                tag_dict[fq_seq] = 1
    sorted_tag_number=sorted(tag_dict.keys(),key=lambda z:tag_dict[z],reverse=True)
    with open(configure['outfile'],'w') as of:
        for n,i in enumerate(sorted_tag_number):
            fq_len = len(i)
            tag_count = tag_dict[i]
            if fq_len in fq_len_frequency_dict:
                fq_len_frequency_dict[fq_len] += tag_count
            else:
                fq_len_frequency_dict[fq_len] = tag_count
            reads_n += tag_count
            tag_number = n+1
            tag_number = 't{0:0>8d}'.format(tag_number)
            if tag_count > configure['cut_off']:
                out_reads_n += tag_count
                of.write('>{0}\t{1:d}\n{2}\n'.format(tag_number,tag_count,i))
    with open(configure['outfile']+'.freq.stat','w') as of:
        for i in fq_len_frequency_dict:
            len_n = fq_len_frequency_dict[i]
            len_freq = len_n/reads_n
            if i <15:
                len_stat_lst[0] += len_freq
            elif 15<=i<30:
                len_stat_lst[1] += len_freq
            elif 30<= i <=40:
                len_stat_lst[2] += len_freq
            else :
                len_stat_lst[3] += len_freq
            of.write('{0}\t{1}\t{2:f}\n'.format(i,len_n,len_freq))
        flag = 'ok' if len_stat_lst[1]+ len_stat_lst[2] > 0.5 else 'no'
        of.write('{0}\t{1}\t{2:d}\t{3:.2f}\t{4:d}\n'.format(flag,'\t'.join([ str('{0:.2f}'.format(j)) for j in len_stat_lst]),out_reads_n, out_reads_n/reads_n, reads_n))

def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        get_analysis_option()
    else :
        configure=get_analysis_option()
        fq2tag(configure)


if __name__ == '__main__':
	main()