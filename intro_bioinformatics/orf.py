# Long Open Reading Frame in DNA
from sys import *
import operator
min_length = 125 # minimum length of long orf (in #codon)
output_file_name = 'personalORF.txt'

# my own prokaryote
fna_files = ['NC_010622.fna',
             'NC_010623.fna',
             'NC_010625.fna',
             'NC_010627.fna'
             ]
ptt_files = ['NC_010622.ptt',
             'NC_010623.ptt',
             'NC_010625.ptt',
             'NC_010627.ptt'
             ]

# community prokaryote
"""fna_files = ['NC_004347.fna',
             'NC_004349.fna'
             ]
ptt_files = ['NC_004347.ptt',
             'NC_004349.ptt'
             ]"""
#****************end of hw*****************

# 06's community prokaryote: Bartonella quintana
"""fna_files = ['NC_005955.fna', 'NC_005955.fna']
ptt_files = ['NC_005955.ptt', 'NC_005955.ptt']"""

#***************end of test****************

# define the start codon, end codon and translation table
start_codon = ['ATG', 'GTG']
stop_codon = ['TAA', 'TAG', 'TGA']
trans = {'A':'T',
         'T':'A',
         'C':'G',
         'G':'C'
         }

# statisticall figures
TP = 0
sTP = 0
FN = 0
FP = 0
B = 0
A = 0


def translate(codon):
    if len(codon)!=3:
        print 'Error: size of codon must be 3'
        return ''
    else:
        for i in [0,1,2]:
            if codon[i] not in trans:
                print 'Error: data contains noise', codon[i]
                return ''
            codon[i] = trans[codon[i]]
        codon.reverse()
        return ''.join(codon)

if __name__ == '__main__':
    # data pre-processing
    # step 0: setup the final output file name
    fout = open(output_file_name, 'w')

    # step 1.1: initialize protein informations included
    #           in the .ptt files
    protein_coding = len(ptt_files) * [{}] # all the positions of proteins listed in .ptt files

    # step 1.2: read in each .ptt file
    for curr_file, ptt in enumerate(ptt_files):
        file_in = open(ptt, 'rt')
        for line_num, line in enumerate(file_in):
            line = line.split()
            # the first 3 lines in .ptt files are useless, discard
            if line_num<3:
                continue
            start_loc = ''
            stop_loc = ''
            if (line[1] == '+'):
                start_loc = ''.join(line[0].split('..')[0])
                stop_loc = ''.join(line[0].split('..')[1])
            else:
                stop_loc = ''.join(line[0].split('..')[0])
                start_loc = ''.join(line[0].split('..')[1])
            protein_coding[curr_file][stop_loc] = start_loc
    print '.ptt data pre-processing finished'

    # step 2.1: initialize DNA sequence
    #         add to each sequence a non-sense symbol
    #         to make the index consistent with .ptt files
    sequence = len(fna_files) * [['#']]

    # step 2.2: read in a DNA file and do the calculation
    for curr_file, fna in enumerate(fna_files):
        # store the data in main memory
        file_in = open(fna, 'rt')
        for line_num, line in enumerate(file_in):
            # the first line is use less, discard it
            if line_num < 1:
                continue
            sequence[curr_file] += line[:-1]
        print 'DNA:%d, .fna data pre-processing finished' % (curr_file+1)
        size = len(sequence[curr_file])

        # setup the long open reading frame collector
        lorf = []

        # find long open reading frame in plus strand
        for start_idx in [1,2,3]:
            idx = start_idx
            # plus strand open read frame collector
            stat = ['null', # current state: null, +, -, both
                    0,      # length of current orf in plus strand (in #codon)
                    0       # start position in plus strand
                    ]

            while idx+2<size:
                # plus strand
                codon = ''.join(sequence[curr_file][idx:idx+3])

                # case 1: no current orf in plus strand
                if stat[0] == 'null':
                    # start new orf in plus strand
                    if codon in start_codon:
                        stat[0] = 'start'
                        stat[1] = 1
                        stat[2] = idx

                # case 2: in some orf in some open frame
                else:
                    if codon in stop_codon:
                        if stat[1]+1 >= min_length:
                            start_loc = str(stat[2])
                            stop_loc = str(stat[2]+3*stat[1]+2)
                            B += 1
                            if stop_loc in protein_coding[curr_file]:
                                if protein_coding[curr_file][stop_loc] == start_loc:
                                    TP += 1
                                else:
                                    sTP += 1
                            # add entry to long open reading frame collector
                            lorf += [[int(start_loc), int(stop_loc), '+']]

                        stat[0] = 'null'
                        stat[1] = 0
                        stat[2] = 0
                    else:
                        stat[1] += 1

                # start next iteration
                idx += 3
            print 'DNA:%d, +strand, offset:%d, finished' %( curr_file+1, start_idx )

        # find long open reading frame in minus strand
        for start_idx in [size-1, size-2, size-3]:
            idx = start_idx
            # minus strand open read frame collector
            stat = ['null', # current state
                    0,      # length of current orf in minus strand (in #codon)
                    0       # start position in minus strand
                    ]

            while idx-2>0:
                # minus strand
                codon = translate(sequence[curr_file][idx-2:idx+1])

                # case 1: no current orf in minus strand
                if stat[0] == 'null':
                    # start new orf in minus strand
                    if codon in start_codon:
                        stat[0] = 'start'
                        stat[1] = 1
                        stat[2] = idx

                # case 2: in some orf in minus strand
                else:
                    if codon in stop_codon:
                        if stat[1]+1 >= min_length:
                            stop_loc =  str(stat[2]-3*stat[1]-2)
                            start_loc = str(stat[2])
                            B += 1
                            if stop_loc in protein_coding[curr_file]:
                                if protein_coding[curr_file][stop_loc] == start_loc:
                                    TP += 1
                                else:
                                    sTP += 1
                            # add entry to long open reading frame collector
                            lorf += [[int(stop_loc), int(start_loc), '-']]

                        stat[0] = 'null'
                        stat[1] = 0
                        stat[2] = 0
                    else:
                        stat[1] += 1

                # start next iteration
                idx -= 3
            print 'DNA:%d, -strand, offset:%d, finished' % ( curr_file+1, size-start_idx)


        # step 3.2: sort the long open reading frame and output the result
        #         with format: NC+004347 13076 13744 +
        lorf.sort(key=operator.itemgetter(0))
        curr_file_name = fna[0:-4]
        print curr_file_name
        for line in lorf:
            line[0] = str(line[0])
            line[1] = str(line[1])
            print>>fout, curr_file_name, ' '.join(line)



    # statistical report
    for item in protein_coding:
        A += len(item)
    A, B = float(A), float(B)
    FN = A - TP - sTP
    FP = B - TP - sTP
    print 'Statistics as follows: '
    print ' A:%d\n B:%d\n TP:%d\n sTP:%d\n FN:%d\n FP:%d\n' %(A, B, TP, sTP, FN, FP)
    Sn = TP/A
    sSn = (TP + sTP)/A
    FOR = FN/A
    PPV = TP/B
    sPPV = (TP + sTP)/B
    FDR = FP/B
    print ' Sn:%f\n sSn:%f\n FOR:%f\n PPV:%f\n sPPV:%f\n FDR:%f\n' %(Sn, sSn, FOR, PPV, sPPV, FDR)
    print>>fout, int(TP), int(sTP), int(FN), int(FP), Sn, sSn, FOR, PPV, sPPV, FDR






