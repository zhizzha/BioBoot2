from Bio import SeqIO,AlignIO
from Bio.Align import AlignInfo
from Bio.SubsMat import FreqTable
import argparse as ap
import sys,os
from defs import *
from classes import Frame
#TODO Посчитать консервативность каждого учатска
#TODO (Сит) Провести анализ BLAST
#TODO Красиво оформить :> (GUI)
segments = ""




def createParser ():
    parser = ap.ArgumentParser()
    parser.add_argument ('-f', '--file', default='aln.fasta', help='File name .fasta')
    return parser


if __name__ == '__main__':
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])
    #Processing raw fasta into row segments
    with open(namespace.file) as handle:
        reqs = list(SeqIO.parse(handle, "fasta"))
        segments = list(set([record.id.split('|')[0] for record in reqs]))
        print("\nSegments in file: {}".format('/'.join(segments)))
        print("Segment count: {}".format(len(segments)))
        for segment in segments:
            print("Processing segment {}...".format(segment))
            records = []
            for record in reqs:
                if segment == record.id.split('|')[0]:
                    records.append(record)
            print("Records count: {}\n".format(len(records)))
            SeqIO.write(records,segment+'.fasta','fasta')
    """
    #alaignment
    files = [i + '.fasta' for i in segments]
    print("List of files to align: \n{}\n".format(' | '.join(files)))

    for file in files:
        print("aligning {}".format(file))
        command = 'mafft --quiet --auto --thread -1 {0} > ./aln/aln_{0}'.format(file)
        os.system(command)
    """
    #conservation
    files = ['./aln/aln_{}.fasta'.format(i) for i in segments]
    print("List of files to processing: \n{}\n".format(' | '.join(files)))
    seqs = []
    for file in files:
        aln = AlignIO.read(file,"fasta")
        summary_align = AlignInfo.SummaryInfo(aln)
        summary_align1 = AlignInfo.SummaryInfo(aln)
        for i in range(20,21):
            consensus = summary_align.dumb_consensus()
            c = 0
            for j in range(0,len(consensus) - i):
                frame = consensus[j:j+i]
                if not("X" in frame):
                    expect_freq = FreqTable.FreqTable({"a": 0.25, "g": 0.25, "t": 0.25, "c": 0.25},FreqTable.FREQ)
                    info_content = summary_align.information_content(
                        j, j+i, e_freq_table=expect_freq, chars_to_ignore=["N","X","y","r"]
                    )
                    if info_content/i >= 2:
                        #print('Seq : {0} Score: {1}'.format(frame,info_content/i))
                        print(frame)
                        seqs.append(frame)
                        c +=1
            print("Count for {0} in {1}: {2}".format(file,i,c))
    print(len(seqs))
    seqs = list(set(seqs))
    print(len(seqs))
