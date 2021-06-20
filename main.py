from Bio import SeqIO,AlignIO
from Bio.Align import AlignInfo
from Bio.SubsMat import FreqTable
import argparse as ap
import sys,os
from defs import *
#TODO (Сит) Провести анализ BLAST
#TODO Красиво оформить :> (GUI)
segments = ""




def createParser ():
    parser = ap.ArgumentParser()
    parser.add_argument ('-f', '--file', default='aln.fasta', help='File name .fasta')
    parser.add_argument ('-t', '--treeshold', default='0.98', help='Conservation threeshold')
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
        consensus = summary_align.dumb_consensus()
        pssm = summary_align.pos_specific_score_matrix(consensus, chars_to_ignore=["N"])
        max = len(summary_align.get_column(0))
        cons = []
        th = float(namespace.treeshold)
        for pos in range(aln.get_alignment_length()):
            max_percent = 0
            for letter in pssm[pos].keys():
                percent = pssm[pos][letter]/max
                if percent > max_percent:
                    max_percent = percent
                    l = letter
            if max_percent > th:
                cons.append(l)
            else:
                cons.append("N")
        cv_cons = ''.join(cons)
        arr = cv_cons.split("N")
        c = 0
        for i in range(len(arr)):
            if len(arr[i]) > 27:
                c +=1
                seqs.append((len(arr[i]),arr[i]))

    [print(i) for i in seqs]
