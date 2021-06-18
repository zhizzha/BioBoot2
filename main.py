from Bio import SeqIO
import argparse as ap
import sys
from defs import *
#TODO Выбрать и скомпонавать последовательности по Segment
#TODO Выравнять последовательности друг на друга
#TODO Посчитать консервативность каждого учатска
#TODO (Сит) Провести анализ BLAST
#TODO Красиво оформить :> (GUI)





def createParser ():
    parser = ap.ArgumentParser()
    parser.add_argument ('-f', '--file', default='aln.fasta', help='File name .fasta')
    return parser


if __name__ == '__main__':
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])
    with open(namespace.file) as handle:
        reqs = list(SeqIO.parse(handle, "fasta"))
        segments = list(set([record.id.split('|')[0] for record in reqs]))
        print("\nSegments in file: {}".format('/'.join(segments)))
        print("Segment count: {}\n".format(len(segments)))
        for segment in segments:
            print("Processing segment {}...".format(segment))
            records = []
            for record in reqs:
                if segment == record.id.split('|')[0]:
                    records.append(record)
            print("Records count: {}\n".format(len(records)))
            SeqIO.write(records,segment+'.fasta','fasta')
