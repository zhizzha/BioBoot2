from Bio import SeqIO
import ap as ap
import sys
from .defs import *
#TODO Выбрать и скомпонавать последовательности по Segment
#TODO Выравнять последовательности друг на друга
#TODO Посчитать консервативность каждого учатска
#TODO (Сит) Провести анализ BLAST
#TODO Красиво оформить :> (GUI)





def createParser ():
    parser = ap.ArgumentParser()
    parser.add_argument ('-f', '--file', default='aln.fasta')
    return parser
 
 
if __name__ == '__main__':
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])
    with open(namespace.name) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            print(record.id)   
        print()
        print(groop(SeqIO.parse(handle, "fasta"),"1"))
        
    