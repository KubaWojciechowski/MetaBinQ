import os
import sys
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def split_fragment(sequence, minl, maxl):
    """
    Split sequence into fragments no longer than n
    
    Inputs:
    sequence - sequence
    n - max length of fragment
    
    Output - fragments
    """
    fragments = []
    lengths = []
    i = 0
    while i < len(sequence):

        tmp = random.randint(minl, maxl)
        seq = sequence[i:i + tmp]
        fragments.append(seq)
        lengths.append(len(seq))
        i += len(seq)

    return(fragments, lengths)


def contaminate(contamination, genome, fraction):

    ref_len = len(genome)
    cont_records = []
    cont = 0
    cont_len = 0

    while cont <= fraction:
        
        tmp = random.randint(0,len(contamination))
        rand_len = random.randint(5000, 50000)
        fragment = contamination[tmp:tmp+rand_len]
        cont_len += len(fragment)
        cont_records.append(SeqRecord(fragment, id="contamination"))

        cont = cont_len/(ref_len+cont_len)

    return(cont_records, cont)


def compleatness(genome, level):
    
    frag_list, frag_lengths = split_fragment(genome.seq, 5000, 50000)
    total_len = sum(frag_lengths)

    comp = 1

    while comp >= level:
        tmp = random.randint(0,len(frag_list)-1)
        

        del frag_list[tmp]
        del frag_lengths[tmp]
        
        comp = sum(frag_lengths)/total_len

    # fragments to seqrecords
    records = []
    for i in frag_list:
        records.append(SeqRecord(i, id="fragment"))

    return(records, comp)    



if __name__ == '__main__':

    path = "./genomes/"
    n_bins = 200
    genomes = os.listdir(path)
    log = open("vscont.csv", 'w')
    print("bin,comp,cont", file=log)

    main = random.choices(genomes, k=n_bins)

    n = 0
    while n < n_bins:

        name = path+random.sample(genomes, 1)[0]
        record = list(SeqIO.parse(name, 'fasta'))[0]
        ref_len = len(record.seq)

        if 100000000 > ref_len > 800000:
        
            cont_name = path + random.sample(genomes, 1)[0]
            contamination = list(SeqIO.parse(cont_name, 'fasta'))[0]
            
            if len(contamination.seq) > 10000:

                if contamination.id != record.id:
                   

                    genome_fragments, comp = compleatness(record, 1-random.random()) 

                    if random.random() < 0.5:
                        cont_fragments, cont = contaminate(contamination.seq, record.seq, random.random()*0.3)
                    else:
                        cont_fragments = []
                        cont = 0
                    
                    if comp > 0.6: 

                        print(record.id, contamination.id)
                        all_fragments = genome_fragments + cont_fragments 
                        name = "bins_sparse/bin_"+str(n)+'.fasta'
                        f_bin = open(name, 'w')
                        SeqIO.write(all_fragments, f_bin, 'fasta')
                        f_bin.close()

                        print("%i,%5.4f,%5.4f" % (n,comp,cont), file=log)
                        n += 1
    log.close()
