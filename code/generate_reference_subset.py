#!/usr/bin/env python3

import argparse
import sys
from pyfaidx import Fasta
from random import choice
from random import seed
from random import randint

if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]

    parser = argparse.ArgumentParser(description="Create a randomised subset of a genome.")
    parser.add_argument('-g', '-genome', help='genome file (.fasta, can be gzipped, but must be indexed)')
    parser.add_argument('-a', '-assignments', help="A teb-delimited table with chromosomal assignmensts (scf, chr columns are expected)")
    parser.add_argument('-c', '-chromosome', help="Chromosome to be generated (must be included in -assignments table, defalt: all)", default = 'all')
    parser.add_argument('-n', '-number_of_chromosomes', help="Number of chromosomes (of length -l) to be generated (defalt: 1)", default = 1, type = int)
    parser.add_argument('-l', '-length', help='The total length of newly sampled reference (in Mbp, default: 20)', default = 1, type = int)
    parser.add_argument('-o', '-output', help='output pattern (default: sampled_genome)', default = 'sampled_genome.fasta')
    parser.add_argument('-s', '-seed', help='seed for generating random numbers (default: defined by time)', default = None, type = int)
    args = parser.parse_args(args)

    sys.stderr.write("processing {} assignment file\n".format(args.a))

    scf2asn = dict()
    does_the_chromosome_exist = False

    ##### process assignment table
    # args.a = 'tables/chr_assignments_Afus1.tsv'
    with open(args.a, 'r') as asgn:
        header = asgn.readline().rstrip('\n').split('\t')
        if not ('scf' in header and 'chr' in header):
            sys.stderr.write("Error: The assignment file does not have scf and chr columns in the header\n".format(args.a))
            exit(-1)
        scf_i = [i for i, colname in enumerate(header) if colname == 'scf'][0]
        chr_i = [i for i, colname in enumerate(header) if colname == 'chr'][0]
        for line in asgn:
            scf_info = line.rstrip('\n').split('\t')
            scf = scf_info[scf_i]
            chr = scf_info[chr_i]
            if chr == args.c or args.c == 'all':
                scf2asn[scf] = chr
                does_the_chromosome_exist = True
        sys.stderr.write('loaded {} scaffolds with feature: {}\n'.format(len(scf2asn), args.a))

    if not does_the_chromosome_exist:
        sys.stderr.write('Error: It apprears that {} chromosome assimeng is not in your {} assignment file.\n'.format(args.c, args.a))
        exit(-1)

    ##### Setting up the random number generator
    if args.s == None:
        used_seed = randint(100000, 99999999)
    else:
        used_seed = args.s

    seed(used_seed)
    sys.stderr.write('This sampling is can be regenerated using seed {} (parameter -s)\n'.format(used_seed))

    ##### Now we load the reference genome index and generate the sub sampled genome
    genome = Fasta(args.g)
    sys.stderr.write('loaded {} genome file\n'.format(args.g))

    chromosome = 0
    total_desired = int(args.l * 1e6)
    added_scaffolds = set()

    with open(args.o, 'w') as output_fasta:
        while chromosome < args.n:
            chromosome += 1
            output_fasta.write('>simulated_' + args.c + '_chromosome_' + str(chromosome) + '_with_' + str(used_seed) + '_seed\n')

            sampled_length = 0
            while sampled_length < total_desired:
                picked_scf = choice(list(scf2asn.keys()))
                # if the scaffold was not picked yet AND if it has the desired chromosome assignment
                if not picked_scf in added_scaffolds and (scf2asn[picked_scf] == args.c or args.c == 'all'):
                    added_scaffolds.add(picked_scf)
                    scaffold_record = genome[picked_scf]
                    if len(scaffold_record) + sampled_length > total_desired:
                        seq_to_print = scaffold_record[0:(total_desired - sampled_length)]
                    else:
                        seq_to_print = scaffold_record[0:]
                    sampled_length += len(scaffold_record)

                    output_fasta.write(str(seq_to_print))
            output_fasta.write('\n')

    sys.stderr.write('Done\n')
