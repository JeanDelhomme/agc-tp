#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):

    """
    Opens a fasta file and reads its sequences.

    Parameters
    ----------
        amplicon_file : str
            The name of the fastq file.
    Returns
    -------
        yield
        A yield containing the sequences.
    """

    with gzip.open(amplicon_file, "rt") as  file:
        seq = ""
        
        for line in file:
            if not line.startswith(">"):
                seq += line.strip()
            if line.startswith(">"):
                if len(seq) >= minseqlen:
		            
                    yield seq
                seq = ""

        if len(seq) >= minseqlen:
    	    
            yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):

    """
    Returns a generator of unique sequences and their occurences.

    Parameters
    ----------
        amplicon_file : str
            The name of the fasta file.

        minseqlen : int
            The minimum length of the sequences we want to keep.

        mincount :
            The minimum count of the sequences we want to keep.

    Returns
    -------
        yield
        A yield containing the sequences and their count.
    """

    seq_list = list(read_fasta(amplicon_file, minseqlen))
    dico = Counter(seq_list)

    sorted_dico = sorted(dico)

    for key, value in sorted(dico.items(), key = lambda item: item[1], reverse = True):
        if value >= mincount:
            yield [key, value]


def get_identity(alignment_list):
    
    """
    Calculates an identity score between two aligned sequences.

    Parameters
    ----------
        amplicon_list :
            A list of two aligned sequences.
    Returns
    -------
        id : int
            The identity score in %.
    """

    identical_nuc = 0
    gap_nb = 0
    seq1 = alignment_list[0]
    seq2 = alignment_list[1]

    for AA in range(len(seq1)):
        if seq1[AA] == "-":
            gap_nb += 1

        if seq1[AA] != "-" or seq2[AA] != "-":

            if seq1[AA] == seq2[AA]:
                identical_nuc += 1

    id = identical_nuc/len(seq1)*100

    return id

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    
    """
    Select OTUs from a fasta file containing sequences.

    Parameters
    ----------
        amplicon_file : str
            The name of the fasta file.

        minseqlen : int
            The minimum length of the sequences we want to keep.

        mincount :
            The minimum count of the sequences we want to keep.

        chunk_size :
            NA

        kmer_sizer :
            NA

    Returns
    -------
        OTU_list :
            A list containing tuples of OTUs and their counts.
    """

    candidate_list = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    OTU_list = [candidate_list[0]]

    for seq1 in candidate_list[1:]:

        id = 0

        for seq2 in OTU_list:

            alignment_list = nw.global_align(seq1[0], seq2[0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            id = max(id, get_identity(alignment_list))

        if  id <= 97:
            OTU_list.append(seq1)

    return OTU_list


def write_OTU(OTU_list, output_file):

    """
    Save an OTU list in a file.

    Parameters
    ----------
        OTU_list :
        A list of tuples containing the contigs and their counts. 

        output_file :
        The name of the file where the results are saved.
    """

    with open(output_file, "w") as file:
        for i, (sequence, count) in enumerate(OTU_list):
            file.write(f">OTU_{i+1} occurrence:{count}\n{textwrap.fill(sequence, width = 80)}\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(OTU_list, args.output_file)

#==============================================================
# Chimera removal section
#==============================================================

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def get_chunks(sequence, chunk_size):
    """Split sequences in a least 4 chunks
    """
    pass

def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    pass

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    pass

def detect_chimera(perc_identity_matrix):
    pass

def search_mates(kmer_dict, sequence, kmer_size):
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


if __name__ == '__main__':
    main()
