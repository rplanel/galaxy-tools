#!/usr/bin/env python3
# coding: utf-8

import locale
import sys
import os
import argparse
from itertools import chain, islice, tee, count
import logging

# BioPython
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s : %(levelname)s : %(message)s"
)
logger = logging.getLogger()


def main():
    args = parse_arguments()
    # extract the basename and the file extension
    basename, file_extension = os.path.splitext(args.sequences)
    # extract the filename (with no extension)
    _, filename = os.path.split(basename)
    chunk_size = args.chunk_size
    # split the sequences in chunks
    if chunk_size:
        logger.info('%s = %s', "chunk size asked", chunk_size)
        sequences_record = gen_sequence_record(args.sequences, args.format)
        chunks = gen_get_chunks_by_size(sequences_record, chunk_size)
    else:
        logger.info('%s = %s', "number of chunks asked", args.nb_chunk)
        chunks = gen_get_chunks(args.sequences, args.format, args.nb_chunk)

    # Write the chunks in numbered files.
    write_chunks(chunks, args.output, filename, file_extension, args.format)


def gen_get_chunks(sequences_path, sequences_format, nb_chunk):

    # First record to count the sequences
    sequences_record_to_count = gen_sequence_record(sequences_path, sequences_format)
    # Get the number of sequences
    nb_sequences = get_nb_sequences(sequences_record_to_count)
    logger.info('%s = %i', "Number of sequences", nb_sequences)
    # Second record to that will be splitted
    sequences_to_split = gen_sequence_record(sequences_path, sequences_format)

    # Get the size of the chunks
    chunk_size = int(nb_sequences / nb_chunk) if nb_sequences > nb_chunk else 1
    return gen_get_chunks_by_size(sequences_to_split, chunk_size)


def gen_get_chunks_by_size(iterable, size=10):
    logger.info('%s = %i', "chunk size got", size)
    iterator = iter(iterable)
    for first in iterator:
        yield chain([first], islice(iterator, size - 1))


def gen_sequence_record(sequences_path, sequence_format):
    return SeqIO.parse(sequences_path, sequence_format)


def get_nb_sequences(sequences):
    return sum(1 for _ in sequences)


def write_chunks(iterable, dirname, filename, file_extension, sequence_format):
    for idx, chunk in enumerate(iterable):
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        output_file = os.path.join(
            dirname, filename + "-chunk-" + str(idx + 1) + file_extension
        )
        with open(output_file, mode="w") as output_handle:
            count_seq, seq_to_write = tee(chunk, 2)
            logger.info(
                "%s : number of seuquences = %i", output_file, len(list(count_seq))
            )
            SeqIO.write(seq_to_write, output_handle, sequence_format)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Split fasta/fastq files")
    parser.add_argument(
        "-s", "--sequences", type=str, help="File that contains the sequences"
    )

    parser.add_argument("-f", "--format", type=str, help="File format (fastq, fasta)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-c", "--chunk-size", type=int, help="The number of sequences by chunks."
    )
    group.add_argument("-n", "--nb-chunk", help="Number of chunks", type=int)

    parser.add_argument("-o", "--output", type=str, default="./", help="The output directory where the chunks will be saved")

    return parser.parse_args()


if __name__ == "__main__":
    logger.info("START")
    main()
    logger.info("FINISHED")
