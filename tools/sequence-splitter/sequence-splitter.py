#!/usr/bin/env python2
# coding: utf-8

import argparse
import logging
import os
from itertools import chain, islice, tee

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
        logger.info(
            "%s = %s", "Number of sequences per chunk parameter", chunk_size)
        sequences_record = gen_sequence_record(args.sequences, args.format)
        chunks = gen_chunks_of_size(sequences_record, chunk_size)
    else:
        logger.info("%s = %s", "number of chunks parameter", args.nb_chunk)
        chunks = gen_chunks(args.sequences, args.format, args.nb_chunk)

    # Write the chunks in numbered files.
    write_chunks(chunks, args.output, filename, file_extension, args.format)


def gen_chunks(sequences_path, sequences_format, nb_chunk):
    """[summary]
    Split sequences to have a number of chunks defined by nb_chunks
    Arguments:
        sequences_path {[type]} -- Path to the sequence file
        sequences_format {[type]} -- Format the sequence
        nb_chunk {int} -- Number of chunks we want to obtain

    Returns:
        [type] -- [description]
    """
    # First record to count the sequences
    sequences_record_to_count = gen_sequence_record(
        sequences_path, sequences_format)
    # Get the number of sequences
    nb_sequences = get_nb_sequences(sequences_record_to_count)
    logger.info("%s = %i", "Number of sequences total", nb_sequences)
    # Second record to that will be splitted
    sequences_to_split = gen_sequence_record(sequences_path, sequences_format)

    # Get the size of the chunks
    chunk_size = int(nb_sequences / nb_chunk) if nb_sequences > nb_chunk else 1
    return gen_chunks_of_size(
        sequences_to_split,
        chunk_size,
        nb_chunk,
        nb_sequences
    )


def gen_chunks_of_size(
    iterable, size=10,
    limit_nb_chunk=None,
    nb_sequences=None
):
    """[summary]
    Split sequences by size
    Arguments:
        iterable {[type]} -- [description]

    Keyword Arguments:
        size {int} -- [description] (default: {10})
    """
    iterator = iter(iterable)
    if limit_nb_chunk is not None and nb_sequences is not None:
        nb_chunk_left = limit_nb_chunk
        nb_sequences_left = nb_sequences
        for first in iterator:
            if (size + 1) * nb_chunk_left > nb_sequences_left:
                nb_sequences_left -= size
                nb_chunk_left -= 1
                yield chain([first], islice(iterator, size - 1))
            else:
                nb_chunk_left -= 1
                size += 1
                nb_sequences_left -= size
                yield chain([first], islice(iterator, size - 1))
    else:
        for first in iterator:
            yield chain([first], islice(iterator, size - 1))


def gen_sequence_record(sequences_path, sequence_format):
    return SeqIO.parse(sequences_path, sequence_format)


def get_nb_sequences(sequences):
    """[summary]
        Compute the number of sequences
    Arguments:
        sequences {[type]} -- Iterable of sequences

    Returns:
        [type] -- Number of sequences
    """
    return sum(1 for _ in sequences)


def write_chunks(iterable, dirname, filename, file_extension, sequence_format):
    sequence_total = 0
    for idx, chunk in enumerate(iterable):
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        output_file = os.path.join(
            dirname, filename + "-chunk-" + str(idx + 1) + file_extension
        )

        with open(output_file, mode="w") as output_handle:
            count_seq, seq_to_write = tee(chunk, 2)
            nb_seq = len(list(count_seq))
            logger.info("%s : number of sequences = %i", output_file, nb_seq)
            sequence_total += nb_seq
            SeqIO.write(seq_to_write, output_handle, sequence_format)
    logger.info("%s = %i", "Total number of chunks", idx + 1)
    logger.info("%s = %i", "Number of sequences total", sequence_total)


def positive_integer(str_value):
    """[summary]
    Define a type for argparse in order to enforce integer greater than 0
    Arguments:
        str_value {[type]} -- Value got by argparse

    Raises:
        argparse.ArgumentTypeError: When the value is not an integer > 0

    Returns:
        [type] -- [description]
    """
    value = int(str_value)
    if isinstance(value, int) and value > 0:
        return value
    else:
        msg = "%r is not an integer > 0" % value
        raise argparse.ArgumentTypeError(msg)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Split fasta/fastq files")
    parser.add_argument(
        "-s", "--sequences", type=str, help="File that contains the sequences"
    )

    parser.add_argument("-f", "--format", type=str,
                        help="File format (fastq, fasta)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-c", "--chunk-size", type=positive_integer,
        help="The number of sequences by chunks."
    )
    group.add_argument(
        "-n", "--nb-chunk",
        type=positive_integer,
        help="Number of chunks"
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="./",
        help="The output directory where the chunks will be saved",
    )

    return parser.parse_args()


if __name__ == "__main__":
    logger.info("START")
    main()
    logger.info("FINISHED")
