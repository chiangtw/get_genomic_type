#! /usr/bin/env python

import argparse
import pysam
import re
from itertools import chain, product
from collections import namedtuple
from bisect import bisect_left

import logging

logger = logging.getLogger(__name__)


class IntronDB:
    _anno_data = namedtuple(
        'AnnoData',
        (
            'gene_id',
            'transcript_id',
            'intron_num',
            'intron_start',
            'intron_end',
            'num_of_exons',
            'len_of_transcript'
        )
    )

    def __init__(self, gtf):
        self._gtf = gtf
        self._tbx = pysam.TabixFile(gtf)

        self.empty_output = self._generate_empty_output()

    def _fetch_tabix(self, chr_, start, end):
        return self._tbx.fetch(chr_, start, end, parser=pysam.asGTF())

    def _fetch_raw_data(self, chr_, pos, strand=None):
        sites_data = self._fetch_tabix(chr_, pos - 1, pos)

        if strand is not None:
            sites_data = filter(lambda site: site.strand == strand, sites_data)

        return sites_data

    @staticmethod
    def _get_transcripts(sites_data):
        transcripts_data = filter(lambda site: site.feature == 'transcript', sites_data)
        transcripts_data = sorted(transcripts_data, key=lambda site: site.transcript_id)
        for transcript in transcripts_data:
            yield transcript

    def _get_exons_of_transcript(self, transcript):
        region_data = self._fetch_tabix(
            transcript.contig,
            transcript.start,
            transcript.end
        )

        exons = list(filter(
            lambda region: (region.feature == 'exon') and (region.transcript_id == transcript.transcript_id),
            region_data
        ))

        return exons

    @staticmethod
    def _sum_of_exons(exons):
        return sum(map(lambda e: e.end - e.start + 1, exons))

    @staticmethod
    def _get_intron_loc(exons, pos):
        exons_pos_list = list(
            chain.from_iterable(
                [[e.start, e.end] for e in exons]
            )
        )

        pos_idx = bisect_left(exons_pos_list, pos)

        if (pos_idx > 0) and (pos_idx < len(exons) * 2):
            q, r = divmod(pos_idx, 2)

            if r == 0:
                intron_idx = (q - 1, q)

                return intron_idx

    @staticmethod
    def _get_exon_num(exon):
        m = re.search(r'exon_number \"?([0-9]+)\"?;', exon.attributes)
        if m:
            return m[1]

    def fetch(self, chr_, pos, strand=None):
        sites_data = list(self._fetch_raw_data(chr_, pos, strand))

        if not sites_data:
            # intergenic case
            pass

        else:
            for transcript in self._get_transcripts(sites_data):
                exons = self._get_exons_of_transcript(transcript)
                exons = sorted(exons, key=lambda e: e.start)

                gid = transcript.gene_id
                num_of_exons = len(exons)
                len_of_transcript = self._sum_of_exons(exons)

                intron_idx = self._get_intron_loc(exons, pos)

                if intron_idx:
                    intron_up_down_exons = exons[intron_idx[0]:(intron_idx[1] + 1)]

                    intron_num = ','.join(
                        map(self._get_exon_num, intron_up_down_exons)
                    )
                    intron_start = intron_up_down_exons[0].end + 1
                    intron_end = (intron_up_down_exons[1].start + 1) - 1

                    yield self._anno_data(
                        gid,
                        transcript.transcript_id,
                        intron_num,
                        intron_start,
                        intron_end,
                        num_of_exons,
                        len_of_transcript
                    )

                else:
                    # pos in exon
                    pass

    @classmethod
    def _generate_empty_output(cls):
        return cls._anno_data(*([''] * len(cls._anno_data._fields)))


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-g',
        '--gtf_file',
        help=(
            "The input gtf file should be zipped by \"bgzip\","
            " and there should be a tabix index file(.tbi) in the same directory."
        ),
        required=True
    )
    parser.add_argument(
        'in_file',
        type=argparse.FileType('r'),
        help="To receive data from stdin, use '-' to this field."
    )
    parser.add_argument(
        '--sites',
        dest='is_sites',
        action='store_true',
        help="Use this flag if the input are sites."
    )

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    intron_db = IntronDB(args.gtf_file)

    try:
        title_printed = False
        if args.is_sites:
            for line in args.in_file:

                line_data = line.rstrip('\n').split('\t')
                chr_, pos, strand = line_data[:3]
                pos = int(pos)

                if not title_printed:
                    title = ['chr', 'pos', 'strand']
                    title += [f'c{i}' for i, _ in enumerate(line_data[3:], start=4)]
                    title += list(intron_db._anno_data._fields)
                    print(*title, sep='\t')
                    title_printed = True

                intron_data = list(intron_db.fetch(chr_, pos, strand))

                if len(intron_data) > 0:
                    for data in intron_data:
                        print(*line_data, *data, sep='\t')
                else:
                    print(*line_data, *intron_db.empty_output, sep='\t')
        else:
            for line in args.in_file:
                line_data = line.rstrip('\n').split('\t')
                chr_, pos1, pos2, strand = line_data[:4]
                pos1, pos2 = list(map(int, [pos1, pos2]))

                if not title_printed:
                    title = ['chr', 'pos1', 'pos2', 'strand']
                    title += [f'c{i}' for i, _ in enumerate(line_data[4:], start=4)]

                    raw_titles = ['intron_pos'] + list(intron_db._anno_data._fields)
                    title += list(map(lambda t: t + '_1', raw_titles))
                    title += list(map(lambda t: t + '_2', raw_titles))

                    print(*title, sep='\t')
                    title_printed = True

                if pos1 > pos2:
                    intron_pos_1 = pos1 + 1
                    intron_pos_2 = pos2 - 1
                else:
                    intron_pos_1 = pos1 - 1
                    intron_pos_2 = pos2 + 1

                intron_data_1 = list(intron_db.fetch(chr_, intron_pos_1, strand))
                intron_data_2 = list(intron_db.fetch(chr_, intron_pos_2, strand))

                if len(intron_data_1) > 0:
                    if len(intron_data_2) > 0:
                        for data1, data2, in product(intron_data_1, intron_data_2):
                            print(*line_data, intron_pos_1, *data1, intron_pos_2, *data2, sep='\t')
                    else:
                        for data in intron_data_1:
                            print(*line_data, intron_pos_1, *data, intron_pos_2, *intron_db.empty_output, sep='\t')
                else:
                    if len(intron_data_2) > 0:
                        for data in intron_data_2:
                            print(*line_data, intron_pos_1, *intron_db.empty_output, intron_pos_2, *data, sep='\t')
                    else:
                        print(
                            *line_data,
                            intron_pos_1,
                            *intron_db.empty_output,
                            intron_pos_2,
                            *intron_db.empty_output,
                            sep='\t'
                        )
    except BrokenPipeError:
        pass
