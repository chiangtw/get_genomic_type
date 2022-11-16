#! /usr/bin/env python

import argparse
import pysam
import re
from itertools import groupby, chain
from collections import namedtuple
from bisect import bisect


class AnnotationDB:
    _anno_data = namedtuple(
        'AnnoData',
        (
            'gene_id',
            'transcript_id',
            'exon_number',
            'num_of_exons',
            'len_of_transcript',
            'on_exon_boundary'
        )
    )

    def __init__(self, gtf):
        self._gtf = gtf
        self._tbx = pysam.TabixFile(gtf)

    def _fetch_tabix(self, chr_, start, end):
        return self._tbx.fetch(chr_, start, end, parser=pysam.asGTF())

    def _fetch_raw_data(self, chr_, pos, strand=None):
        sites_data = self._fetch_tabix(chr_, pos - 1, pos)

        if strand is not None:
            sites_data = filter(lambda site: site.strand == strand, sites_data)

        return sites_data

    @staticmethod
    def _get_transcript_groups(sites_data):
        sites_data = filter(lambda site: site.feature != 'gene', sites_data)
        sites_data = sorted(sites_data, key=lambda site: site.transcript_id)
        for tid, gp in groupby(sites_data, key=lambda site: site.transcript_id):
            yield tid, list(gp)

    @staticmethod
    def _get_feature_data(feature_name, gp):
        data = list(filter(
            lambda region: region.feature == feature_name,
            gp
        ))

        if not data:
            return None
        else:
            return data[0]

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
    def _get_attribute(region, attr):
        m = re.search(rf'{attr} \"?([^";]+)\"?;', region.attributes)
        if m:
            return m.group(1)

    @staticmethod
    def _get_psuedo_exon_num_for_intron(exons, pos):
        exons = sorted(exons, key=lambda e: e.start)
        exons_pos_list = list(
            chain.from_iterable(
                [[e.start, e.end] for e in exons]
            )
        )
        exon_num = bisect(exons_pos_list, pos) / 2 + 0.5

        return str(exon_num)

    def fetch(self, chr_, pos, strand=None):
        all_anno_data = []

        sites_data = list(self._fetch_raw_data(chr_, pos, strand))

        if not sites_data:
            # intergenic case
            all_anno_data.append(self._anno_data('', '', '', '', '', 0))
        else:
            for tid, gp in self._get_transcript_groups(sites_data):
                exon = self._get_feature_data('exon', gp)
                transcript = self._get_feature_data('transcript', gp)
                transcript_exons = self._get_exons_of_transcript(transcript)

                gid = transcript.gene_id
                num_of_exons = len(transcript_exons)
                len_of_transcript = self._sum_of_exons(transcript_exons)
                on_exon_boundary = 0

                if exon:
                    exon_number = self._get_attribute(exon, 'exon_number')

                    if (pos == exon.start + 1) or (pos == exon.end):
                        on_exon_boundary = 1
                else:
                    exon_number = self._get_psuedo_exon_num_for_intron(
                        transcript_exons,
                        pos
                    )

                all_anno_data.append(
                    self._anno_data(
                        gid,
                        tid,
                        exon_number,
                        num_of_exons,
                        len_of_transcript,
                        on_exon_boundary
                    )
                )

        return all_anno_data


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-g',
        '--gtf_file',
        help=(
            "The input gtf file should be zipped by \"bgzip\","
            " and there should be a tabix index file(.tbi) in the same directory."
        )
    )
    parser.add_argument(
        'sites_file',
        type=argparse.FileType('r'),
        help="To receive data from stdin, use '-' to this field."
    )

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    anno_db = AnnotationDB(args.gtf_file)

    titles = ['chr', 'pos', 'strand'] + list(anno_db._anno_data._fields)
    print(*titles, sep='\t')

    try:
        for line in args.sites_file:
            chr_, pos, strand = line.rstrip('\n').split('\t')
            pos = int(pos)

            anno_data = anno_db.fetch(chr_, pos, strand)
            for data in anno_data:
                print(chr_, pos, strand, *data, sep='\t')

    except BrokenPipeError:
        pass
