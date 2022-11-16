#! /usr/bin/env python

import argparse
import pysam
from itertools import groupby
from collections import defaultdict
from operator import itemgetter


class GenomicTypeDB:
    TYPE_NAME = (
        'CDS',
        'intron',
        '5p_UTR',
        '3p_UTR',
        'start_codon',
        'stop_codon',
        'noncoding',
        'intergenic'
    )
    FEATURES_ORDER = (
        'gene',
        'transcript',
        'Selenocysteine',
        'exon',
        'CDS',
        'five_prime_utr',
        'three_prime_utr',
        'stop_codon',
        'start_codon'
    )

    def __init__(self, gtf):
        self._gtf = gtf
        self._tbx = pysam.TabixFile(gtf)

    def _fetch_raw_data(self, chr_, pos, strand=None):
        sites_data = self._tbx.fetch(chr_, pos - 1, pos, parser=pysam.asGTF())

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
    def _uniq_features(group_data):
        features = sorted(set(map(lambda site: site.feature, group_data)))
        return features

    @staticmethod
    def _is_intron(features):
        return ('transcript' in features) and (len(features) == 1)

    @staticmethod
    def _is_CDS(features):
        return ('transcript' in features) and ('CDS' in features)

    @staticmethod
    def _is_start_codon(features):
        return ('transcript' in features) and ('start_codon' in features)

    @staticmethod
    def _is_stop_codon(features):
        return ('transcript' in features) and ('stop_codon' in features)

    def _get_region_features(self, transcript, region_start, region_end):
        region_data = self._tbx.fetch(chr_, region_start, region_end, parser=pysam.asGTF())
        region_data = filter(lambda region: region.feature != 'gene', region_data)
        region_data = filter(lambda region: region.transcript_id == transcript.transcript_id, region_data)
        features = self._uniq_features(region_data)

        return features

    def _is_5p_UTR(self, features, data):
        if 'transcript' in features:
            if 'five_prime_utr' in features:
                return True
            elif 'UTR' in features:
                transcript = next(filter(lambda region: region.feature == 'transcript', data))
                exon = next(filter(lambda region: region.feature == 'exon', data))
                UTR = next(filter(lambda region: region.feature == 'UTR', data))

                strand = transcript.strand

                if strand == '+':
                    if UTR.end < exon.end:
                        return True
                    else:
                        if UTR.start == exon.start:
                            if UTR.start == transcript.start:
                                return True
                            else:
                                if UTR.end < transcript.end:
                                    LHS_region_len = UTR.start - transcript.start + 1
                                    RHS_region_len = transcript.end - UTR.end + 1
                                    if LHS_region_len <= RHS_region_len:

                                        LHS_features = self._get_region_features(transcript, transcript.start, UTR.start)

                                        if 'CDS' not in LHS_features:
                                            return True

                                    else:

                                        RHS_features = self._get_region_features(transcript, UTR.end, transcript.end)

                                        if 'CDS' in RHS_features:
                                            return True

                elif strand == '-':
                    if UTR.start > exon.start:
                        return True
                    else:
                        if UTR.end == exon.end:
                            if UTR.end == transcript.end:
                                return True
                            else:
                                if UTR.start > transcript.start:
                                    LHS_region_len = UTR.start - transcript.start + 1
                                    RHS_region_len = transcript.end - UTR.end + 1
                                    if LHS_region_len <= RHS_region_len:

                                        LHS_features = self._get_region_features(transcript, transcript.start, UTR.start)

                                        if 'CDS' in LHS_features:
                                            return True

                                    else:

                                        RHS_features = self._get_region_features(transcript, UTR.end, transcript.end)

                                        if 'CDS' not in RHS_features:
                                            return True

        return False

    def _is_3p_UTR(self, features, data):
        if 'transcript' in features:
            if 'three_prime_utr' in features:
                return True
            elif 'UTR' in features:
                transcript = next(filter(lambda region: region.feature == 'transcript', data))
                exon = next(filter(lambda region: region.feature == 'exon', data))
                UTR = next(filter(lambda region: region.feature == 'UTR', data))

                strand = transcript.strand

                if strand == '+':
                    if UTR.start > exon.start:
                        return True
                    else:
                        if UTR.end == exon.end:
                            if UTR.end == transcript.end:
                                return True
                            else:
                                if UTR.start > transcript.start:
                                    LHS_region_len = UTR.start - transcript.start + 1
                                    RHS_region_len = transcript.end - UTR.end + 1
                                    if LHS_region_len <= RHS_region_len:

                                        LHS_features = self._get_region_features(transcript, transcript.start, UTR.start)

                                        if 'CDS' in LHS_features:
                                            return True

                                    else:

                                        RHS_features = self._get_region_features(transcript, UTR.end, transcript.end)

                                        if 'CDS' not in RHS_features:
                                            return True

                elif strand == '-':
                    if UTR.end < exon.end:
                        return True
                    else:
                        if UTR.start == exon.start:
                            if UTR.start == transcript.start:
                                return True
                            else:
                                if UTR.end < transcript.end:
                                    LHS_region_len = UTR.start - transcript.start + 1
                                    RHS_region_len = transcript.end - UTR.end + 1
                                    if LHS_region_len <= RHS_region_len:

                                        LHS_features = self._get_region_features(transcript, transcript.start, UTR.start)

                                        if 'CDS' not in LHS_features:
                                            return True

                                    else:

                                        RHS_features = self._get_region_features(transcript, UTR.end, transcript.end)

                                        if 'CDS' in RHS_features:
                                            return True

        return False

    @staticmethod
    def _is_noncoding(features):
        return ('transcript' in features) and ('exon' in features) and (len(features) == 2)

    @classmethod
    def _to_gt_binary_array(cls, gt_list):
        gt_binary_dict = defaultdict(int)
        for gt in gt_list:
            gt_binary_dict[gt] = 1

        gt_binary_array = list(
            map(
                lambda gt: gt_binary_dict[gt],
                cls.TYPE_NAME
            )
        )

        return gt_binary_array

    def fetch(self, chr_, pos, strand=None,
              return_binary_array=False, show_isoform=False):
        genomic_types_data = []

        sites_data = list(self._fetch_raw_data(chr_, pos, strand))

        if not sites_data:
            genomic_types_data.append(('', 'intergenic'))

        else:
            for tid, gp in self._get_transcript_groups(sites_data):
                gp = list(gp)
                features = self._uniq_features(gp)

                if self._is_intron(features):
                    genomic_types_data.append((tid, 'intron'))

                if self._is_CDS(features):
                    genomic_types_data.append((tid, 'CDS'))

                if self._is_start_codon(features):
                    genomic_types_data.append((tid, 'start_codon'))

                if self._is_stop_codon(features):
                    genomic_types_data.append((tid, 'stop_codon'))

                if self._is_5p_UTR(features, gp):
                    genomic_types_data.append((tid, '5p_UTR'))

                if self._is_3p_UTR(features, gp):
                    genomic_types_data.append((tid, '3p_UTR'))

                if self._is_noncoding(features):
                    genomic_types_data.append((tid, 'noncoding'))

        if show_isoform:
            all_genomic_types = []
            genomic_types_data = sorted(genomic_types_data, key=itemgetter(0))
            for tid, gp in groupby(genomic_types_data, key=itemgetter(0)):
                genomic_types = sorted(set(map(itemgetter(1), gp)))

                if return_binary_array:
                    genomic_types = self._to_gt_binary_array(genomic_types)

                all_genomic_types.append([tid, genomic_types])

            return all_genomic_types

        else:
            genomic_types = sorted(set(map(itemgetter(1), genomic_types_data)))

            if return_binary_array:
                genomic_types = self._to_gt_binary_array(genomic_types)

            return genomic_types

    def _debug_show(self, chr_, pos, strand=None):
        sites_data = list(self._fetch_raw_data(chr_, pos, strand))

        for tid, gp in self._get_transcript_groups(sites_data):
            print(tid)
            for site in gp:
                print(site)
            print()

        for tid, gp in self._get_transcript_groups(sites_data):
            print(tid, list(self._uniq_features(gp)))


def create_parser():
    parser = argparse.ArgumentParser(
        usage=(
            "\n       %(prog)s -g [GTF_FILE] [SITES_FILE] > [OUT_FILE]\n"
            "  or\n"
            "       cat [SITE_FILE] | %(prog)s -g [GTF_FILE] - > [OUT_FILE]\n\n"
            "  The data format of sites_file is consist of three columns: chr_name, pos(1-base), and strand\n"
            "  eg. \n"
            "      1   633567  -\n"
            "      1   634095  -\n"
            "      1   732017  -\n"
            "      1   733251  +\n"
            "      1   746695  -"
        )
    )
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
    parser.add_argument(
        '-b',
        '--return_binary_array',
        action='store_true',
        help="Show results with binary format."
    )
    parser.add_argument(
        '--show-isoform',
        action='store_true',
        help="Show results in isoform level."
    )

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    gdb = GenomicTypeDB(args.gtf_file)

    titles = ['chr', 'pos', 'strand']

    if args.return_binary_array:
        titles.extend(list(GenomicTypeDB.TYPE_NAME))
    else:
        titles.append('genomic_types')

    if args.show_isoform:
        titles.append('transcript_id')

    print(*titles, sep='\t')

    try:
        for line in args.sites_file:
            chr_, pos, strand = line.rstrip('\n').split('\t')
            pos = int(pos)

            gt_res = gdb.fetch(
                chr_,
                pos,
                strand,
                return_binary_array=args.return_binary_array,
                show_isoform=args.show_isoform
            )

            if args.return_binary_array:
                if args.show_isoform:
                    for tid, gt in gt_res:
                        print(chr_, pos, strand, *gt, tid, sep='\t')
                else:
                    print(chr_, pos, strand, *gt_res, sep='\t')
            else:
                if args.show_isoform:
                    for tid, gt in gt_res:
                        print(chr_, pos, strand, '|'.join(gt), tid, sep='\t')
                else:
                    print(chr_, pos, strand, '|'.join(gt_res), sep='\t')

    except BrokenPipeError:
        pass
