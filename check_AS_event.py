#! /usr/bin/env python

import argparse
import pysam
import re
from itertools import groupby
from collections import namedtuple

import logging


class GTFDB:
    def __init__(self, gtf):
        self._gtf = gtf
        self._tbx = pysam.TabixFile(gtf)

    def _fetch_tabix(self, chr_, start, end):
        return self._tbx.fetch(chr_, start, end, parser=pysam.asGTF())

    def fetch_site(self, chr_, pos, strand=None):
        sites_data = self._fetch_tabix(chr_, pos - 1, pos)

        if strand is not None:
            sites_data = filter(lambda site: site.strand == strand, sites_data)

        return sites_data

    @staticmethod
    def group_by_gene(sites_data):
        sites_data = sorted(sites_data, key=lambda site: site.gene_id)
        for gid, gp in groupby(sites_data, key=lambda site: site.gene_id):
            yield gid, list(gp)

    @staticmethod
    def group_by_transcript(sites_data):
        sites_data = filter(lambda site: site.feature != 'gene', sites_data)
        sites_data = sorted(sites_data, key=lambda site: site.transcript_id)
        for tid, gp in groupby(sites_data, key=lambda site: site.transcript_id):
            yield tid, list(gp)

    @staticmethod
    def get_feature_data(feature_name, gp):
        data = list(filter(
            lambda region: region.feature == feature_name,
            gp
        ))
        return data

    @staticmethod
    def get_attribute(region, attr):
        m = re.search(rf'{attr} \"?([^\";]+)\"?;', region.attributes)
        if m:
            return m.group(1)

    def get_exons_of_transcript(self, transcript):
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


class SplicingSitesDB(GTFDB):
    _site = namedtuple('Site', ['chr_', 'pos', 'strand', 'attrs'])

    def _get_additional_attrs(self, exon):
        attrs = {
            'gene_id': exon.gene_id,
            'transcript_id': exon.transcript_id,
            'exon_start_end': (exon.start + 1, exon.end),
            'exon_number': int(self.get_attribute(exon, 'exon_number')),
        }
        return attrs

    def fetch(self, chr_, pos, strand=None):
        splicing_sites = []

        sites_data = list(self.fetch_site(chr_, pos, strand))

        if sites_data:
            logging.debug('There are "sites_data"!')
            for gid, g_gp in self.group_by_gene(sites_data):
                for tid, t_gp in self.group_by_transcript(g_gp):
                    transcript = self.get_feature_data('transcript', t_gp)[0]

                    if (pos != transcript.start + 1) and (pos != transcript.end):

                        exon = self.get_feature_data('exon', t_gp)

                        if exon:
                            logging.debug('Exon!')
                            exon = exon[0]
                            logging.debug((exon.start, exon.end))
                            logging.debug(exon)

                            if (pos == exon.start + 1) or (pos == exon.end):
                                logging.debug('Splicing site!')
                                attrs = self._get_additional_attrs(exon)
                                site = self._site(chr_, pos, strand, attrs)
                                logging.debug(site)

                                splicing_sites.append(site)
                        else:
                            # intron
                            logging.debug("Intron!")
                            pass
        else:
            # intergenic
            logging.debug('Intergenic!')
            pass

        return splicing_sites, sites_data


class ASEventChecker:
    RESULT_TITLES = (
        'is_annotated_splicing_site',
        'has_AS_event'
    )

    def __init__(self, gtf):
        self._ssdb = SplicingSitesDB(gtf)

    @staticmethod
    def _get_gene_group(sites_data, gene_id):
        gene_gp = list(
            filter(
                lambda site: site.gene_id == gene_id,
                sites_data
            )
        )
        return gene_gp

    def check(self, chr_, pos, strand=None):
        logging.info(f'Checking the site: ({chr_}, {pos}, {strand})')

        splicing_sites, sites_data = self._ssdb.fetch(chr_, pos, strand)

        logging.debug('==============================')
        logging.debug(splicing_sites)
        for site in sites_data:
            logging.debug(site)

        if splicing_sites:

            for splicing_site in splicing_sites:

                # Case 1: genomic type changed from exon to intron
                # Possible AS event:
                #  - exon skipping,
                #  - mutually exclusive exons,
                #  - alternative 5' donor site,
                #  - alternative 3' donor site
                gene_gp = self._get_gene_group(
                    sites_data,
                    splicing_site.attrs['gene_id']
                )
                for tid, t_gp in GTFDB.group_by_transcript(gene_gp):
                    t_gp = list(t_gp)
                    if len(t_gp) == 1:
                        logging.info('Case 1!')
                        return (1, 1)
                
                # Case 2: genomic type changed from exon to exon
                # Possible AS event:
                #  - alternative 5' donor site,
                #  - alternative 3' donor site,
                #  - intron retention
                if splicing_site.pos == splicing_site.attrs['exon_start_end'][0]:
                    neighbor_site_pos = splicing_site.pos - 1
                elif splicing_site.pos == splicing_site.attrs['exon_start_end'][1]:
                    neighbor_site_pos = splicing_site.pos + 1

                neighbor_sites_data = list(
                    self._ssdb.fetch_site(
                        splicing_site.chr_,
                        neighbor_site_pos,
                        splicing_site.strand
                    )
                )
                nb_gene_gp = self._get_gene_group(
                    neighbor_sites_data,
                    splicing_site.attrs['gene_id']
                )
                nb_site_exons = self._ssdb.get_feature_data('exon', nb_gene_gp)

                if len(nb_site_exons) > 0:
                    logging.info('Case 2!')
                    return (1, 1)
            else:
                logging.info('There are no AS events on this site!')
                return (1, 0)
        else:
            logging.info('The site is not an annotated splicing site!')
            return (0, 0)


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
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--log_file', default='check_AS_event.log')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    if args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        filename=args.log_file,
        encoding='utf-8',
        level=level
    )

    AS_checker = ASEventChecker(args.gtf_file)

    titles = ['chr', 'pos', 'strand'] + list(ASEventChecker.RESULT_TITLES)
    print(*titles, sep='\t')

    try:
        for line in args.sites_file:
            chr_, pos, strand = line.rstrip('\n').split('\t')
            pos = int(pos)

            AS_check_results = AS_checker.check(chr_, pos, strand)

            print(chr_, pos, strand, *AS_check_results, sep='\t')

    except BrokenPipeError:
        pass
