#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import random
import sys
import argparse
from traceback import print_exc
from collections import defaultdict
import gffutils
import numpy
import requests
import json
import sys
import time
from typing import List, Dict, Optional
import argparse
import logging
from dataclasses import dataclass


# Intron length distribution (internal v. not interal)
# terminal exon lengths
# introns w GTAG, etc.
# DS with 100% data, 90%, 80%, and so on



@dataclass
class ParalogInfo:
    """Data class to store paralog information."""
    target_id: str
    target_species: str
    target_protein_id: str
    source_id: str
    source_species: str
    source_protein_id: str
    percent_identity: float
    dnds_ratio: Optional[float]
    taxonomy_level: str
    type: str


class EnsemblParalogFinder:
    """Class to handle retrieval of paralogous genes from Ensembl."""

    def __init__(self):
        self.base_url = "https://rest.ensembl.org"
        self.headers = {
            "Content-Type": "application/json",
            "User-Agent": "ParalogFinder/1.0"
        }

        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)

    def get_paralogs(self, gene_id: str, species: str = "human") -> List[ParalogInfo]:
        """
        Retrieve paralogous genes for the given gene ID.

        Args:
            gene_id: The Ensembl gene ID
            species: Species name (default: human)

        Returns:
            List of ParalogInfo objects containing paralog information
        """
        url = f"{self.base_url}/homology/id/{species}/{gene_id}"
        params = {
            "content-type": "application/json",
            "type": "paralogues"
        }

        try:
            response = requests.get(url, headers=self.headers, params=params)

            if response.status_code != 200:
                self.logger.error(f"Failed to retrieve paralogs: {response.status_code}")
                return []

            data = response.json()
            return self._parse_paralog_response(data)

        except requests.exceptions.RequestException as e:
            self.logger.error(f"Error retrieving paralogs: {e}")
            return []
        except json.JSONDecodeError as e:
            self.logger.error(f"Error parsing response: {e}")
            return []

    def _parse_paralog_response(self, response_data: Dict) -> List[ParalogInfo]:
        """
        Parse the Ensembl API response and extract paralog information.

        Args:
            response_data: JSON response from Ensembl API

        Returns:
            List of ParalogInfo objects
        """
        paralogs = []

        try:
            if 'data' not in response_data or not response_data['data']:
                return paralogs

            for homology_data in response_data['data'][0].get('homologies', []):
                # Get target and source information
                target = homology_data.get('target', {})
                source = homology_data.get('source', {})

                paralog = ParalogInfo(
                    target_id=target.get('id', ''),
                    target_species=target.get('species', ''),
                    target_protein_id=target.get('protein_id', ''),
                    source_id=source.get('id', ''),
                    source_species=source.get('species', ''),
                    source_protein_id=source.get('protein_id', ''),
                    percent_identity=float(target.get('perc_id', 0)),
                    dnds_ratio=homology_data.get('dn_ds'),
                    taxonomy_level=homology_data.get('taxonomy_level', ''),
                    type=homology_data.get('type', '')
                )

                paralogs.append(paralog)

        except Exception as e:
            self.logger.error(f"Error parsing paralog data: {e}")

        return paralogs


def format_output(paralogs: List[ParalogInfo], format_type: str = 'text') -> str:
    """
    Format the paralog information for output.

    Args:
        paralogs: List of ParalogInfo objects
        format_type: Output format ('text' or 'json')

    Returns:
        Formatted string of paralog information
    """
    if format_type == 'json':
        # Convert to dictionary for JSON output
        paralog_data = [vars(p) for p in paralogs]
        return json.dumps(paralog_data, indent=2)

    # Text format
    output = []
    output.append(f"Found {len(paralogs)} paralogs:")
    output.append("=" * 80)

    for i, p in enumerate(paralogs, 1):
        output.append(f"Paralog {i}:")
        output.append("-" * 40)
        output.append(f"Type: {p.type}")
        output.append(f"Taxonomy Level: {p.taxonomy_level}")
        output.append("\nTarget Gene:")
        output.append(f"  ID: {p.target_id}")
        output.append(f"  Species: {p.target_species}")
        output.append(f"  Protein ID: {p.target_protein_id}")
        output.append("\nSource Gene:")
        output.append(f"  ID: {p.source_id}")
        output.append(f"  Species: {p.source_species}")
        output.append(f"  Protein ID: {p.source_protein_id}")
        output.append("\nStatistics:")
        output.append(f"  Percent Identity: {p.percent_identity:.2f}%")
        if p.dnds_ratio is not None:
            output.append(f"  dN/dS Ratio: {p.dnds_ratio}")
        output.append("=" * 80)

    return "\n".join(output)


def load_genes(inf):
    gene_list = []
    for l in open(inf):
        v = l.strip().split()
        gene_name = v[0]
        if v[0].find(".") != -1:
            continue
        gene_list.append(gene_name)

    return gene_list


def load_genedb(genedb):
    gffutils_db = gffutils.FeatureDB(genedb)
    gene_dict = {}
    gene_id2name = {}
    print("Loading genedb from %s" % genedb)
    for g in gffutils_db.features_of_type('gene'):
        if "gene_name" in g.attributes:
            gene_name = g.attributes["gene_name"][0]
        else:
            continue
        if "gene_type" in g.attributes:
            gene_type = g.attributes["gene_type"][0]
        else:
            continue
        g_id = g.id.split('.')[0]
        gene_id2name[g_id] = gene_name
        exon_counts = []
        tlens = []
        for t in gffutils_db.children(g, featuretype=("transcript", "mRNA")):
            ex = 0
            tlen = 0
            for e in gffutils_db.children(t, featuretype="exon"):
                ex += 1
                tlen += e.end - e.start + 1
            exon_counts.append(ex)
            tlens.append(tlen)
        gene_dict[gene_name] = (gene_type, len(exon_counts), numpy.mean(exon_counts), max(tlens), g_id)
    return gene_dict, gene_id2name


def count_stats(gene_list, gene_dict):
    total = 0
    unspliced = 0
    gene_lengths = []
    spliced_gene_lengths = []
    #pcg = []

    gene_types = defaultdict(int)

    for g in gene_list:
        if g not in gene_dict: continue
        total += 1
        gene_info = gene_dict[g]
        gene_types[gene_info[0]] += 1
        #if gene_info[0] == 'protein_coding':
        #    pcg.append(g)
        gene_lengths.append(gene_info[3])
        if gene_info[2] == 1:
            unspliced += 1
        else:
            spliced_gene_lengths.append(gene_info[3])

    print("Total genes found: %d" % total)
    print("Of them unspliced: %d" % unspliced)
    print(gene_types)
    print("Mean gene len %.2f" % numpy.mean(gene_lengths))
    print("Mean spliced gene len %.2f" % numpy.mean(spliced_gene_lengths))


def count_paralogs(list1, list2, gene_dict, gene_id2name, finder):
    paralog_dict = defaultdict(list)
    paralogs1 = set()
    paralogs2 = set()
    pc_paralogs1 = set()
    pc_paralogs2 = set()

    for g in list1:
        # Get paralogs
        g_id = gene_dict[g][4]
        paralogs = finder.get_paralogs(g_id, 'human')
        for p in paralogs:
            if p.target_species == 'human':
                paralog_dict[gene_id2name[p.target_id]].append(g)

    for g in list2:
        if g in paralog_dict:
            paralogs2.add(g)
            if g in gene_dict and gene_dict[g][0] == 'protein_coding':
                pc_paralogs2.add(g)
        for g1 in paralog_dict[g]:
            paralogs1.add(g1)
            if g1 in gene_dict and gene_dict[g1][0] == 'protein_coding':
                pc_paralogs1.add(g1)

    print("Paralogs in 1: %d" % len(paralogs1))
    print("Paralogs in 2: %d" % len(paralogs2))
    print("PC Paralogs in 1: %d" % len(pc_paralogs1))
    print("PC Paralogs in 2: %d" % len(pc_paralogs1))

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--de_genes", "-d", nargs='+', type=str, help="2 list of genes (up/down regulated)", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils genedb", required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    gene_dict, gene_id2name = load_genedb(args.genedb)
    gene_lists = []
    for f in args.de_genes:
        gene_list = load_genes(f)
        gene_lists.append(gene_list)
        count_stats(gene_list, gene_dict)

    if len(args.de_genes) == 2:
        finder = EnsemblParalogFinder()
        count_paralogs(gene_lists[0], gene_lists[1], gene_dict, gene_id2name, finder)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
