import sys
from traceback import print_exc
import time
import socket
try:
    from urllib2 import urlopen
    import urllib
except:
    from urllib.request import urlopen
    import urllib.request as urllib

import xml.etree.ElementTree as ET


ORGANISM = "Mus%20musculus"
DESCRIPTION = "splicing factor"


def try_send_request(url):
    attempts = 0
    response = None
    connection_errors = 0
    while attempts < 3:
        try:
            request = urlopen(url)
            connection_errors = 0
            response = request.read()
            if not isinstance(response, str):
                response = response.decode('utf-8')
            if response is None or 'ERROR' in response:
                request.close()
                raise Exception
            break
        except Exception:
            attempts += 1
            if attempts >= 3:
                connection_errors += 1
                if connection_errors >= 3:
                    print('Cannot establish connection to NCBI')
                return None
            # NCBI recommends users post no more than three URL requests per second, so adding artificial 1-sec delay
            time.sleep(1)
    return response


def get_gene_ids(gene_name):
    query = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=({gene_name}[gene])%20AND%20({organism}[orgn])%20AND%20alive[prop]%20NOT%20newentry[gene]".\
            format(gene_name=gene_name, organism=ORGANISM)
    response = try_send_request(query)

    if response is None:
        print('Failed to retrieve gene information')
        return None

    xml_tree = ET.fromstring(response)

    if xml_tree.find('Count').text == '0':  # Organism is not found
        print("Gene " + gene_name + " was not found")
        return None

    gene_id_list = xml_tree.find('IdList').findall('Id')
    gene_ids = []
    for gid in gene_id_list:
        gene_ids.append(gid.text)

    return gene_ids


def check_gene_description(ncbi_gene_id, gene_name):
    query = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={gene_id}".format(gene_id=ncbi_gene_id)
    response = try_send_request(query)
    
    gene_block_found = False
    description = ""
    locus = ""
    in_syn_block = False
    synonyms = ""
    for l in response.split('\n'):
        if not gene_block_found and l.find("gene {") != -1:
            gene_block_found = True
            continue
        if not gene_block_found:
            continue
        #print(l)
            
        if l.find("desc") != -1:
            description = l

        if l.find("locus") != -1:
            locus = l

        if l.find("syn") != -1:
            in_syn_block = True
            continue

        if in_syn_block:
            if l.find('}') != -1:
                in_syn_block = False
                continue
            synonyms += l.strip()

        if len(description) > 0 and len(locus) > 0  and len(synonyms) > 0 and not in_syn_block:
            break

    return description.lower().find(DESCRIPTION.lower()) != -1 and (locus.lower().find(gene_name.lower()) != -1 or synonyms.lower().find(gene_name.lower()) != -1)
    

def main():
    if len(sys.argv) != 2:
        print("Usage: " + sys.argv[0] + " <file with gene list one per line>  > <genes that have " + DESCRIPTION + " in their description>")
        exit(0)

    for l in open(sys.argv[1]):
        gene_name = l.strip()
        if len(gene_name) == 0:
            continue
        ids = get_gene_ids(gene_name)
        for gid in ids:
            if check_gene_description(gid, gene_name):
                print(gene_name)
    

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
