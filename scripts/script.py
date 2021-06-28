import pandas as pd
from Bio import Entrez
from Bio import Medline
from Bio.Entrez import efetch, read
import os
import json

def search(query, retmax = 100):
    Entrez.email = 'ex.aa@mail.edu'
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax=retmax,
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    pmids = results['IdList']
    return pmids

def fetch_pubmed(id_list):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.efetch(db='pubmed',
                            rettype="medline",
                            id=id_list)
    results = Medline.parse(handle)
    results = list(results)
    return results

def artical_parser(document):
    
    keys_document = document.keys()
    parsed_document = dict()
    if 'PMID' in keys_document:
        parsed_document['PMID'] = document['PMID']
    if 'TI' in keys_document:
        parsed_document['TI'] = document['TI']
    if 'AB' in keys_document:
        parsed_document['AB'] = document['AB']
    if 'JT' in keys_document:
        parsed_document['JT'] = document['JT']
    if 'FAU' in keys_document:
        parsed_document['FAU'] = document['FAU']
    
    return parsed_document


def medline_parser(documents, path_for_saved_document='data/'):
    for document in documents:
        pmid = document['PMID']
        document_dictionary = artical_parser(document)
        with open(f'{path_for_saved_document}{pmid}.json', "w") as file:
            json.dump(document_dictionary,file)

def main(search_term):
    pmids = search(search_term)
    fetche_pubmeds = fetch_pubmed(pmids)
    medline_parser(fetche_pubmeds)