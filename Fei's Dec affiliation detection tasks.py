import pandas as pd
import pubmed_parser as pp
import json
import re
import csv
from .utils import *
"""
Extract the pubmed id and doi of the articles with missing affiliation and by targeted affiliations 
12/4/2019
"""

target_uni = ["Johns Hopkins University", 
              "Duke University", 
              "Harvard University", 
              "National Institutes of Health", 
              "North Carolina State University"]

# import files
pubmed_dic  = pp.parse_medline_xml(path, author_list=True)  # here path is the xml file path, a placeholder here
#print(json.dumps(pubmed_dic, indent = 4))  #explore the parsed data if needed

"""
Basic idea:
Medline xml in, through pmc search and process pmc_oa_xml

12/7: 
aff reference id is not consistent and all pmcids are not captured 
Here show two types of oa xml formats https://www.ncbi.nlm.nih.gov/pmc/pmcdoc/tagging-guidelines/article/style.html
"""

miss_affi_article_info = []
target_affi_article_info = []
artcile_count = len(pubmed_dic)
for article in pubmed_dic: # a list of article info
    pmid = article['pmid']
    doi = article['doi']  # Could be empty
    pmcid = article['pmc_id']  # Could be empty
    #first_author_affiliation = article['authors'][0]['affiliation']
    i_miss_affi = 0
    target_affi = []
    for author in article['authors']:
        # Task_1
        if author['affiliation'] == '':
            i_miss_affi = i_miss_affi+1
        # Task_2  better separate from task 1 
        if any(j in author['affiliation'] for j in target_uni):
            target_affi.append(author['affiliation'])
    if i_miss_affi>0:
        print("here lacks ", i_miss_affi)
        if pmcid!='':
            article['authors'] = retrieve_affil(pmcid) # here do the search and return affiliation
            print("can add and add the following affiliation: ", article['authors'])
        else:
            print("cannot add because pmcid is empty")
            miss_affi_article_info.append((pmid, doi, 'first_author_affiliation'))   
    if len(target_affi)!= 0:
        target_affi_article_info.append((pmid, doi, target_affi))

# export as csv files
with open('319task.csv','w', encoding='utf-8') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['pubmed_id','doi', 'target_affiliation'])
    for row in target_affi_article_info:
        csv_out.writerow(row)
with open('234task.csv','w', encoding='utf-8') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['pubmed_id','doi', 'first_author affiliation'])
    for row in miss_affi_article_info:
        csv_out.writerow(row)
