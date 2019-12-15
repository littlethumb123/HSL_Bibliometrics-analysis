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

'''
Basic idea:
Medline xml in, through pmc search and process pmc_oa_xml

12/7: 
aff reference id is not consistent and all pmcids are not captured 
Here show two types of oa xml formats https://www.ncbi.nlm.nih.gov/pmc/pmcdoc/tagging-guidelines/article/style.html

12/14 updates: just add the missing authors
remove ('forename': '', 'firstname': '', 'lastname': '', 'affiliation': '')

'''
miss_affi_article_info = []  
target_affi_article_info = {key:[] for key in target_uni} # add target affiliation as key and 
artcile_count = len(pubmed_dic)

# code testing data
complete_add_article = []
fail_add_article = []
miss_article_before_add = set()
for article in pubmed_dic: # a list of article info
    pmid = article['pmid']
    doi = article['doi']  # Could be empty
    pmcid = article['pmc_id']  # Could be empty
    #first_author_affiliation = article['authors'][0]['affiliation']
    target_affi = []
    article_affi = []  # particularly for task2 to collect affiliation information list
    for author in article['authors']:
        if author['lastname']!="" and author['forename']!="":
            article_affi.append(author['affiliation'])
        else:
            article['authors'].remove(author)
    if "" in article_affi:
        miss_article_before_add.add(pmid)
        if pmcid!='':
            article['authors'] = retrieve_affil(pmcid) # here do the search and retreive affiliation information
            article_affi.clear()
            article_affi = [author['affiliation'] for author in article['authors']]
            complete_add_article.append(pmid)
        else:
            print("{} cannot perform pmc search because pmcid isnot available".format(pmid))
            # Task_1
            fail_add_article.append(pmid)
            miss_affi_article_info.append((pmid, doi))
    
    # Task_2  better separate from task 1 can be performed after affiliation completion
    for j in target_uni:
        for i in article_affi:
            if type(i) == str:  # deal with string type affiliation
                if j in i and pmid not in target_affi_article_info[j]:
                    target_affi_article_info[j].append(pmid)
            else: # deal with list-like affiliation
                for ii in i:
                    if j in ii and pmid not in target_affi_article_info[j]:
                        target_affi_article_info[j].append(pmid)
                

# # export as csv files
# with open('319task.csv','w', encoding='utf-8') as out:
#     csv_out=csv.writer(out)
#     csv_out.writerow(['pubmed_id','doi', 'target_affiliation'])
#     for row in target_affi_article_info:
#         csv_out.writerow(row)

