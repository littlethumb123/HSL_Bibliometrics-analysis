import pandas as pd
import pubmed_parser as pp
import json
import re
import csv

'''
Extract the pubmed id and doi of the articles with missing affiliation and by targeted affiliations 
12/4/2019
'''

target_uni = ["Johns Hopkins University", 
              "Duke University", 
              "Harvard University", 
              "National Institutes of Health", 
              "North Carolina State University"]

# import files
pubmed_dic  = pp.parse_medline_xml(path, author_list=True)  # here path is the xml file path, a placeholder here
#print(json.dumps(pubmed_dic, indent = 4))  #explore the parsed data if needed

miss_affi_article_info = []
target_affi_article_info = []
for article in pubmed_dic: # a list of article info
    pmid = article['pmid']
    doi = article['doi']
    first_author_affiliation = article['authors'][0]['affiliation']
    i_miss_affi = 0
    target_affi = []
    for author in article['authors']: # a list of author info in a particular article
        # Task_1
        if '' in author['affiliation']:
            i_miss_affi = i_miss_affi+1
        # Task_2
        if any(j in author['affiliation'] for j in target_uni):
            target_affi.append(author['affiliation'])
    if i_miss_affi>0: 
        miss_affi_article_info.append((pmid, doi, first_author_affiliation))   
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