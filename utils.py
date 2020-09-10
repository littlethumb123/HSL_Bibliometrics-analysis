from Bio import Entrez
import xml.etree.ElementTree as ET
from lxml import html
import requests
import time
import random
from collections import deque
import io
import pandas as pd
# reference: https://github.com/titipata/pubmed_parser/blob/master/pubmed_parser/pubmed_oa_parser.py

Entrez.email = "hslexample@gmail.com"
user_agent_list = [
   #Chrome
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.113 Safari/537.36',
    'Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.90 Safari/537.36',
    'Mozilla/5.0 (Windows NT 5.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.90 Safari/537.36',
    'Mozilla/5.0 (Windows NT 6.2; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.90 Safari/537.36',
    'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/44.0.2403.157 Safari/537.36',
    'Mozilla/5.0 (Windows NT 6.3; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.113 Safari/537.36',
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/57.0.2987.133 Safari/537.36',
    'Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/57.0.2987.133 Safari/537.36',
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/55.0.2883.87 Safari/537.36',
    'Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/55.0.2883.87 Safari/537.36',
    #Firefox
    'Mozilla/4.0 (compatible; MSIE 9.0; Windows NT 6.1)',
    'Mozilla/5.0 (Windows NT 6.1; WOW64; Trident/7.0; rv:11.0) like Gecko',
    'Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; WOW64; Trident/5.0)',
    'Mozilla/5.0 (Windows NT 6.1; Trident/7.0; rv:11.0) like Gecko',
    'Mozilla/5.0 (Windows NT 6.2; WOW64; Trident/7.0; rv:11.0) like Gecko',
    'Mozilla/5.0 (Windows NT 10.0; WOW64; Trident/7.0; rv:11.0) like Gecko',
    'Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.0; Trident/5.0)',
    'Mozilla/5.0 (Windows NT 6.3; WOW64; Trident/7.0; rv:11.0) like Gecko',
    'Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; Trident/5.0)',
    'Mozilla/5.0 (Windows NT 6.1; Win64; x64; Trident/7.0; rv:11.0) like Gecko',
    'Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.1; WOW64; Trident/6.0)',
    'Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.1; Trident/6.0)',
    'Mozilla/4.0 (compatible; MSIE 8.0; Windows NT 5.1; Trident/4.0; .NET CLR 2.0.50727; .NET CLR 3.0.4506.2152; .NET CLR 3.5.30729)'
]

def XML_request(pmids):
    """
    Request data from entrez using request method
    param: pmid or a list of pmid
    :return: a dictionary of pmid: html_tree
    """
    PMD_LINK = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id={}"
    full_xml = {}
    if isinstance(pmids, list) and pmids:
        dq = deque(pmids)
        while dq:
            pmid = dq.popleft()
            try:
                resp = requests.get(PMD_LINK.format(pmid), headers={'user-agent': random.choice(user_agent_list)})
                full_xml[pmid] = html.fromstring(resp.content)
            except Exception as e:
                if not isinstance(e, requests.exceptions.RequestException):
                    print("ERROR {}: {}".format(pmid, e))
                dq.append(pmid)
        print('%d records left\r'%len(dq), end = "")
        time.sleep(1)
    elif isinstance(pmids, str):
        resp = requests.get(PMD_LINK.format(pmids), headers={'user-agent': random.choice(user_agent_list)})
        full_xml[pmids] = html.fromstring(resp.content)
    else:
        raise TypeError("Please import single pmid or pmid list")

    return full_xml


def parse_handle(handle):
    """
    parse handle returned by efetch with xml format
    :param: handle, entrez handle object
    :return: df, dataframe having publication data with a column ["PMID", "TI", "AB", "AU", "AF", "PD", "KW"]

    """
    column = ["PMID", "TI", "AB", "AU", "AF", "PD", "KW"]

    if not isinstance(handle, io.TextIOWrapper):
        raise TypeError("Input should be handle object")

    fetchrecords = Entrez.read(handle)
    if not fetchrecords['PubmedArticle']:
        raise ValueError("File is empty")

    df = pd.DataFrame(columns=column)
    for i in fetchrecords['PubmedArticle']:

        rec = i['MedlineCitation']
        rec_art = rec['Article']

        au_list = [au_af['ForeName'] + " " + au_af['LastName']
                   if au_af['ForeName'] and au_af['ForeName'] else "name NA"
                   for au_af in rec_art['AuthorList']]
        aff_list = [au_af['AffiliationInfo'][0]['Affiliation']
                    if au_af['AffiliationInfo'][0] else "affiliation NA"
                    for au_af in rec_art['AuthorList']]

        pmid = "PMID NA"
        if rec['PMID']:
            pmid = rec['PMID']

        title = "title NA"
        if rec_art['ArticleTitle']:
            title = rec_art['ArticleTitle']

        date = "Date NA"
        if rec_art['ArticleDate']:
            date = rec_art['ArticleDate'][0]['Month'] + \
                   "/" + rec_art['ArticleDate'][0]['Day'] + \
                   "/" + rec_art['ArticleDate'][0]['Year']

        abstract = "Abstract NA"
        if 'Abstract' in rec_art.keys():
            abstract = rec_art['Abstract']['AbstractText']

        keywords = "Keywords NA"
        if rec['KeywordList']:
            keywords = "; ".join(j for j in rec['KeywordList'][0])
        rec_dict = {
            "PMID": pmid,
            "TI": title,
            "AB": abstract,
            "PD": date,
            "AU": "; ".join(j for j in au_list),
            "AF": "; ".join(j for j in aff_list),
            "KW": keywords
        }

        df = df.append(rec_dict, ignore_index=True)
        print("Data are ready to export")
    return df



def retrieve_affil(pmcid):
    '''
    This is for retrieve lacking author information in pubmed xml, from pubmed central databases. 
    Noted that it is just for processing single article. 
    
    para: input the pubmed central id of the article which miss affiliation information
    return: return a dictionary, e.g., 
            {"forename": "ret_from_pmc",
             "firstname": firstname,
             "lastname": lastname,
             "affiliation": a list of affiliation of the current author}``
    
    Note: limits to access api
    Any site (IP address) posting more than 3 requests per second to the E-utilities 
    without an API key will receive an error message. 
    By including an API key, a site can post up to 10 requests per second by default.
    
    '''
    # process single xml
    print("In side ret_affili and pmcid is: ", pmcid)
    try:
        fetch_handle = Entrez.efetch(db="pmc", id = pmcid, retmode="xml") 
        record = fetch_handle.read()  # a dictionary with xml structure
    except ValueError:
        raise ValueError("The pmc_id is empty")
    
    tree = ET.fromstring(record)
    tree.findall(".")
    
    #obtain affiliation information
    aff_id_address_dic = {"no_id_available":[]} # map aff_id with aff_name
    for aff in tree.iter('aff'):
        if aff.attrib.get('id') is not None: # cound be empty, in pubmed xml example
            aff_id = aff.attrib.get('id') # here is aff_id, originally is {id:aff1}
        else:
            aff_id = "no_id_available"
        
        if aff.find('addr-line') is not None:
            aff_address = aff.find('addr-line').text
        else:
            aff_address = "".join(aff.itertext())
            #print("aff.text", aff_address)
        if aff_id == "no_id_available":
            # no reference id is available, then assign all affiliation to every single author
            aff_id_address_dic[aff_id].append(aff_address)
        else:
            # any reference id then map id to the affiliation
            aff_id_address_dic.update({aff_id:aff_address})

    #print("here is aff and id mapping:", aff_id_address_dic)

    
    # Obtain author information
    tree.findall(".")
    tree_author = tree.findall('.//contrib-group//contrib[@contrib-type="author"]')
    author_aff_mapping_list = list()
    for author in tree_author:
        aff_of_this_auth = []
        if len(author.findall('xref[@ref-type="aff"]'))!=0:
            # in this case where affiliation and author are not matched
            author_aff_xref = author.findall('xref[@ref-type="aff"]')
            author_affi_refid = [a.attrib['rid'] for a in author_aff_xref]  # all affiliations of an author
            
            # map author reference to corresponding
            for i in author_affi_refid:
                t = aff_id_address_dic.get(i)
                aff_text = t.replace('\n', ' ').replace('\t', ' ').strip()
                aff_of_this_auth.append(aff_text)
            #print("here are ids reference and: ", aff_of_this_auth)
        elif author.find('aff') is not None:
            # in this case where affiliation information is attached with author's information
            aff_text = author.find('aff').text.replace('\n', ' ').replace('\t', ' ').strip()
            aff_of_this_auth = [aff_text]
            #print("here no id reference mapping but this author affi is: ", aff_of_this_auth)
        else: 
            # in this case affiliation is neither mapped or attached to author, then use the previously extracted affiliation information:
            try:
                aff_of_this_auth.append(aff_id_address_dic["no_id_available"])
            except KeyError as error:
                print("  Here is the key error raised in no_id ", error)
#                 print("  Affiliation information is not available in pmc", pmcid)
        try:
            author_aff_mapping_list.append(
                                        {"firstname": "ret_from_pmc",
                                        "forename": author.find('name/given-names').text,
                                        "lastname": author.find('name/surname').text, 
                                        "affiliation": aff_of_this_auth})

        except AttributeError as error:
             print("     There might be something wrong with Nonetype find", error)
#             print("     author.find('name/given-names').text: ", author.find('name/given-names'))
#             print("     author.find('name/surname').text: ", author.find('name/surname'),)
            author_aff_mapping_list.append({"firstname": "ret_from_pmc",
                                    "forename": "ret_from_pmc",
                                    "lastname": "ret_from_pmc", 
                                    "affiliation": "ret_from_pmc"})
            
    #print("{} Add the following list {}". format(pmcid, author_aff_mapping_list))
    return author_aff_mapping_list
