from tqdm import tqdm
import pandas as pd
from Bio import Medline
import re
import requests
import time
import random
from collections import deque
from utils import user_agent_list  # import browser to iteratively access pubmed database server
from bs4 import BeautifulSoup

# Regex expression to find R; or add other targets regular expression
RE_MATCH = r"([^\s]{0,40}[ \(\{]?)(R: A[ \n][Ll]anguage[ \n]and[ \n][Ee]nvironment| R[,]?[ \n](?:[(]?[Vv](?:ersion)?)?[ (]?\d\.\d{1,2}\.\d|R[\n ][Ss]tatistical software|R[\n ]Foundation|R[- \n][Ss]tudio|R[- \n][Ll]ibrary|R[- \n][Ll]ibraries|R[- \n][Pp]ackage| R[- \n][Pp]rogram| R[- \n][Ss]cript| R[- \n][Ss]oftware|R[- \n][Cc]ore [Tt]eam|R[- \n][Cc]ode\)[A-Za-z ,]*|www\.R-project\.org|RStudio)([ \(\{]?[^\s]{0,40})"
PMC_PATH = 'https://www.ncbi.nlm.nih.gov/pmc/articles/'
AFFILIATION = ["Chapel Hill"]   # The target affiliations

# the path to medline.text files
FILE_PATH = 'data_sample/pubmed-Rprogram-set.txt'



def main():
    # Here is the targeted affiliation
    '''
    The function will detect the use of R in a publications from particular affiliations
    :return: None; a dataframe will be exported as R_check_output.xlsx
    '''

    # create record xml list
    rec_list_ = []

    # create pmc id list
    pmc_list_ = []


    # read medline file exported from pubmed
    with open(FILE_PATH) as handle:
        for record in Medline.parse(handle):
            pmc_list_.append(record['PMC'])
            rec_list_.append(record)

    # for demonstration, we get first 10 elements
    pmc_list = deque(pmc_list_[:10])
    rec_list = rec_list_[:10]

    # create pandas for R and SAS detect
    df_R = pd.DataFrame(index=pmc_list, columns=["Target_Aff", "use_R", "R"])

    # store the scraped text
    full_xml = {}
    p_bar = tqdm(total=len(pmc_list), position=0, leave=True)
    # scrape full text from pubmed using requests. This will take a long time to request content from PubMed.
    while pmc_list:
        pmc = pmc_list.popleft()
        print('%d records left\r' % len(pmc_list), end="")
        try:
            resp = requests.get(PMC_PATH + pmc, headers={'user-agent': random.choice(user_agent_list)})
            xml_text = BeautifulSoup(resp.text, features="html.parser").get_text(" ", strip=True)
            full_xml[pmc] = xml_text
            p_bar.update(1)
        except Exception as e:
            if not isinstance(e, requests.exceptions.RequestException):
                print("ERROR {}: {}".format(pmc, e))
            pmc_list.append(pmc)
        time.sleep(1)


    # find R with regex and, if any, corresponding affiliation information

    # get pmcid that has no author and did not use R
    no_info = []
    # iterate the medline records and find usage of R and corresponding affiliation information
    print("Checking R use...")
    for rec in tqdm(rec_list, position=0, leave=True):
        pmc = rec.get('PMC')
        try:
            # get the matched text with R showing up
            xml_text = full_xml[pmc]
            match = re.findall(RE_MATCH, xml_text)
            if len(match)>0:
                # get the markers of R
                df_R.at[pmc, "R"] = [i[1] for i in match]
                df_R.at[pmc, "use_R"] = 1
            else:
                df_R.at[pmc, "use_R"] = 0
                no_info.append(pmc)

            # get affiliation information
            if 'AD' in rec.keys():
                aff_list = rec.get('AD')
                for index, aff in enumerate(aff_list):
                    if any([i in aff for i in AFFILIATION]):
                        df_R.at[pmc, "Target_Aff"] = 1
                    else:
                        no_info.append(pmc)
            else:
                # no author information
                no_info.append(pmc)

            # alternatively, find the information of missing pmcid in scraped text
            for pmcid in set(no_info):
                if any([i in xml_text for i in AFFILIATION]):
                    df_R.at[pmcid, "Target_Aff"] = 1
                else:
                    df_R.at[pmcid, "Target_Aff"] = 0

        except Exception as e:
            print ("ERROR {}: {}".format(pmc, e))

    # print(df_R)
    df_R.to_excel("data_sample/R_check_output.xlsx")

if __name__ == "__main__":
    main()