from tqdm import tqdm
from pandas import pd
from Bio import Medline
import re
import requests
import time
import random
from collections import deque

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
# Regex expression to find R; or add other targets regular expression
RE_MATCH = r"([A-Za-z ]*(<(?:kbd|em)>R<(?:\/em|\/kbd)>|R: A language and environment for statistical computing" \
           r"| R[,]?[ \n](?:[(]?[Vv](?:ersion)?)?[ (]?[0-9.]{1,5}|R[\n ][Ss]tatistical software|R[\n ]Foundation" \
           r"|R[- \n][Ss]tudio|R[- \n][Ll]ibrary|R[- \n][Ll]ibraries|R[- \n][Pp]ackage|R[- \n][Pp]rogram|R[- \n][Ss]cript|" \
           r"R[- \n][Ss]oftware|R[- \n][Cc]ore [Tt]eam|R[- \n][Cc]ode)[A-Za-z ,]*)"

PMC_PATH = 'https://www.ncbi.nlm.nih.gov/pmc/articles/'
AFFILIATION = ["Chapel Hill"]   # The target affiliations

def main():
    # Here is the targeted affiliation
    '''
    The function will detect the use of R in a publications from particular affiliations
    :return: None; a dataframe will be exported as R_check_output.xlsx
    '''

    # Create record xml list
    rec_list = []

    # Create pmc id list
    pmc_list = deque()

    # read medline file exported from pubmed
    with open('data_sample/pmc_result_medline.txt') as handle:
        for record in Medline.parse(handle):
            pmc_list.append(record['PMC'])
            rec_list.append(record)

    # create pandas for R and SAS detect
    df_R = pd.DataFrame(index=pmc_list, columns=["Target_Aff", "use_R", "R"])

    # store the scraped text
    full_xml = {}

    # scrape full text from pubmed using requests. This will take a long time to request content from PubMed.
    while pmc_list:
        pmc = pmc_list.popleft()
        try:
            resp = requests.get(PMC_PATH + pmc, headers={'user-agent': random.choice(user_agent_list)})
            xml_text = str(resp.text)
            full_xml[pmc] = xml_text
        except Exception as e:
            if not isinstance(e, requests.exceptions.RequestException):
                print("ERROR {}: {}".format(pmc, e))
            pmc_list.append(pmc)

        print('%d records left\r' % len(pmc_list), end="")
        time.sleep(1)


    # find R with regex

    # get pmcid that has no author and did not use R
    no_info = []
    # iterate the medline records and find usage of R and corresponding affiliation information
    for rec in tqdm(rec_list):
        pmc = rec.get('PMC')
        try:
            # get the matched text with R showing up
            xml_text = full_xml[pmc]
            match = re.findall(RE_MATCH, xml_text)
            if len(match)>0:
                # get the markers of R
                df_R.at[pmc, "R"] = [i[0] for i in match]
                df_R.at[pmc, "use_R"] = 1
            else:
                df_R.at[pmc, "use_R"] = 0
                no_info.append(pmc)
            if 'AD' in rec.keys():
                aff_list = rec.get('AD')
                for index, aff in enumerate(aff_list):
                    if any([i in aff for i in AFFILIATION]):
                        df_R.at[pmc, "Target"] = 1
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
    # df_R.to_excel("R_check_output.xlsx")

if __name__ == "__main__":
    main()