from tqdm import tqdm
from pandas import pd
from Bio import Medline
import re

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
# Regex expression to find R 
re_match = r"([A-Za-z ]*(<(?:kbd|em)>R<(?:\/em|\/kbd)>|R: A language and environment for statistical computing| R[,]?[ \n](?:[(]?[Vv](?:ersion)?)?[ (]?[0-9.]{1,5}|R[\n ][Ss]tatistical software|R[\n ]Foundation|R[- \n][Ss]tudio|R[- \n][Ll]ibrary|R[- \n][Ll]ibraries|R[- \n][Pp]ackage|R[- \n][Pp]rogram|R[- \n][Ss]cript|R[- \n][Ss]oftware|R[- \n][Cc]ore [Tt]eam|R[- \n][Cc]ode)[A-Za-z ,]*)"
PMC_PATH = 'https://www.ncbi.nlm.nih.gov/pmc/articles/'


rec_list = []
pmc_list = []
# read medline file for parsing
with open('file/path') as handle:
    for record in Medline.parse(handle):
        pmc_list.append(record['PMC'])
        rec_list.append(record)


# create pandas for R and SAS detect
df_R = pd.DataFrame(index = pmc_list, columns = ["UNC", "use_R", "R"])  
#df_SAS = pd.DataFrame(index = pmc, columns = ["UNC", "use_SAS", "SAS"])
# store the scraped text
full_text = {}  

# scrape full text from pubmed
while 1:
    error_id = []
    for pmc in tqdm(pmc_list):
        #pmc = rec.get('PMC')
        try:
            resp = requests.get(PMC_PATH + pmc, headers = random_header(user_agent_list))
            xml_text = str(resp.text)
            full_text[pmc] = xml_text
        except Exception as e:
            print ("ERROR {}: {}".format(pmc, e))
            error_id.append(pmc)
        time.sleep(1)
    if len(error_id) == 0:
        break
    else:
        redo_id = error_id

# find R with regex
no_AD = []  # get pmcid that has no author info
no_R = []   # get pmcid that has no R used
no_unc = [] # get pmcid that is not UNC authored
error = []  # get pmcid with errors

for rec in tqdm(rec_list):
    pmc = rec.get('PMC')
    try:
        xml_text = full_text[pmc]
        match = re.findall(re_match, xml_text)
        if len(match)>0:
            df_R.at[pmc, "R"] = [i[0] for i in match]
            df_R.at[pmc, "use_R"] = 1
        else:
            df_R.at[pmc, "use_R"] = 0
            no_R.append(pmc)
        if 'AD' in rec.keys():
            aff_list = rec.get('AD')
            add_auth = []
            for index, aff in enumerate(aff_list):
                if any([i in aff for i in AFFILIATION]):
                    df_R.at[pmc, "UNC"] = 1
                else:
                    no_unc.append(pmc)
        else:
            no_AD.append(pmc)
        
        # find missing pmcid in scraped text
        missunc = no_AD+no_unc
        for pmcid in set(missunc):
            if any([i in xml_text for i in AFFILIATION]):
                df_R.at[pmcid, "UNC"] = 1
            else:
                df_R.at[pmcid, "UNC"] = 0
    
    except Exception as e:
        print ("ERROR {}: {}".format(pmc, e))
        error.append(pmc)
df_R.to_excel("R_UNCcheck_output_check.xlsx")
