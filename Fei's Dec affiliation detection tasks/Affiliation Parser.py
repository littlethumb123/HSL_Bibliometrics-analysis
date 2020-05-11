import pandas as pd
import re

path = "/your/file/path.xlsx"
addre_df = pd.read_excel(path)
addparser = AddressParser()
ORG_COL = ['unrecognizable', 'univ', 'asct', 'co.', 'hosp', 'sch', 'ctr', 'clg', 'inst', 'dpt_pgm', 'div', 'region', 'country']
addre_df = pd.concat([addre_df, pd.DataFrame('', columns = ORG_COL, index = addre_df.index)], axis=1)
SPECIAL_GROUP = ['mayo clinic', 'inserm', 'cepheid', 'trento']
AVOID_Prep = ['on', 'at', 'in', '', ' ']   # Avoid preprosition in the address tokens
ENTITY_Typo = {'brasil': 'brazil',
               'pr china': 'china',
               'prchina': 'china', 
               'peoples republic of china': 'china', 
               'p r china': 'china', 
               ', uk': ', united kingdom',
               ', usa': ', united states',
               ', us': ', united states',
               'belgique': 'belgium'}
US_states_abbrev = {
    'AK': 'Alaska',
    'AL': 'Alabama',
    'AR': 'Arkansas',
    'AS': 'American Samoa',
    'AZ': 'Arizona',
    'CA': 'California',
    'CO': 'Colorado',
    'CT': 'Connecticut',
    'DC': 'District of Columbia',
    'DE': 'Delaware',
    'FL': 'Florida',
    'GA': 'Georgia',
    'GU': 'Guam',
    'HI': 'Hawaii',
    'IA': 'Iowa',
    'ID': 'Idaho',
    'IL': 'Illinois',
    'IN': 'Indiana',
    'KS': 'Kansas',
    'KY': 'Kentucky',
    'LA': 'Louisiana',
    'MA': 'Massachusetts',
    'MD': 'Maryland',
    'ME': 'Maine',
    'MI': 'Michigan',
    'MN': 'Minnesota',
    'MO': 'Missouri',
    'MP': 'Northern Mariana Islands',
    'MS': 'Mississippi',
    'MT': 'Montana',
    'NA': 'National',
    'NC': 'North Carolina',
    'ND': 'North Dakota',
    'NE': 'Nebraska',
    'NH': 'New Hampshire',
    'NJ': 'New Jersey',
    'NM': 'New Mexico',
    'NV': 'Nevada',
    'NY': 'New York',
    'OH': 'Ohio',
    'OK': 'Oklahoma',
    'OR': 'Oregon',
    'PA': 'Pennsylvania',
    'PR': 'Puerto Rico',
    'RI': 'Rhode Island',
    'SC': 'South Carolina',
    'SD': 'South Dakota',
    'TN': 'Tennessee',
    'TX': 'Texas',
    'UT': 'Utah',
    'VA': 'Virginia',
    'VI': 'Virgin Islands',
    'VT': 'Vermont',
    'WA': 'Washington',
    'WI': 'Wisconsin',
    'WV': 'West Virginia',
    'WY': 'Wyoming'
}
CA_province_abbrev = {
    'AB': 'Alberta',
    'BC': 'British Columbia',
    'MB': 'Manitoba',
    'NB': 'New Brunswick',
    'NL': 'Newfoundland and Labrador',
    'NT': 'Northwest Territories',
    'NS': 'Nova Scotia',
    'NU': 'Nunavut',
    'ON': 'Ontario',
    'PE': 'Prince Edward Island',
    'QC': 'Quebec',
    'SK': 'Saskatchewan',
    'YT': 'Yukon'
}

# regex: remove unnecessary marks
noise_regex = r"province|(?<!co)\.|\ball in\b|\bboth in\b|\ball in the\b|\bboth in the\b|ltd|&quot|[\'|\"|\+|\-]|\w*[0-9]\w*|[0-9-]+|\([^\(\)]+\)|\[[^\[\]]+\]|grid.+$"
def clean_address(address_str):
    '''
    Arg: Address string
    return clean address string
    '''
    basic_clean = re.sub(noise_regex, '', address_str)
    pattern = re.compile(r'\b(' + '|'.join(ENTITY_Typo.keys()) + r')\b')
    result = pattern.sub(lambda x: ENTITY_Typo[x.group()], basic_clean)
    
    return result

delimiter = r';|,|:| {2,}'
def token_address(address_string):
    '''
    Arg: Address string
    return token list
    '''
    temp = re.split(delimiter, clean_address(address_string))
    return [i.strip() for i in temp if i not in AVOID_Prep]

	
def check_state_abbrev(string, ret = 0):
    '''
    check 2 digits and return states or country
    Just for US and Canada
    '''
    state = ''
    country = ''
    string = string.upper().strip()
    if string.upper() in US_states_abbrev.keys():
        state = US_states_abbrev[string.upper()]
        country = 'united states'
    elif string.upper() in CA_province_abbrev.keys():
        state = CA_province_abbrev[string.upper()]
        country = 'canada'
    else:
        state = False
        country = False
    
    if ret == 0:
        return country
    else:
        return state
    

def match_org(address_token):
    '''
    Will not use regex here because it is more complex
    use predefined label to identify organization
    
    Arg: address token
    return:
    1. addre_dict: a dict{'unrecognizable
                            univ: abc, 
                            sch: abc, 
                            clg: abc, 
                            dpt: abc,
                            ctr: abc, 
                            inst: abc,
                            hosp: abc,
                            div: abc
                            co.: abc
                            dpt_pgm: abc
                            region: abc,
                            counry: abc
                            }
    2. find_org: bool, True is any field has value, otherwise false (all None)
    3. unrecognized list: 
    
    '''
   
    unrecognize_list = []   # Mark the ones that may not exist
    keyList = ORG_COL
    addre_dict = {key: '' for key in keyList}
    
    find_org = False
    try:
        # iterate the address component in the address token list
        for r, i in enumerate(address_token):
            i = i.strip()
            if any(x in i for x in ['univ', 'university','universita', 'universite', 'universidade', 'academy', 'universidad']):
                addre_dict['univ'] = i
            elif any(x in i for x in ['society', 'association', 'agency', 'organization', 'bureau', 'minister', 'ministry', 'government', 'ministerio']):
                addre_dict['asct'] = i
            elif any(x in i for x in ['hospital', 'system', 'hospices', 'hopital']+SPECIAL_GROUP):
                addre_dict['hosp'] = i
            elif any(x in i for x in ['co.','company', 'cooperation', 'foundation']):
                addre_dict['co.'] = i
            elif any(x in i for x in ['school', 'professor', 'faculty', 'faculdade', 'faculte']):
                addre_dict['sch'] = i
            elif any(x in i for x in ['center','centre', 'centro']):
                addre_dict['ctr'] = i
            elif 'college' in i:
                addre_dict['clg'] = i
            elif any(x in i for x in ['institute', 'institue', 'instituto', 'institut']):
                addre_dict['inst'] = i
            elif any(x in i for x in ['department', 'dipartimento', 'program', 'programme', 'dept.', 'dept']):
                addre_dict['dpt_pgm'] = i
            elif any(x in i for x in ['group','division', 'divison', 'laboratory']):
                addre_dict['div'] = i
            else:
                if i != '':
                    unrecognize_list.append(i)    # (location, content)
        
        # Check if find any org in at least one field
        if len(set(addre_dict.values())) > 1:
            # found at least one organization, in addition to ''
            find_org = True
        
        # Check if any city and state information here. Better not check city with api
        
        # if more than one items not recognize
        if len(unrecognize_list)>0:
            # pop out the last value of this list, which is the last position in such address tokens, as country
            addre_dict['country'] = unrecognize_list.pop(-1)
            
            if len(unrecognize_list)>0:
                addre_dict['region'] = unrecognize_list.pop(-1)
            
        # if still have item, then unrecognized
        if len(unrecognize_list)>0:
            # if still have values in unrecognized list, put it in unrecognized, else this field is equal to ''
            addre_dict['unrecognizable'] = ', '.join(unrecognize_list)            
            
    except Exception as e:
        print("cannot address: ", address_token)
        traceback.print_exc()
        pass
    
    return addre_dict, unrecognize_list, find_org


# Example
def main():

	address = "birmingham acute care research group, institute of inflammation and ageing, college of medical and dental sciences, university of birmingham, birmingham, uk"
	addr_token = token_address(address)
	address_dict, _, _ = match_org(addr_token)
	print(address_dict)

if __name__=="__main__": 
    main() 

# Output: 
# {'unrecognizable': '', 
	# 'univ': 'university of birmingham', 
	# 'asct': '', 
	# 'co.': '', 
	# 'hosp': '', 
	# 'sch': '', 
	# 'ctr': '', 
	# 'clg': 'college of medical and dental sciences', 
	# 'inst': 'institute of inflammation and ageing', 
	# 'dpt_pgm': '', 
	# 'div': 'birmingham acute care research group', 
	# 'region': 'birmingham', 
	# 'country': 'united kingdom'}








