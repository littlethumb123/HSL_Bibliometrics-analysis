from Bio import Entrez
import xml.etree.ElementTree as ET
# reference: https://github.com/titipata/pubmed_parser/blob/master/pubmed_parser/pubmed_oa_parser.py

Entrez.email = "poppincorngofurther@gmail.com"

def retrieve_affil(pmcid):
    '''
    This is just processing single article. 
    
    para: input the pubmed central id of the article which miss affiliation information
    return: return a list, e.g., 
        ["Lastname, firstname", [aff1, aff2...]], ["Lastname, firstname", [aff1, aff2...]]
    
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
        print("The pmc_id is empty")
    
    tree = ET.fromstring(record)
    tree.findall(".")
    
    #obtain affiliation information
    aff_id_address_dic = {} # map aff_id with aff_name
    for aff in tree.iter('aff'):
        if aff.attrib.get('id') is not None: # cound be empty, in pubmed xml example
            aff_id = aff.attrib.get('id') # here is aff_id, originally is {id:aff1}
        else:
            aff_id = ""  
        if aff.find('addr-line') is not None:
            aff_address = aff.find('addr-line').text
        else:
            aff_address = "".join(aff.itertext())
            #print("aff.text", aff_address)
        aff_id_address_dic.update({aff_id: aff_address})

    #print("here is aff and id mapping:", aff_id_address_dic)


    # Obtain author information
    tree.findall(".")
    tree_author = tree.findall('.//contrib-group//contrib[@contrib-type="author"]')
    author_aff_mapping_list = list()
    for author in tree_author:

        if len(author.findall('xref[@ref-type="aff"]'))!=0:
            # in this case where affiliation and author information is separated
            author_aff_xref = author.findall('xref[@ref-type="aff"]')
            author_affi_refid = [a.attrib['rid'] for a in author_aff_xref]  # all affiliations of an author
            aff_of_this_auth = []
            for i in author_affi_refid:
                t = aff_id_address_dic.get(i)
                aff_text = t.replace('\n', ' ').replace('\t', ' ').strip()
                aff_of_this_auth.append(aff_text)
            #print("here are ids reference and: ", aff_of_this_auth)
        elif author.find('aff') is not None:
            aff_text = author.find('aff').text.replace('\n', ' ').replace('\t', ' ').strip()
            aff_of_this_auth = [aff_text]
            #print("here no id reference mapping but this author affi is: ", aff_of_this_auth)
        else: 
            aff_of_this_auth = []
            print("Affiliation information may not be available in pmc")
        try:
            author_aff_mapping_list.append(
                                    {"forename": "ret_from_pmc",
                                    "firstname": author.find('name/given-names').text,
                                    "lastname": author.find('name/surname').text, 
                                    "affiliation": aff_of_this_auth})
        except AttributeError:
            print("Name information may not be available in pmc")
            author_aff_mapping_list.append({"forename": "ret_from_pmc",
                                    "firstname": "ret_from_pmc",
                                    "lastname": "ret_from_pmc", 
                                    "affiliation": "ret_from_pmc"})
    return author_aff_mapping_list
