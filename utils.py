from Bio import Entrez
import xml.etree.ElementTree as ET
# reference: https://github.com/titipata/pubmed_parser/blob/master/pubmed_parser/pubmed_oa_parser.py

Entrez.email = "poppincorngofurther@gmail.com"

def retrieve_affil(pmcid):
    '''
    This is just processing single article. 
    
    para: input the pubmed central id of the article which miss affiliation information
    return: return a dictionary, e.g., 
        {"forename": "ret_from_pmc",
         "firstname": firstname,
         "lastname": lastname, 
         "affiliation": a list of affiliation of the current author}
    
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
