## Retrieving Records from PubMed with Python

Two ways to access to Entrez database [1](#_ftn1) with python. E-utilities and entrezpy. 

### E-utilities 

E-utility is a set of programs that interface Entrez database for users in National Center 
for Biotechnology Information (NCBI), including PubMed. The programs lies in Biopython library 
and execute functions like EInfo (get database statistics),  ESearch (perform text searches), 
EPost (uploads UID), ESummary (downloads document summary), EFetch (data record downloads), 
ELink (Links related records in Entrez) (see [official video](https://www.youtube.com/watch?v=BCG-M5k-gvE), 
[official tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf]), [API document](https://biopython.org/DIST/docs/api/)). 
Under the hood, how these programs works? When you use the PubMed web interface to search "prostate cancer", 
you may notice "https://pubmed.ncbi.nlm.nih.gov/?term=prostate+cancer" showing up in URL bar. You can simply 
interpret it as the website creating a request for you to servers to retrieve data from PubMed database. 
Similarly, the programs post and retrieve data by sending a **URL** to the server, but with the programs you use.

The E-utitlies return records information in a particular format, xml. see example. Tools are needed to parse these files 
to obtain a readable format (see [XML element description and attributes](https://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html#authorlist)) 
and [xml example](https://www.ncbi.nlm.nih.gov/pmc/pmcdoc/tagging-guidelines/article/style.html).  




### Entrezpy

Entrezpy is a recently released python library to dynamically interact with the NCBI Entrez databases. 
This library, similar to Biopython, can interact with Entrez database via E-utility. 
Specifically, it supports the above-mentioned e-utility programs and, beyond that, enables an all-in-one "combo" 
for users to address xml retrieval and response from the server (see [Publication](https://academic.oup.com/bioinformatics/article/35/21/4511/5488119), 
[PyPi](https://pypi.org/project/entrezpy/), [Official document](https://entrezpy.readthedocs.io/en/master/)). 



------

[[1\]](#_ftnref1) There are 39 databases in Entrez databases and pubmed is one of it. See details (https://www.ncbi.nlm.nih.gov/books/NBK3837/)
