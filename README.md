# General description
Creates WOS paper (or any other entry type) level datasets from .bgz formatted raw WOS file. Outputs are lists. Disambiguated authors for WOS is imported form MAG. Provides matched MAG and WOS paper IDs as well. Applies filtering variables on papers. Python files can be run as standalone scripts or as a snakemake process. 

# process_wos.py
Provides various paper/publication level information in lists from raw WOS file. Enables filtering of records:
- collections_filter: List of collections to include. Science Citation Index is "WOS.SCI". If empty, no filtering.
- docutypes_filter: List of document types such as "Article","Letter", "Proceedings paper". If empty, no filtering.
- pubtypes_filter: List of publication types such as "Journal". If empty, no filtering.
- language_filter: List of languages to include. If empty, no filtering.
- nauthors_filter: Minimum number of authors.
- zeroyear_filter: First year to include.
- endyear_filter: Last year to include.
- minRef_filter: Minimum number of references by the paper needed to include the paper. If a paper has too few references it is filtered out.
The outputs are the followings:
- wosid2index: Dictionary of WOS publication ids and list indices. The key is the WOS id, and the value is the index pertaining to that record. The index identifies the position of the publication record in the rest of the output (lists).  The number of records in the BGZ files is N. The indices run from 0 to N-1.
- index2issn: ISSN
- index2year: Publication year.
- index2source: Typically a journal if it is an article.
- index2doi: DOI
- index2authors: List of author names. Last name and first name(s) separated by "_".
- index2nref: Number of references.
- index2title: Title of the publication.
- index2orcid: Authors' ORCID numbers. 
- index2languages: List of language(s) of the publication.
- index2authors_no: A list of integers for each author, linking authors to addresses. If the author does not have an address, the value is 0. Follows the ordering of names in index2authors. 
- index2orgs: Organizations listed for addresses. Stored as a dictionary where the key is the address number, and the value is an organization.
- index2countries: Countries listed for addresses. Stored as a dictionary where the key is the address number, and the value is a country.
- index2orcid: List of strings. Each string contains the last name + first names + ORCID, separated by "_".
- index2rauthors: A dictionary for reprint authors and address numbers. The keys are address numbers and the values are authors on that address. 
- index2rorgs: Reprint author address organizations. Stored as a dictionary where the key is the address number, and the value is an organization.
- index2rcountries: Reprint author countries. Stored as a dictionary where the key is the address number, and the value is a country.
- index2heading: List of level 1 disciplinary categories.
- index2subheading: List of level 2 disciplinary categories. 
- index2tradSC: List of level 3 disciplinary categories. Looks like subject categories.
- index2extSC: List of level 3 disciplinary categories. Looks like subject categories.
- index2agency: List of funding agencies.
- index2grantn: List of lists of official grant numbers for agencies. Follow the agency ordering as in index2agency.
- index2refs: References are publication indices. These are internal references, meaning that a reference is only recorded here if the referenced publication has a unique WOS id and it is indexed in this process (i.e they are not filtered out). The reference extracting function (get_references) also provides a bibliometric id for all the references, if it has enough bibliometric information recorded, not just for the ones with WOS ids. These IDs are constructed as follows: First author's last name + first author's first name's first character + publication year + volume + page number. This output is not recorded currently. Older papers tend to have more references to publications without WOS ids.
- duplicates: List of duplicated wosids and their respective titles
- filters: The values of the filters used.

# mag_wos_match.py
Matching algorithm for papers between MAG and WOS. Matching is based on DOIs, titles, and periodical information. The periodical information is used to create a bibliographic ID for each paper, when its possible. This approach originated in the practice of Eugene Garfield, a founder of WOS. These IDs are constructed as follows: First author's last name + first author's first name's first character + publication year + volume + page number. The matching algorithm is the following:
1) If DOI match is found, and only one ID is associated with the DOI in both datasets, use that.
2) If periodical information (bibliographic info) match is found, and only one ID is associated with that bibliographic ID in both datasets, use that.
3) If title match is found, and only one ID is associated with the title in both datasets, use that.
- WOS2MAG: WOS paper ids matched to MAG paper ids. Dictionary.
- WOSpaper2MAGauthor: WOS paper ids as keys, and MAG author ids as values. Authors maybe repeated, if author has several affiliations in MAG data.
- WOSpaper2MAGcountry: WOS paper ids as keys, and country affiliation of authors as values. "XX" represents missing data. Authors maybe repeated, if author has several affiliations in MAG data.

