
import sys
import dbgz
from tqdm.auto import tqdm
import numpy as np
import pickle
from itertools import chain
import gc

#
# I/O files
#

if "snakemake" in sys.modules:
    # INPUTS
    wos_data_file = snakemake.input["wos_data_file"]
    # PARAMS
    collections_filter = snakemake.params["collections_FILTER"]
    docutypes_filter = snakemake.params["docutypes_FILTER"]
    pubtypes_filter = snakemake.params["pubtypes_FILTER"]
    language_filter = snakemake.params["language_FILTER"]
    nauthors_filter = snakemake.params["nauthors_FILTER"]
    zeroyear_filter = snakemake.params["zeroyear_FILTER"]
    endyear_filter = snakemake.params["endyear_FILTER"]
    minRef_filter = snakemake.params["minRef_FILTER"]
    # OUTPUTS
    wosid2index_file = snakemake.output["wosid2index_file"]
    index2issn_file = snakemake.output["index2issn_file"]
    index2year_file = snakemake.output["index2year_file"]
    index2source_file = snakemake.output["index2source_file"]
    index2doi_file = snakemake.output["index2doi_file"]
    index2authors_file = snakemake.output["index2authors_file"]
    index2nref_file = snakemake.output["index2nref_file"]
    index2title_file = snakemake.output["index2title_file"]
    index2languages_file = snakemake.output["index2languages_file"]
    index2orgs_file = snakemake.output["index2orgs_file"]
    index2authors_no_file = snakemake.output["index2authors_no_file"]
    index2countries_file = snakemake.output["index2countries_file"]
    index2orcid_file = snakemake.output["index2orcid_file"]
    index2rorgs_file = snakemake.output["index2rorgs_file"]
    index2rauthors_no_file = snakemake.output["index2rauthors_no_file"]
    index2rcountries_file = snakemake.output["index2rcountries_file"]
    index2heading_file = snakemake.output["index2heading_file"]
    index2subheading_file = snakemake.output["index2subheading_file"]
    index2tradSC_file = snakemake.output["index2tradSC_file"]
    index2extSC_file = snakemake.output["index2extSC_file"]
    index2agency_file = snakemake.output["index2agency_file"]
    index2grantn_file = snakemake.output["index2grantn_file"]
    index2refs_file = snakemake.output["index2refs_file"]
    duplicates_file = snakemake.output["duplicates_file"]
    filters_file = snakemake.output["filters_file"]
else:
    # INPUTS
    wos_data_file = "/gpfs/sciencegenome/WOS/WoS_2022_All.dbgz"
    # PARAMS
    collections_filter = ['WOS.SCI']
    docutypes_filter = ["ARTICLE","LETTER", "PROCEEDINGS PAPER"]
    pubtypes_filter = []
    language_filter = ["ENGLISH"]
    nauthors_filter = 1
    zeroyear_filter = 1900
    endyear_filter = 2021
    minRef_filter = 0
    # OUTPUTS
    out_folder = "/gpfs/sciencegenome/Protodisciplines/WOS_SCI_2020"
    wosid2index_file = out_folder + "/wosid2index.pickle"
    index2issn_file = out_folder + "/index2issn.pickle"
    index2year_file = out_folder + "/index2year.pickle"
    index2source_file = out_folder + "/index2source.pickle"
    index2doi_file = out_folder + "/index2doi.pickle"
    index2authors_file = out_folder + "/index2authors.pickle"
    index2nref_file = out_folder + "/index2nref.pickle"
    index2title_file = out_folder + "/index2title.pickle"
    index2languages_file = out_folder + "/index2languages.pickle"
    index2orgs_file = out_folder + "/index2orgs.pickle"
    index2authors_no_file = out_folder + "/index2authors_no.pickle"
    index2countries_file = out_folder + "/index2countries.pickle"
    index2orcid_file = out_folder + "/index2orcid.pickle"
    index2rorgs_file = out_folder + "/index2rorgs.pickle"
    index2rauthors_no_file = out_folder + "/index2rauthors_no.pickle"
    index2rcountries_file = out_folder + "/index2rcountries.pickle"
    index2heading_file = out_folder + "/index2heading.pickle"
    index2subheading_file = out_folder + "/index2subheading.pickle"
    index2tradSC_file = out_folder + "/index2tradSC.pickle"
    index2extSC_file = out_folder + "/index2extSC.pickle"
    index2agency_file = out_folder + "/index2agency.pickle"
    index2grantn_file = out_folder + "/index2grantn.pickle"
    index2refs_file = out_folder + "/index2refs.pickle"
    duplicates_file = out_folder + ["duplicates.pickle"]
    filters_file = out_folder + ["filters.pickle"]





                    #
                    # Data processed from raw WoS
                    #







def getentry(data, keys):
    for key in keys:
        if key in data:
            data = data[key]
        else:
            data = ""
            break
    return data

def remove_Nones(totest):
    if totest != None:
        return(totest)
    else:
        return("")

def ref_biblioID(author_entry, year, vol, page, maxchar):
    # Parse surnames
    author = author_entry.split(",")
    author = author[0].split(" ")
    while len(author) > 1:
        author = author[0].split(" ")
    author = author[0].split(".")
    while len(author) > 1:
        author = author[0].split(".")
    surname=[c for c in author[0] if c.isalpha() == True]
    surname = "".join(surname).upper().rstrip()
    # Parse firstnames
    author = author_entry.split(",")
    author = author[-1].split(" ")
    while len(author) > 1:
        author = author[-1].split(" ")
    author = author[-1].split(".")
    while len(author) > 1:
        author = author[-1].split(".")
    firstname = [c for c in author[0] if c.isalpha() == True]
    firstname = "".join(firstname).upper().rstrip()
    if len(surname) > 0 and len(firstname) > 0 and len(year) > 0 and len(page) > 0 and len(vol) > 0 and "Anonymous" not in surname and "ANONYMOUS" not in surname:
        ID = surname[0:maxchar] + '_' + firstname[0] + '(' + year.rstrip() + ')' + vol.rstrip() + ':' + page.rstrip()
    else:
        ID = ""
    return(ID)

def arange_wos_entries(WOSentry):
    paperID = WOSentry["UID"]
    # disaggregatums
    sdata = getentry(WOSentry, ["data", "static_data"])
    summary = getentry(sdata, ["summary"])
    pubinfo = getentry(summary, ["pub_info"])
    ddata = getentry(WOSentry, ["data","dynamic_data"])
    fullrec = getentry(sdata, ["fullrecord_metadata"])
    # misc paper meta data
    if pubinfo:
        Pubtype = getentry(pubinfo, ["@pubtype"])
        Year = getentry(pubinfo, ["@pubyear"])
    Doctype = getentry(fullrec, ["normalized_doctypes","doctype"])
    if isinstance(Doctype, list) == False and isinstance(Doctype, tuple) == False:
        Doctype = [Doctype]
    # Edition
    editions = getentry(summary, ["EWUID","edition"])
    if isinstance(editions, list) == False and isinstance(editions, tuple) == False:
        editions = [editions]
    edition = [entry['@value'] for entry in editions]
    # ISSN and doi
    issn = ""
    doi = ""
    identifiers = getentry(ddata, ["cluster_related","identifiers","identifier"])
    if identifiers:
        if isinstance(identifiers, list) == False and isinstance(identifiers, tuple) == False:
            identifiers = [identifiers]
        identifiers = {entry["@type"]: entry["@value"] for entry in identifiers}
        if("issn" in identifiers):
            issn = identifiers["issn"]
        if("xref_doi" in identifiers):
            doi = identifiers["xref_doi"]
        if("doi" in identifiers):
            doi = identifiers["doi"]
    # authors, and address identifier numbers
    Authors = getentry(summary, ["names", "name"])
    Authors_full = []
    author_seq = []
    if Authors:
        if isinstance(Authors, list) == False and isinstance(Authors, tuple) == False:
            Authors = [Authors]
        for au in Authors:
            if "last_name" in au and "first_name" in au:
                addr_no = ""
                if '@addr_no' in au:
                    addr_no = str(remove_Nones(au['@addr_no']))
                lastname, firstname = remove_Nones(au["last_name"]), remove_Nones(au["first_name"])
                Authors_full.append(lastname + "_" + firstname + "^" + addr_no)
                author_seq.append(au["@seq_no"])
    Authors_full = [Authors_full[i] for i in np.argsort(author_seq)] # use seq_n to arrange names
    # addresses, locations and organizations
    Organizations = []
    Countries = []
    address_name = getentry(fullrec, ["addresses","address_name"])
    if address_name:
        if isinstance(address_name, list) == False and isinstance(address_name, tuple) == False:
            address_name = [address_name]
        for address_info in address_name:
            address = address_info["address_spec"]
            addr_no = str(remove_Nones(address['@addr_no']))
            if 'organizations' in address:
                organization = remove_Nones(address["organizations"]["organization"])
                if type(organization) != str:
                    if isinstance(organization, list) == False and isinstance(organization, tuple) == False:
                        organization = [organization]
                    organization = organization[-1]["#text"]
                Organizations.append(remove_Nones(organization) + "^" + addr_no)
            if "country" in address:
                country = address["country"]
                country = remove_Nones(country)
                Countries.append(country + "^" + addr_no)
    # reprint_addresses
    Reprint_org = []
    Reprint_authors_full = []
    Reprint_country = []
    reprint_addresses = getentry(fullrec, ["reprint_addresses", 'address_name'])
    if reprint_addresses:
        if isinstance(reprint_addresses, list) == False and isinstance(reprint_addresses, tuple) == False:
            reprint_addresses = [reprint_addresses]
        for reprint_author in reprint_addresses:
            address = reprint_author['address_spec']
            addr_no = str(remove_Nones(address['@addr_no']))
            if 'organizations' in address:
                organization = remove_Nones(address["organizations"]["organization"])
                if type(organization) != str:
                    if isinstance(organization, list) == False and isinstance(organization, tuple) == False:
                        organization = [organization]
                    organization = organization[-1]["#text"]
                Reprint_org.append(remove_Nones(organization) + "^" + addr_no)
            if "country" in address:
                country = remove_Nones(address["country"])
                Reprint_country.append(country + "^" + addr_no)
            names = getentry(reprint_author, ["names"])
            if names:
                name = reprint_author['names']['name']
                if isinstance(name, list) == False and isinstance(name, tuple) == False:
                    name = [name]
                for au in name:
                    if "last_name" in au and "first_name" in au and '@addr_no' in au:
                        addr_no = str(au['@addr_no'])
                        lastname, firstname = remove_Nones(au["last_name"]), (au["first_name"])
                        Reprint_authors_full.append(lastname + "_" + firstname + "^" + str(addr_no))
    # ORCID
    orcids = []
    contributors = getentry(sdata, ['contributors', 'contributor'])
    if contributors:
        if isinstance(contributors, list) == False and isinstance(contributors, tuple) == False:
            contributors = [contributors]
        for au in contributors:
            contrib = au["name"]
            if "last_name" in contrib and "first_name" in contrib and '@orcid_id' in contrib:
                lastname = remove_Nones(contrib["last_name"])
                firstname = remove_Nones(contrib["first_name"])
                orcid = remove_Nones(contrib['@orcid_id'])
                orcids.append(lastname + "_" + firstname + "_" + orcid)
    # paper title
    Title = getentry(summary, ["titles", "title"])[-1]["#text"]
    # journal title
    Source = getentry(summary, ["titles", "title"])[0]["#text"]
    # no of references
    totref = int(getentry(fullrec, ["references", "@count"]))
    # languages
    Language = ""
    languages = getentry(fullrec, ["languages", "language"])
    if languages:
        if isinstance(languages, list) == False and isinstance(languages, tuple) == False:
            languages = [languages]
        for language in languages:
            if language["@type"] == "primary":
                Language = language["#text"]
    # Missing data cosmetics
    if not Year: Year = 9999
    # Outpuut
    return(
        paperID,
        int(Year),
        issn,
        int(totref),
        edition, # list
        Doctype, # list
        Pubtype,
        Authors_full, # list
        Organizations, # list
        Countries, # list
        Reprint_authors_full, # list
        Reprint_country, # list
        Reprint_org, # list
        orcids, # list
        Title,
        doi,
        Source,
        Language) # list

def get_categories(wosentry):
    subcategs, maincategs, subjtrad, subjext, = [],[],[],[]
    categorization = getentry(wosentry, ["data", "static_data", "fullrecord_metadata", "category_info"])   
    if categorization:
        maincategs = getentry(categorization, ['headings', "heading"])
        subcategs = getentry(categorization, ['subheadings', "subheading"])
        subjects = getentry(categorization, ['subjects', "subject"])
        if isinstance(subjects, list) == False and isinstance(subjects, tuple) == False:
            subjects = [subjects]
        for subject in subjects:
            subj = subject["#text"]
            ascatype = subject["@ascatype"]
            if ascatype == "traditional":
                subjtrad.append(subj)
            else:
                subjext.append(subj)
        # cosmetics
        if isinstance(subcategs, list) == False and isinstance(subcategs, tuple) == False:
            subcategs = [subcategs]
        if isinstance(maincategs, list) == False and isinstance(maincategs, tuple) == False:
            maincategs = [maincategs]
    return(subcategs, maincategs, subjtrad, subjext)

def get_funding(wosentry):
    agencies, gnums = [],[]
    fund_ack = getentry(wosentry, ["data", "static_data", "fullrecord_metadata","fund_ack"])
    if fund_ack:
        grants = getentry(fund_ack, ["grants","grant"])
        grants = remove_Nones(grants)
        if grants:
            if isinstance(grants, list) == False and isinstance(grants, tuple) == False:
                grants = [grants]
            for grant in grants:
                grantagency = ""
                gids = []
                if "grant_agency" in grant:
                    grantagency = remove_Nones(grant["grant_agency"])
                    if isinstance(grantagency, list) == False and isinstance(grantagency, tuple) == False:
                        grantagency = [grantagency]
                    if "#text" in grantagency[-1]:
                        grantagency = grantagency[-1]["#text"]
                if "grant_ids" in grant:
                    grantids = remove_Nones(grant["grant_ids"])
                    if "grant_ids" in grantids:
                        gids = grant["grant_ids"]["grant_id"]
                agencies.append(grantagency)
                gnums.append(gids)
    return(agencies, gnums)

def get_references(wosentry):
    wosid_references = []
    biblio_references = []
    references = getentry(wosentry, ["data", "static_data", "fullrecord_metadata", "references", "reference"])   
    if isinstance(references, list) == False and isinstance(references, tuple) == False:
        references = [references]
    for ref in references:
        paperwosID = ""
        biblioID = ""
        if "uid" in ref:
            WoS_ID = ref['uid']
            if "." not in WoS_ID:
                paperwosID = WoS_ID
        if 'citedAuthor' in ref and 'year' in ref and "volume" in ref and "page" in ref:
            Authorname, Year, Volume, Page = remove_Nones(ref['citedAuthor']), remove_Nones(ref['year']), remove_Nones(ref['volume']), remove_Nones(ref['page'])
            biblioID = ref_biblioID(Authorname, Year, Volume, Page, 9999)
        if len(paperwosID) > 0 or len(biblioID) > 0:
            wosid_references.append(paperwosID)
            biblio_references.append(biblioID) 
    return(wosid_references, biblio_references)


#
# Set up outputs
#

# paper meta data
wosid2index = {}
index2issn = []
index2year = []
index2source = []
index2doi = []
index2authors = []
index2nref = []
index2title = []
index2languages = []
# Addresses and orcid
index2authors_no = []
index2orgs = []
index2countries = []
index2orcid = []
# Reprint addresses
index2rauthors_no = []
index2rorgs = []
index2rcountries = []
# categorization information
index2heading = []
index2subheading = []
index2tradSC = []
index2extSC = []
# Funding 
index2agency = []
index2grantn = []
# references
index2refs = []
# misc
duplicates = []
filters = {}

#
# Read data
#

reference_lists = []

collections_filter = set([i.upper() for i in collections_filter])
docutypes_filter = set([i.upper() for i in docutypes_filter])
pubtypes_filter = set([i.upper() for i in pubtypes_filter])
language_filter = set([i.upper() for i in language_filter])
zero = set([0])
index = 0


with dbgz.DBGZReader(wos_data_file) as fd:
    for wosentry in tqdm(fd.entries, total=fd.entriesCount):
        ### access and arange meta data
        wos_meta = arange_wos_entries(wosentry)
        WOSID,Year,issn,totref,edition,Doctype, Pubtype,Authors_full, Organizations, Countries, Reprint_authors_full,Reprint_country, Reprint_org, orcids, Title,doi,Source,Language = wos_meta
        edition = [s.upper().strip() for s in edition]
        Doctype = [s.upper().strip() for s in Doctype]
        Pubtype = Pubtype.upper().strip()
        ### implement filter
        if totref >= minRef_filter and Year >= zeroyear_filter and Year <= endyear_filter and len(WOSID) > 0 and len(Authors_full) >= nauthors_filter:
            edf = False
            docf = False
            pubf = False
            lang = False
            if len(set(edition) & collections_filter) > 0 or len(collections_filter) == 0: edf = True
            if len(set(Doctype) & docutypes_filter) > 0 or len(docutypes_filter) == 0: docf = True
            if len(set([Pubtype]) & pubtypes_filter) > 0 or len(pubtypes_filter) == 0: pubf = True
            if len(set([Language]) & language_filter) > 0 or len(language_filter) == 0: lang = True
            if edf == True and docf == True and pubf == True and lang == True:
                if WOSID not in wosid2index:
                    wosid2index[WOSID] = index
                    ### arange into dictionaries
                    index2issn.append(issn)
                    index2year.append(Year)
                    index2source.append(Source.strip().upper())
                    index2doi.append(doi.strip().upper())
                    index2nref.append(totref)
                    index2title.append(Title.strip())
                    index2languages.append(Language)
                    # Authors and address numbers
                    authornames = []
                    author_no = []
                    for a in Authors_full:
                        authorinfo = a.split("^")
                        if len(authorinfo[1]) > 0:
                            author_no.append([int(i) for i in authorinfo[1].split(" ") if i.isnumeric()])
                        else:
                            author_no.append(0)
                        authornames.append(authorinfo[0].strip().upper())
                    index2authors.append(authornames)
                    index2authors_no.append(author_no)
                    # Orcids
                    index2orcid.append(orcids)
                    # Organizational affiliation and address numbers
                    orgnames = []
                    org_no = []
                    if len(Organizations) > 0:
                        for o in Organizations:
                            orginfo = o.split("^")
                            if len(orginfo[1]) > 0 and orginfo[1].isnumeric():
                                org_no.append(int(orginfo[1]))
                            else:
                                org_no.append(0)
                            orgnames.append(orginfo[0].strip().upper())
                    index2orgs.append({org_no[i]:orgnames[i] for i in range(len(orgnames))})
                    # Country affiliation and address numbers
                    cnames = []
                    c_no = []
                    if len(Countries) > 0:
                        for c in Countries:
                            cinfo = c.split("^")
                            if len(cinfo[1]) > 0 and cinfo[1].isnumeric():
                                c_no.append(int(cinfo[1]))
                            else:
                                c_no.append(0)
                            cnames.append(cinfo[0].strip().upper())
                    index2countries.append({c_no[i]:cnames[i] for i in range(len(cnames))})
                    # REPRINT Authors and address numbers
                    authornames = []
                    author_no = []
                    indexdict = {}
                    if len(Reprint_authors_full) > 1:
                        for a in Reprint_authors_full:
                            authorinfo = a.split("^")
                            if len(authorinfo) > 1:
                                author_no.append([int(i) for i in authorinfo[1].split(" ") if i.isnumeric()])
                            else:
                                author_no.append([0])
                            authornames.append(authorinfo[0].strip().upper())
                        numbers = set(list(chain(*author_no))) - zero
                        indexdict = {n:[] for n in numbers}
                        for i in range(len(author_no)):
                            for j in author_no[i]:
                                indexdict[j].append(authornames[i])
                    index2rauthors_no.append(indexdict)
                    # REPRINT Organizational affiliation and address numbers
                    orgnames = []
                    org_no = []
                    if Reprint_org and Reprint_authors_full:
                        for o in Reprint_org:
                            orginfo = o.split("^")
                            if len(orginfo[1]) > 0 and orginfo[1].isnumeric():
                                org_no.append(int(orginfo[1]))
                            else:
                                org_no.append(0)
                            orgnames.append(orginfo[0].strip().upper())
                    index2rorgs.append({org_no[i]:orgnames[i] for i in range(len(orgnames))})
                    # REPRINT Country affiliation and address numbers
                    cnames = []
                    c_no = []
                    if len(Reprint_country) > 0:
                        for c in Reprint_country:
                            cinfo = c.split("^")
                            if len(cinfo[1]) > 0 and cinfo[1].isnumeric():
                                c_no.append(int(cinfo[1]))
                            else:
                                c_no.append(0)
                            cnames.append(cinfo[0].strip().upper())
                    index2rcountries.append({c_no[i]:cnames[i] for i in range(len(cnames))})
                    ### access and arange funding information and categorizations
                    subcategs, maincategs, subjtrad, subjext = get_categories(wosentry)
                    agencies, gnums = get_funding(wosentry)
                    ### arange into dictionaries
                    index2heading.append(maincategs)
                    index2subheading.append(subcategs)
                    index2tradSC.append(subjtrad)
                    index2extSC.append(subjext)
                    index2agency.append(agencies)
                    index2grantn.append(gnums)
                    ### access references
                    wos_references, biblio_references = get_references(wosentry)
                    reference_lists.append([r for r in wos_references if r])
                    index += 1
                else:
                    duplicates.append((WOSID, Title))


#
# References
#


for references in tqdm(reference_lists):
    index2refs.append([wosid2index[wosid] for wosid in references if wosid in wosid2index])


#
# Save results
#

print("saving1")


with open(wosid2index_file, 'wb') as f:
    pickle.dump(wosid2index, f)

with open(index2issn_file, 'wb') as f:
    pickle.dump(index2issn, f)

with open(index2year_file, 'wb') as f:
    pickle.dump(index2year, f)

with open(index2source_file, 'wb') as f:
    pickle.dump(index2source, f)

with open(index2doi_file, 'wb') as f:
    pickle.dump(index2doi, f)

with open(index2authors_file, 'wb') as f:
    pickle.dump(index2authors, f)

with open(index2nref_file, 'wb') as f:
    pickle.dump(index2nref, f)

with open(index2title_file, 'wb') as f:
    pickle.dump(index2title, f)

with open(index2languages_file, 'wb') as f:
    pickle.dump(index2languages, f)


# Addresses
print("saving2")
with open(index2orgs_file, 'wb') as f:
    pickle.dump(index2orgs, f)

with open(index2authors_no_file, 'wb') as f:
    pickle.dump(index2authors_no, f)

with open(index2countries_file, 'wb') as f:
    pickle.dump(index2countries, f)

with open(index2orcid_file, 'wb') as f:
    pickle.dump(index2orcid, f)


# Reprint addresses
print("saving3")
with open(index2rorgs_file, 'wb') as f:
    pickle.dump(index2rorgs, f)

with open(index2rauthors_no_file, 'wb') as f:
    pickle.dump(index2rauthors_no, f)

with open(index2rcountries_file, 'wb') as f:
    pickle.dump(index2rcountries, f)


# Categorization
print("saving4")
with open(index2heading_file, 'wb') as f:
    pickle.dump(index2heading, f)

with open(index2subheading_file, 'wb') as f:
    pickle.dump(index2subheading, f)

with open(index2tradSC_file, 'wb') as f:
    pickle.dump(index2tradSC, f)

with open(index2extSC_file, 'wb') as f:
    pickle.dump(index2extSC, f)


# Funding 
print("saving5")
with open(index2agency_file, 'wb') as f:
    pickle.dump(index2agency, f)

with open(index2grantn_file, 'wb') as f:
    pickle.dump(index2grantn, f)


# flush memory
reference_lists = 0
wosid2index = 0
index2issn = 0
index2year = 0
index2source = 0
index2doi = 0
index2authors = 0
index2nref = 0
index2title = 0
index2languages = 0
index2authors_no = 0
index2orgs = 0
index2countries = 0
index2orcid = 0
index2rauthors_no = 0
index2rorgs = 0
index2rcountries = 0
index2heading = 0
index2subheading = 0
index2tradSC = 0
index2extSC = 0
index2agency = 0
index2grantn = 0
gc.collect()

# references
print("saving6")
with open(index2refs_file, 'wb') as f:
    pickle.dump(index2refs, f)



# misc

filters = {
    "wos_data_file" : "/gpfs/sciencegenome/WOS/WoS_2022_All.dbgz",
    "collections_filter" : ['WOS.SCI'],
    "docutypes_filter" : ["ARTICLE","LETTER", "PROCEEDINGS PAPER"],
    "pubtypes_filter" : [],
    "language_filter" : ["ENGLISH"],
    "nauthors_filter" : 1,
    "zeroyear_filter" : 1900,
    "endyear_filter" : 2021,
    "minRef_filter" : 0
}

print("saving7")
with open(filters_file, 'wb') as f:
    pickle.dump(filters, f)

with open(duplicates_file, 'wb') as f:
    pickle.dump(duplicates, f)





