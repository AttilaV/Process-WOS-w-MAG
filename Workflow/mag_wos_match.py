
import sys
from tqdm.auto import tqdm
import pickle
import numpy as np
import dbgz
ar = np.array

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

def arange_authors(last_name, first_name):
    # First name field: may contains initials or full names. 
    # First name field:Initials may divided by space, comma, or nothing.
    surname = last_name.upper().strip()
    if len(surname) == 0:
        surname = ""
    firstname = first_name.upper().strip() 
    if len(firstname) == 0: 
        init1 = ""
        init2 = ""
    if len(firstname) == 1: 
        init1 = firstname
        init2 = ""
    if len(firstname) > 1: 
        init1 = firstname[0]
        init2 = firstname[1]
    return(surname, init1, init2)



#
# I/O files
#

if "snakemake" in sys.modules:
    Papers_file = snakemake.input["Papers_file"]
    Journals_file = snakemake.input["Journals_file"]
    Authors_file = snakemake.input["Authors_file"]
    PaperAuthorAffiliations_file = snakemake.input["PaperAuthorAffiliations_file"]
    Affiliations_file = snakemake.input["Affiliations_file"]
    wos_data_file= snakemake.input["wos_data_file"]
    wos_magauthors_file = snakemake.output["wos_magauthors_file"]
    wos_mag_matching_file = snakemake.output["wos_mag_matching_file"]
    wos_country_file = snakemake.output["wos_country_file"]
else:
    # MAG
    mag_data_path = '/gpfs/sciencegenome/MAG/mag-2021-08-02/mag/'
    Papers_file = mag_data_path + "Papers.txt"
    Journals_file = mag_data_path + "Journals.txt"
    Authors_file = mag_data_path + "Authors.txt"
    PaperAuthorAffiliations_file = mag_data_path + "PaperAuthorAffiliations.txt"
    Affiliations_file = mag_data_path + "Affiliations.txt"
    # WOS
    wos_data_file = "/gpfs/sciencegenome/WOS/WoS_2022_All.dbgz"
    # Outputs
    outpath = "/gpfs/sciencegenome/WOS2MAG/WOS_MAG_match"
    wos_magauthors_file = outpath + "WOSpaper2MAGauthor.pickle"
    wos_mag_matching_file = outpath + "WOS2MAG.pickle"
    wos_country_file = outpath + "WOSpaper2MAGcountry.pickle"



#
# MAG data
#

# Papers ~ 38 mins
print("Mag papers")
MAG2year= {}
MAG2doi = {}
MAG2title = {}
MAG2vol = {}
MAG2first_page = {}
count = 0
with open(Papers_file) as f:
    pbar = tqdm(total= 265150698)
    for line in f:
        count+=1
        #pbar.update(1)
        line = line.split("\t")
        paperID, doi, doctype, title, year, journalid, vol, issue, first_page = int(line[0]), line[2], line[3], line[4], line[7], line[11], line[14], line[15], line[16]
        title = title.upper().strip()
        title = ''.join(ch for ch in title if ch.isalnum())
        MAG2doi[paperID] = doi
        MAG2title[paperID] = title
        MAG2vol[paperID] = vol
        MAG2first_page[paperID] = first_page
        MAG2year[paperID] = year

# Authors ~ 7
print("MAG author names")
MAGauthor2name = {}
count = 0
with open(Authors_file) as f:
    pbar = tqdm(total= 274862239)
    for line in f:
        count+=1
        #pbar.update(1)
        line = line.split("\t")
        authorID, name = int(line[0]), line[2]
        MAGauthor2name[authorID] = name

# Papers 2 authors/affiliation ~ 22
print("MAG paper to MAG authors")
MAG2author_unstructured = {}
count = 0
with open(PaperAuthorAffiliations_file) as f:
    pbar = tqdm(total= 717851883)
    for line in f:
        count+=1
        pbar.update(1)
        line = line.split("\t")
        paperID, authorID, affiliationID, sequence = int(line[0]), int(line[1]), line[2], int(line[3])
        if affiliationID:
            affiliationID = int(affiliationID)
        else:
            affiliationID = 0
        if paperID in MAG2author_unstructured:
            MAG2author_unstructured[paperID].append((authorID, sequence, affiliationID))
        else:
            MAG2author_unstructured[paperID] = [(authorID, sequence, affiliationID)]

# Affiliations
print("MAG affiliations")
affiliationID2country = {}
count = 0
with open(Affiliations_file) as f:
    pbar = tqdm(total= 27043)
    for line in f:
        count+=1
        #pbar.update(1)
        line = line.split("\t")
        affiliationID, country = int(line[0]), line[10]
        affiliationID2country[affiliationID] = country

print("MAG author dictionary")
MAG2author = {}
for magid in tqdm(MAG2author_unstructured):
    authorinfo = MAG2author_unstructured[magid]
    sequence = [i[1] for i in authorinfo]
    authors = [i[0] for i in authorinfo]
    affiliationIDs = [i[2] for i in authorinfo]
    ordered_authors = [authors[i] for i in np.argsort(sequence)]
    ordered_countries = [affiliationID2country[affiliationIDs[i]] if affiliationIDs[i] != 0  else "XX" for i in np.argsort(sequence)]
    author_and_country = [ordered_authors, ordered_countries]
    MAG2author[magid] = author_and_country


# Bibliographic IDs for MAG
print("MAG bibliographic id")
MAG2biblio = {}
for paperID in tqdm(MAG2year):
    if paperID in MAG2author:
        year = MAG2year[paperID]
        vol = MAG2vol[paperID]
        page = MAG2first_page[paperID]
        firstau = MAG2author[paperID][0][0]
        if year and vol and page and firstau:
            fa_name = MAGauthor2name[firstau]
            if fa_name:
                if " " in fa_name:
                    fa_name = fa_name.split(" ")
                    family = fa_name[-1]
                    initial = fa_name[0][0]
                    MAG2biblio[paperID] = family.upper() + "_" + initial.upper() + "(" + str(year) + ")" + vol.strip() + ":" + page.strip()







        #
        # WOS data presets
        # 

print("WOS data")
WOS2biblio = {}
WOS2doi = {}
WOS2title = {}
ALLWOSIDS = set()
with dbgz.DBGZReader(wos_data_file) as fd:
    for WOSentry in tqdm(fd.entries, total=fd.entriesCount):
        paperID = WOSentry["UID"]
        # disaggregatums
        sdata = getentry(WOSentry, ["data", "static_data"])
        summary = getentry(sdata, ["summary"])
        pubinfo = getentry(summary, ["pub_info"])
        ddata = getentry(WOSentry, ["data","dynamic_data"])
        # year
        if pubinfo:
            Year = getentry(pubinfo, ["@pubyear"])
        # doi
        doi = ""
        identifiers = getentry(ddata, ["cluster_related","identifiers","identifier"])
        if identifiers:
            if isinstance(identifiers, list) == False and isinstance(identifiers, tuple) == False:
                identifiers = [identifiers]
            identifiers = {entry["@type"]: entry["@value"] for entry in identifiers}
            if("xref_doi" in identifiers):
                doi = identifiers["xref_doi"]
            if("doi" in identifiers):
                doi = identifiers["doi"]
        # first author
        Authors = getentry(summary, ["names", "name"])
        author_seq = []
        Authors_full = []
        if Authors:
            if isinstance(Authors, list) == False and isinstance(Authors, tuple) == False:
                Authors = [Authors]
            for au in Authors:
                if "last_name" in au and "first_name" in au:
                    lastname, firstname = remove_Nones(au["last_name"]), remove_Nones(au["first_name"])
                    Authors_full.append(lastname + "_" + firstname)
                    author_seq.append(au["@seq_no"])
        Authors_full = [Authors_full[i] for i in np.argsort(author_seq)] # use seq_n to arrange names
        if Authors_full:
            firstauthor = Authors_full[0]
            firstauthor = firstauthor.split("_")
            if firstauthor[0] and firstauthor[1]:
                last, i1, i2 = arange_authors(firstauthor[0], firstauthor[1])
                firstauthor = last + "_" + i1
            else:
                firstauthor = ""
            # volume
            volume = getentry(pubinfo, ["@vol"])
            # page number
            page = getentry(pubinfo, ["page", "@begin"])
            # paper title
            title = getentry(summary, ["titles", "title"])[-1]["#text"]
            title = title.upper().strip()
            title = ''.join(ch for ch in title if ch.isalnum())
            if title:
                WOS2title[paperID] = title
            if doi:
                WOS2doi[paperID] = doi
            if firstauthor and Year and volume and page:
                WOS2biblio[paperID] = firstauthor + "(" + Year + ")" + volume.strip() + ":" + page.strip()
            ALLWOSIDS.add(paperID)




                #
                # Matching
                #

"""
Matching algorithm
1) If doi match is found, and only one ID is associated with the doi in both dataset, use that.
2) If periodical information (bibliographic info) match is found, and only one ID is associated with that bibliographic ID in both dataset, use that.
3) If title match is found, and only one ID is associated with the title in both dataset, use that.
"""

print("title2MAG")
title2MAG = {} #~ MAG 4.7% is ambiguous
for paperid in MAG2title:
    title = MAG2title[paperid].upper().strip()
    if title not in title2MAG:
        title2MAG[title] = [paperid]
    else:
        title2MAG[title].append(paperid)

print("doi2MAG")
doi2MAG = {} #~ MAG 0.7% is ambiguous in WOS its 0.2%
for paperid in MAG2doi:
    doi = MAG2doi[paperid].upper().strip()
    if doi not in doi2MAG:
        doi2MAG[doi] = [paperid]
    else:
        doi2MAG[doi].append(paperid)

print("biblio2MAG")
biblio2MAG = {} #~ MAG 1.7% is ambiguous in WOS its 1%
for paperid in MAG2biblio:
    biblio = MAG2biblio[paperid]
    if biblio not in biblio2MAG:
        biblio2MAG[biblio] = [paperid]
    else:
        biblio2MAG[biblio].append(paperid)

print("title2WOS")
title2WOS= {} 
for paperid in WOS2title:
    title = WOS2title[paperid].upper().strip()
    if title not in title2WOS:
        title2WOS[title] = [paperid]
    else:
        title2WOS[title].append(paperid)

print("doi2WOS")
doi2WOS = {} 
for paperid in WOS2doi:
    doi = WOS2doi[paperid].upper().strip()
    if doi not in doi2WOS:
        doi2WOS[doi] = [paperid]
    else:
        doi2WOS[doi].append(paperid)

print("biblio2WOS")
biblio2WOS = {}
for paperid in WOS2biblio:
    biblio = WOS2biblio[paperid]
    if biblio not in biblio2WOS:
        biblio2WOS[biblio] = [paperid]
    else:
        biblio2WOS[biblio].append(paperid)

#
# Matching
#

print("Match MAG and WOS")
WOS2MAG = {}
dois = 0
biblios = 0
titles = 0
total = 0
for wosid in ALLWOSIDS:
    doi = "NA"
    biblio = "NA"
    title = "NA"
    # doi
    if wosid in WOS2doi:
        wos_doi = WOS2doi[wosid].upper().strip()
        if wos_doi in doi2MAG:
            if len(doi2MAG[wos_doi]) == 1 and len(doi2WOS[wos_doi]) == 1:
                doi = wos_doi
    # title
    if wosid in WOS2title:
        wos_title = WOS2title[wosid]
        if wos_title in title2MAG:
            if len(title2MAG[wos_title]) == 1 and len(title2WOS[wos_title]) == 1:
                title = wos_title
    # biblio
    if wosid in WOS2biblio:
        wos_biblio = WOS2biblio[wosid]
        if wos_biblio in biblio2MAG:
            if len(biblio2MAG[wos_biblio]) == 1 and len(biblio2WOS[wos_biblio]) == 1:
                biblio = wos_biblio
    # matching statistics
    if doi != "NA": dois += 1
    if biblio != "NA": biblios += 1
    if title != "NA": titles += 1
    if doi != "NA" or biblio != "NA" or title != "NA": total += 1
    #
    if doi != "NA": 
        WOS2MAG[wosid] = doi2MAG[doi][0]
    else:
        if biblio != "NA": 
            WOS2MAG[wosid] = biblio2MAG[biblio][0]
        else:
            if title != "NA":
                WOS2MAG[wosid] = title2MAG[title][0]

#
# Assign authors to WOS papers.
#


print("Assign authors to WOS papers.")
WOSpaper2MAGauthor = {}
WOSpaper2MAGaffiliation = {}
for wosid in WOS2MAG:
    magid = WOS2MAG[wosid]
    if magid in MAG2author:
        authors = MAG2author[magid][0]
        countries = MAG2author[magid][1]
        if len(authors) > 0:
            WOSpaper2MAGauthor[wosid] = authors
            if len([i for i in countries if i != "XX"]) > 0: # has any affiliation information
                WOSpaper2MAGaffiliation[wosid] = countries

#
# Save
#

with open(wos_mag_matching_file, 'wb') as f:
    pickle.dump(WOS2MAG, f)

with open(wos_magauthors_file, 'wb') as f:
    pickle.dump(WOSpaper2MAGauthor, f)

with open(wos_country_file, 'wb') as f:
    pickle.dump(WOSpaper2MAGaffiliation, f)



"""


#
# Utilities
#

# why titles mismatch
k = ar(list(WOS2MAG.keys()))
sample = np.random.choice(k, 1000, replace = False)
for wosid in sample:
    magid = WOS2MAG[wosid]
    wtit= WOS2title[wosid]
    mtit =MAG2title[magid]
    wtit = ''.join(ch for ch in wtit if ch.isalnum()).upper().strip()
    mtit = ''.join(ch for ch in mtit if ch.isalnum()).upper().strip()
    if  wtit != mtit:
        print(wtit)
        print(mtit)
        print("")
   
# why doi dont match  
k = ar(list(WOS2MAG.keys()))
sample = np.random.choice(k, 20000, replace = False)
for wosid in sample:
    magid = WOS2MAG[wosid]
    if wosid in WOS2doi and magid in MAG2doi:
        wd= WOS2doi[wosid]
        md =MAG2doi[magid]
        if  wd.upper().strip() != md.upper().strip():
            print(wd)
            print(md)
            print("")
       
"""






