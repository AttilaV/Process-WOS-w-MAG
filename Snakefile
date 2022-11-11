from os.path import join as j

configfile: "Workflow/config.yaml"


#
# Directories
#

WOS_DATA_DIR = config["wos_data_dir"]
MAG_DATA_DIR = config["mag_data_dir"]
WOS_PROCESSED_DIR = config["wos_processed_dir"]
WOS_MAG_MATCH_DIR = config["wos_mag_match_dir"]


#
# aggregated bgz
#

wos_data_FILE = j(WOS_DATA_DIR, "WoS_2022_All.dbgz")


#
# WOS processing files
#

# procssed wos
wosid2index_FILE = j(WOS_PROCESSED_DIR, "wosid2index.pickle")
index2issn_FILE = j(WOS_PROCESSED_DIR, "index2issn.pickle")
index2year_FILE = j(WOS_PROCESSED_DIR, "index2year.pickle")
index2source_FILE = j(WOS_PROCESSED_DIR, "index2source.pickle")
index2doi_FILE = j(WOS_PROCESSED_DIR, "index2doi.pickle")
index2authors_FILE = j(WOS_PROCESSED_DIR, "index2authors.pickle")
index2nref_FILE = j(WOS_PROCESSED_DIR, "index2nref.pickle")
index2title_FILE = j(WOS_PROCESSED_DIR, "index2title.pickle")
index2languages_FILE = j(WOS_PROCESSED_DIR, "index2languages.pickle")
index2orgs_FILE = j(WOS_PROCESSED_DIR, "index2orgs.pickle")
index2authors_no_FILE = j(WOS_PROCESSED_DIR, "index2authors_no.pickle")
index2countries_FILE = j(WOS_PROCESSED_DIR, "index2countries.pickle")
index2orcid_FILE = j(WOS_PROCESSED_DIR, "index2orcid.pickle")
index2rorgs_FILE = j(WOS_PROCESSED_DIR, "index2rorgs.pickle")
index2rauthors_no_FILE = j(WOS_PROCESSED_DIR, "index2rauthors_no.pickle")
index2rcountries_FILE = j(WOS_PROCESSED_DIR, "index2rcountries.pickle")
index2heading_FILE = j(WOS_PROCESSED_DIR, "index2heading.pickle")
index2subheading_FILE = j(WOS_PROCESSED_DIR, "index2subheading.pickle")
index2tradSC_FILE = j(WOS_PROCESSED_DIR, "index2tradSC.pickle")
index2extSC_FILE = j(WOS_PROCESSED_DIR, "index2extSC.pickle")
index2agency_FILE = j(WOS_PROCESSED_DIR, "index2agency.pickle")
index2grantn_FILE = j(WOS_PROCESSED_DIR, "index2grantn.pickle")
index2refs_FILE = j(WOS_PROCESSED_DIR, "index2refs.pickle")


#
# Misc. files
#

duplicates_FILE = j(WOS_PROCESSED_DIR, "duplicates.pickle")
filters_FILE = j(WOS_PROCESSED_DIR, "filters.pickle")


#
# MAG WOS matching files
#

# MAG
Papers_FILE = j(MAG_DATA_DIR, "Papers.txt")
Journals_FILE = j(MAG_DATA_DIR, "Journals.txt")
Authors_FILE = j(MAG_DATA_DIR, "Authors.txt")
PaperAuthorAffiliations_FILE = j(MAG_DATA_DIR, "PaperAuthorAffiliations.txt")
Affiliations_FILE = j(MAG_DATA_DIR, "Affiliations.txt")

# Outputs
wos_magauthors_FILE = j(WOS_MAG_MATCH_DIR, "WOSpaper2MAGauthor.pickle")
wos_mag_matching_FILE = j(WOS_MAG_MATCH_DIR, "WOS2MAG.pickle")
wos_country_FILE = j(WOS_MAG_MATCH_DIR, "WOSpaper2MAcountry.pickle")




#
# Process WOS
#

rule process_wos:
	input:
		wos_data_file = wos_data_FILE,
	params:
		nauthors_FILTER = 1,
		zeroyear_FILTER = 1900,
		endyear_FILTER = 2021,
		minRef_FILTER = 3,
		collections_FILTER = ['WOS.SCI'],
		docutypes_FILTER = ["ARTICLE","LETTER", "PROCEEDINGS PAPER"],
		pubtypes_FILTER = [],
		language_FILTER = ["ENGLISH"],
	output:
		wosid2index_file = wosid2index_FILE,
		index2issn_file = index2issn_FILE,
		index2year_file = index2year_FILE,
		index2source_file = index2source_FILE,
		index2doi_file = index2doi_FILE,
		index2authors_file = index2authors_FILE,
		index2nref_file = index2nref_FILE,
		index2title_file = index2title_FILE,
		index2languages_file = index2languages_FILE,
		index2orgs_file = index2orgs_FILE,
		index2authors_no_file = index2authors_no_FILE,
		index2countries_file = index2countries_FILE,
		index2orcid_file = index2orcid_FILE,
		index2rorgs_file = index2rorgs_FILE,
		index2rauthors_no_file = index2rauthors_no_FILE,
		index2rcountries_file = index2rcountries_FILE,
		index2heading_file = index2heading_FILE,
		index2subheading_file = index2subheading_FILE,
		index2tradSC_file = index2tradSC_FILE,
		index2extSC_file = index2extSC_FILE,
		index2agency_file = index2agency_FILE,
		index2grantn_file = index2grantn_FILE,
		index2refs_file = index2refs_FILE,
		duplicates_file = duplicates_FILE,
		filters_file = filters_FILE,
	script:
		"Workflow/process_wos.py"



#
# Match with MAG
#

rule mag_wos_match:
	input:
		Papers_file = Papers_FILE,
		Journals_file = Journals_FILE,
		Authors_file = Authors_FILE,
		PaperAuthorAffiliations_file = PaperAuthorAffiliations_FILE,
		Affiliations_file = Affiliations_FILE,
		wos_data_file = wos_data_FILE,
	output:
		wos_magauthors_file = wos_magauthors_FILE,
		wos_mag_matching_file = wos_mag_matching_FILE,
		wos_country_file = wos_country_FILE,
	script:
		"Workflow/mag_wos_match.py"
