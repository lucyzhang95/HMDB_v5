# Data parser for [Human Metabolome Database (HMDB)]("https://www.hmdb.ca/")

HMDB 5.0 increased annotated metabolite entries to 217,920 and unannotated derivatized metabolite entries for GC-MS to 1,581,537.
The metabolite descriptions are generated, corrected, and expanded manually or semi-manually.
According to the publication, 9,445,375 highly quality predicted spectral datasets and other experimental ‘observables’ have been added to the database including:
- 312,980 predicted 1H and 13C NMR spectra
- 1,752,677 more predicted GC–MS spectra
- 1,440,324 predicted LC–MS/MS spectra
- 5,067,714 predicted retention indices
- 871,680 predicted collision cross-section values

## Stats
- Total metabolites found: 217,920

## Microbe-Metabolite 
- **Total Microbe-Metabolite records: 922**
- Unique Metabolite: 154
- Unique Microbes: 349
- Unique publication: 88
- Records with publication with duplication: 728 (include duplicates for different microbe-metabolite associations with the same reference)

### Taxon 
- Total unique microbes in HMDB v5: 387
- ETE3 mapped unique microbes: 335
- No hit names after ETE3 mapping: 52
- Entrez mapped unique microbes: 3
- No hit names after ETE3 and Entrez mapping: 49
- text2term mapped unique microbes: 23
- text2term mapped unique microbes are wrong: 8
- text2term manually override microbes: 7
- text2term manually deleted microbe: 1
- text2term final mapped unique microbes: 22
- No hit names after ETE3, Entrez, and text2term: 27
- Manually mapped taxon names: 25
- Total mapped microbes: 360 (360/387 = ~93%) → *Due to manual map, some of the taxon names get broader from strains in species, so there are 24 duplicates, the unique microbes without these duplicates are 360
- Unique microbes used for getting taxon info via biothings: 360
- NCIT mapped descriptions: 141 (54 are species, 74 are genus, 13 are taxon groups, and other ranks)
- Microbe **description** in records (with duplications): 646
- Unique microbe **description** in records: 135

```json
{
   "genus":842,
   "species":626,
   "phylum":12,
   "family":12,
   "subspecies":6,
   "no rank":3,
   "kingdom":2,
   "class":2,
   "order":1
}
```

### Metabolite
- Metabolite description in records (with duplications): 922
- Unique metabolite description in records: 154

```json
{"PUBCHEM.COMPOUND": 911, "INCHIKEY": 11}
```

### Publication
```json
{"PMID": 645, "JournalArticle (url)": 81, "doi": 2} 
```


