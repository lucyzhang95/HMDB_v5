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
Last Modified: 05/28/2025
- **Total Microbe-Metabolite records: 922, 35 records fewer than the full data (957)** <br>
**~ 96.3% of the entire records** <br>
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
- text2term manually overrides microbes: 7
- text2term manually deleted microbe: 1
- text2term final mapped unique microbes: 22
- No hit names after ETE3, Entrez, and text2term: 27
- Manually mapped taxon names: 25
- Total mapped microbes: 360 (360/387 = ~93%)
**Due to manual mapping, some taxon names were broadened from strain level to species level, resulting in 24 duplicates. After removing these, the total number of unique microbes is 360.*
- Unique microbes used for getting taxon info via biothings: 360
- NCIT mapped descriptions: 141 (54 are species, 74 are genus, 13 are taxon groups, and other ranks)
- Microbial **description** in records (with duplications): 646
- Unique microbial **description** in records: 135

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
- Metabolite **description** in records (with duplications): 922
- Unique metabolite **description** in records: 154

```json
{"PUBCHEM.COMPOUND": 911, "INCHIKEY": 11}
```

### Publication
```json
{"PMID": 645, "JournalArticle (url)": 81, "doi": 2} 
```

## Example output
```json
{
   "_id":"a29b50df-da09-41e9-8ee9-39fda90b08b6",
   "association":{
      "predicate":"biolink:OrganismTaxonToChemicalEntityAssociation",
      "infores":"hmdb_v5",
      "publication":{
	       "id":"PMID:22815244",
         "pmid":[
            22815244,
            16346691
         ],
         "type":"biolink:Publication"
      }
   },
   "object":{
      "id":"PUBCHEM.COMPOUND:14925",
      "name":"1,2,3-propanetricarboxylic acid",
      "synonym":[
         "3-carboxyglutaric acid",
         "3-carboxypentanedioic acid",
         "beta-carboxyglutaric acid",
         "carballylic acid",
         "carboxymethylsuccinic acid",
         "tricarballylate",
         "3-carboxyglutarate",
         "3-carboxypentanedioate",
         "b-carboxyglutarate",
         "b-carboxyglutaric acid",
         "beta-carboxyglutarate",
         "β-carboxyglutarate",
         "β-carboxyglutaric acid",
         "carballylate",
         "carboxymethylsuccinate",
         "tricarballylic acid",
         "1,2,3-propanetricarboxylate",
         "1,2,3-tripropanetricarboxylic acid",
         "propane 1,2,3-tricarboxylic acid",
         "tricarballylic acid, trisodium salt",
         "propane-1,2,3-tricarboxylate",
         "tricarballylic acid, sodium salt",
         "1,2,3-propanetricarboxylic acid"
      ],
      "description":"1,2,3-Propanetricarboxylic acid is found in corn. 1,2,3-Propanetricarboxylic acid is isolated from plants e.g. sugarbeet sap, sap of Acer saccharinum (maple syrup). Propane-1,2,3-tricarboxylic acid, also known as tricarballylic acid, carballylic acid, and beta-carboxyglutaric acid, is a tricarboxylic acid that has three carboxylic acid functional groups. The compound is an inhibitor of the enzyme aconitase and interferes with the Krebs cycle. 1,2,3-Propanetricarboxylic acid can be produced by Bacteroides, Butyrivibrio, Megasphaera, Wolinella and fungi Nectriaceae (PMID:22815244; PMID:16346691). It is also associated with Fumonisins. Fumonisins are fungal toxins produced by Fusarium verticilloides. Detection of this compound indicates presence of fumonisins in gastrointestinal tract. Corn intake or corn contaminated with fumonisins can lead to increased levels of tricarballylic acid (PMID:22815244).",
      "chemical_formula":"C6H8O6",
      "molecular_weight":{
         "average_molecular_weight":176.1241,
         "monoisotopic_molecular_weight":176.032087988
      },
      "state":"solid",
      "water_solubility":"3.32E+05 mg/L @ 18C (exp)",
      "logp":"-1.420 (est)",
      "melting_point":"166 °C",
      "type":"biolink:SmallMolecule",
      "xrefs":{
         "inchikey":"INCHIKEY:KQTIIICEAUMSDG-UHFFFAOYSA-N",
         "drugbank":"DRUGBANK:DB04562",
         "chebi":"CHEBI:45969",
         "hmdb":"HMDB:HMDB0031193",
         "cas":"CAS:99-14-9",
         "kegg":"KEGG.COMPOUND:C19806",
         "chemspider":"chemspider:14220",
         "foodb":"foodb.compound:FDB003213"
      }
   },
   "subject":{
      "id":"taxid:844",
      "taxid":844,
      "name":"wolinella succinogenes",
      "parent_taxid":843,
      "lineage":[
         844,
         843,
         72293,
         213849,
         3031852,
         29547,
         3379134,
         2,
         131567,
         1
      ],
      "rank":"species",
      "description":"A species of anaerobic, Gram negative, rod shaped bacteria assigned to the phylum Proteobacteria. This species is motile, oxidase positive, peroxidase and urease negative and cannot metabolize carbohydrates. W. succinogenes is a commensal organism of the gastrointestinal tract and may be an emerging pathogen. [NCIT]",
      "xrefs":{
         "ncit":"C86851"
      },
      "original_name":"wolinella succinogenes",
      "type":"biolink:OrganismTaxon"
   }
}
```

## Metabolite-Disease
Last Modified: 