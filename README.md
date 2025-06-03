# Data parser for [Human Metabolome Database (HMDB)]("https://www.hmdb.ca/")
***

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

<br><br>

## Microbe-Metabolite 
***
Last Modified: 05/28/2025
- **Total Microbe-Metabolite records: 922, 35 records fewer than the full data (957)** <br>
**~ 96.3% of the entire records** <br>
- Unique Metabolite: 154
- Unique Microbes: 349
- Unique publication: 88
- Records with publication with duplication: 728 (include duplicates for different microbe-metabolite associations with the same reference)

### Taxon (subject)
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

### Metabolite (object)
- Metabolite **description** in records (with duplications): 922
- Unique metabolite **description** in records: 154

```json
{"PUBCHEM.COMPOUND": 911, "INCHIKEY": 11}
```

### Publication (association)
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
<br><br>

## Metabolite-Disease
***
Last Modified: 06/02/2025
- **Total Metabolite-Disease records: 92,892**
- Unique metabolite: 22,569
- Unique disease names: 657
- Unique pmids: 1,189

### Metabolite (subject)
- `description`: 92,892

- `id` count
```json
{"PUBCHEM.COMPOUND":82580, "INCHIKEY":10312}
```

- `id xrefs` count
```json
{
   "HMDB":92892,
   "foodb.compound":89356,
   "INCHIKEY":82580,
   "CHEBI":70601,
   "CAS":69634,
   "chemspider":66317,
   "KEGG.COMPOUND":65692,
   "VMH":58656,
   "DRUGBANK":47304,
   "BIGG.METABOLITE":46086,
   "PDB":1936,
   "KEGG.DRUG":1
}
```

- metabolite `property` count
```json
{
   "synonym":92768,
   "description":92892,
   "molecular_weight":92842,
   "state":92166,
   "water_solubility":52177,
   "logp":44438,
   "melting_point":62739
}
```

### Disease (object)

- Mapped disease names to identifiers: 650 (total is 657)
- Unique disease name with OMIM ids: 402
- Unique diseases do not have identifiers: 255
- Disease with description: 73,925

- BT mapped: 135 (120 no hit)
```ruby
92 input query terms found dup hits
120 input query terms found no hit
```

- text2term mapped using MONDO ontology: 165 (90)
```json
{"MONDO": 146, "HP": 14, "GO": 3, "ECTO": 2}
```

- text2term mapped using EFO ontology: 167 (88) → did not use this mapping
 
- UMLS mapped: 47
- Manual mapped: 36
- HMDB mapped with OMIM: 402
- Total mapped unique diseases: 650 (165 + 47 + 36 + 402)
- Removed disease names: 7 
```ruby
['thymidine treatment', 'sodium nitrate consumption', 'tetrahydrofuran exposure', 'anephric patients', 'supradiaphragmatic malignancy', 'pregnene hydroxylation deficiency', 'eucalyptol exposure']
```

- Disease `id` count:
```json
{
   "MONDO":85078,
   "UMLS":4427,
   "OMIM":1749,
   "GO":1118,
   "HP":473,
   "ECTO":47
}
```

- Disease `id xrefs` count:
```json
{
   "OMIM":66032,
   "MONDO":18203,
   "UMLS":352,
   "HP":258
}
```


## Example Output

```json
{
   "_id":"5dd15728-0cfc-41c0-a063-eda42cb4e67e",
   "association":{
      "predicate":"biolink:ChemicalEntityToDiseaseAssociation",
      "infores":"hmdb_v5",
      "publication":{
         "id":"PMID:7482520",
         "pmid":[
            7482520,
            19006102,
            19678709,
            20156336,
            21773981,
            22148915,
            23940645,
            24424155,
            25037050,
            25105552,
            27015276,
            27107423,
            27275383,
            28587349
         ],
         "type":"biolink:Publication"
      }
   },
   "object":{
      "id":"MONDO:0005575",
      "name":"colorectal cancer",
      "description":"A primary or metastatic malignant neoplasm that affects the colon or rectum. Representative examples include carcinoma, lymphoma, and sarcoma. [NCIT:C4978]",
      "type":"biolink:Disease",
      "xrefs":{
         "omim":"OMIM:114500"
      },
      "original_name":"colorectal cancer"
   },
   "subject":{
      "id":"PUBCHEM.COMPOUND:66141",
      "name":"n-acetylproline",
      "synonym":[
         "(2s)-1-acetylpyrrolidine-2-carboxylic acid",
         "acetylproline",
         "(2s)-1-acetylpyrrolidine-2-carboxylate"
      ],
      "description":"N-Acetyl-L-proline or N-Acetylproline, belongs to the class of organic compounds known as N-acyl-alpha amino acids. N-acyl-alpha amino acids are compounds containing an alpha amino acid which bears an acyl group at its terminal nitrogen atom. N-Acetylproline can also be classified as an alpha amino acid or a derivatized alpha amino acid. Technically, N-Acetylproline is a biologically available N-terminal capped form of the proteinogenic alpha amino acid L-proline. N-acetyl amino acids can be produced either via direct synthesis of specific N-acetyltransferases or via the proteolytic degradation of N-acetylated proteins by specific hydrolases. N-terminal acetylation of proteins is a widespread and highly conserved process in eukaryotes that is involved in protection and stability of proteins (PMID: 16465618).  About 85% of all human proteins and 68% of all yeast proteins are acetylated at their N-terminus (PMID: 21750686). Several proteins from prokaryotes and archaea are also modified by N-terminal acetylation. The majority of eukaryotic N-terminal-acetylation reactions occur through N-acetyltransferase enzymes or NAT‚Äôs (PMID: 30054468).  These enzymes consist of three main oligomeric complexes NatA, NatB, and NatC, which are composed of at least a unique catalytic subunit and one unique ribosomal anchor. The substrate specificities of different NAT enzymes are mainly determined by the identities of the first two N-terminal residues of the target protein. The human NatA complex co-translationally acetylates N-termini that bear a small amino acid (A, S, T, C, and occasionally V and G) (PMID: 30054468). NatA also exists in a monomeric state and can post-translationally acetylate acidic N-termini residues (D-, E-). NatB and NatC acetylate N-terminal methionine with further specificity determined by the identity of the second amino acid.  N-acetylated amino acids, such as N-acetylproline can be released by an N-acylpeptide hydrolase from peptides generated by proteolytic degradation (PMID: 16465618). In addition to the NAT enzymes and protein-based acetylation, N-acetylation of free proline can also occur.  Many N-acetylamino acids, including N-acetylproline are classified as uremic toxins if present in high abundance in the serum or plasma (PMID: 26317986; PMID: 20613759). Uremic toxins are a diverse group of endogenously produced molecules that, if not properly cleared or eliminated by the kidneys, can cause kidney damage, cardiovascular disease and neurological deficits (PMID: 18287557).",
      "chemical_formula":"C7H11NO3",
      "molecular_weight":{
         "average_molecular_weight":157.1671,
         "monoisotopic_molecular_weight":157.073893223
      },
      "type":"biolink:SmallMolecule",
      "xrefs":{
         "inchikey":"INCHIKEY:GNMSLDIYJOSUSW-LURJTMIESA-N",
         "drugbank":"DRUGBANK:DB03360",
         "chebi":"CHEBI:21560",
         "hmdb":"HMDB:HMDB0094701",
         "cas":"CAS:68-95-1",
         "pdb":"PDB:N7P"
      }
   }
}
```