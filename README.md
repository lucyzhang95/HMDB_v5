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

<br><br>
***

## Microbe-Metabolite
Last Modified: 05/28/2025
- **Total Microbe-Metabolite records: 922, 35 records fewer than the full data (957)** <br>
**~ 96.3% of the entire records** <br>
- Unique Metabolite: 154
- Unique Microbes: 349
- Unique publication: 88
- Records with publication with duplication: 728 (include duplicates for different microbe-metabolite associations with the same reference)

<details>
<summary>Details</summary>

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
</details>

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
***

## Metabolite-Disease
Last Modified: 06/02/2025
- **Total Metabolite-Disease records: 92,892**
- **Total unique Metabolite-Disease records: 27,575** <br>
**Many records have the same subject and object but different publication references*
- Unique metabolite: 22,569
- Unique disease names: 657
- Unique pmids: 1,189

<details>
<summary>Details</summary>

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
</details>

## Example Output

```json
{
   "_id":"5dd15728-0cfc-41c0-a063-eda42cb4e67e",
   "association":{
      "predicate":"biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation",
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
<br><br>
***

## Metabolite-Protein
Last Modified: 06/30/2025
- **Total records with potential duplicates: 863,759**
- **Total unique records:  863,333 (426 fewer)**

## Example Output
```json
{
   "_id":"6baa1455-9ece-4ee3-8ab3-d142e3fa55d3",
   "association":{
      "predicate":"biolink:ChemicalGeneInteractionAssociation",
      "infores":"hmdb_v5"
   },
   "object":{
      "id":"UniProtKB:Q9NPD5",
      "name":"SLCO1B3",
      "full_name":"Solute carrier organic anion transporter family member 1B3",
      "description":"Mediates the Na(+)-independent uptake of organic anions (PubMed:10779507, PubMed:15159445, PubMed:17412826). Shows broad substrate specificity, can transport both organic anions such as bile acid taurocholate (cholyltaurine) and conjugated steroids (17-beta-glucuronosyl estradiol, dehydroepiandrosterone sulfate (DHEAS), and estrone 3-sulfate), as well as eicosanoid leukotriene C4, prostaglandin E2 and L-thyroxine (T4) (PubMed:10779507, PubMed:11159893, PubMed:12568656, PubMed:15159445, PubMed:17412826, PubMed:19129463). Hydrogencarbonate/HCO3(-) acts as the probable counteranion that exchanges for organic anions (PubMed:19129463). Shows a pH-sensitive substrate specificity towards sulfated steroids, taurocholate and T4 which may be ascribed to the protonation state of the binding site and leads to a stimulation of substrate transport in an acidic microenvironment (PubMed:19129463). Involved in the clearance of bile acids and organic anions from the liver (PubMed:22232210). Can take up bilirubin glucuronides from plasma into the liver, contributing to the detoxification-enhancing liver-blood shuttling loop (PubMed:22232210). Transports coproporphyrin I and III, by-products of heme synthesis, and may be involved in their hepatic disposition (PubMed:26383540). May contribute to regulate the transport of organic compounds in testes across the blood-testis-barrier (Probable). Can transport HMG-CoA reductase inhibitors (also known as statins) such as pitavastatin, a clinically important class of hypolipidemic drugs (PubMed:15159445). May play an important role in plasma and tissue distribution of the structurally diverse chemotherapeutic drugs methotrexate and paclitaxel (PubMed:23243220). May also transport antihypertension agents, such as the angiotensin-converting enzyme (ACE) inhibitor prodrug enalapril, and the highly selective angiotensin II AT1-receptor antagonist valsartan, in the liver (PubMed:16624871, PubMed:16627748)",
      "type":"biolink:Protein",
      "protein_type":"transporter",
      "xrefs":{
         "hmdbp":"HMDBP:HMDBP00446",
         "uniprotkb":"UniProtKB:Q9NPD5"
      }
   },
   "subject":{
      "id":"PUBCHEM.COMPOUND:3080612",
      "name":"12-ketodeoxycholic acid",
      "synonym":[
         "12-ketodeoxycholate",
         "(3a,5b)-3-hydroxy-12-oxo-cholan-24-oate",
         "(3a,5b)-3-hydroxy-12-oxo-cholan-24-oic acid",
         "(3alpha,5beta)-3-hydroxy-12-oxo-cholan-24-oate",
         "(3alpha,5beta)-3-hydroxy-12-oxo-cholan-24-oic acid",
         "12-dehydrodeoxycholate",
         "12-dehydrodeoxycholic acid",
         "12-ketolithocholate",
         "12-ketolithocholic acid",
         "12-oxo-3a-hydroxy-5b-cholanate",
         "12-oxo-3a-hydroxy-5b-cholanic acid",
         "12-oxolithocholate",
         "12-oxolithocholic acid",
         "3a-hydroxy-12-oxo-5b-cholan-24-oate",
         "3a-hydroxy-12-oxo-5b-cholan-24-oic acid",
         "3a-hydroxy-12-oxo-5b-cholanate",
         "3a-hydroxy-12-oxo-5b-cholanic acid",
         "3a-hydroxy-12-oxo-5b-cholanoate",
         "3a-hydroxy-12-oxo-5b-cholanoic acid",
         "3-hydroxy-12-ketocholanoic acid",
         "12-ketolithocholic acid, (3beta,5beta)-isomer",
         "3 alpha-hydroxy-12-keto-5 beta-cholanoic acid",
         "(4r)-4-[(1s,2s,5r,7r,10r,11s,14r,15r)-5-hydroxy-2,15-dimethyl-16-oxotetracyclo[8.7.0.0²,⁷.0¹¹,¹⁵]heptadecan-14-yl]pentanoate"
      ],
      "description":"12-Ketodeoxycholic acid is a bile acid. Bile acids are steroid acids found predominantly in the bile of mammals. The distinction between different bile acids is minute, depending only on the presence or absence of hydroxyl groups on positions 3, 7, and 12. Bile acids are physiological detergents that facilitate excretion, absorption, and transport of fats and sterols in the intestine and liver. Bile acids are also steroidal amphipathic molecules derived from the catabolism of cholesterol. They modulate bile flow and lipid secretion, are essential for the absorption of dietary fats and vitamins, and have been implicated in the regulation of all the key enzymes involved in cholesterol homeostasis. Bile acids recirculate through the liver, bile ducts, small intestine and portal vein to form an enterohepatic circuit. They exist as anions at physiological pH and, consequently, require a carrier for transport across the membranes of the enterohepatic tissues. The unique detergent properties of bile acids are essential for the digestion and intestinal absorption of hydrophobic nutrients. Bile acids have potent toxic properties (e.g. membrane disruption) and there are a plethora of mechanisms to limit their accumulation in blood and tissues (PMID: 11316487, 16037564, 12576301, 11907135).",
      "chemical_formula":"C24H38O4",
      "molecular_weight":{
         "average_molecular_weight":390.5561,
         "monoisotopic_molecular_weight":390.277009704
      },
      "state":"solid",
      "type":"biolink:SmallMolecule",
      "xrefs":{
         "pubchem":"PUBCHEM.COMPOUND:3080612",
         "inchikey":"INCHIKEY:CVNYHSDFZXHMMJ-VPUMZWJWSA-N",
         "chebi":"CHEBI:803536",
         "hmdb":"HMDB:HMDB0000328",
         "cas":"CAS:5130-29-0",
         "chemspider":"chemspider:2338365",
         "foodb":"foodb.compound:FDB021953"
      }
   }
}
```
<br><br>
***

## Metabolite-Pathway
Last Modified: 07/01/2025
- **Total records with potential duplicates: **

## Example Output
```json
{
   "_id":"92135_participates_in_SMP0000071",
   "association":{
      "predicate":"biolink:ChemicalToPathwayAssociation",
      "infores":"hmdb_v5"
   },
   "object":{
      "id":"SMPDB:SMP0000071",
      "name":"ketone body metabolism",
      "description":"Ketone bodies are consisted of acetone, beta-hydroxybutyrate and acetoacetate. In liver cells' mitochondria, acetyl-CoA can synthesize acetoacetate and beta-hydroxybutyrate; and spontaneous decarboxylation of acetoacetate will form acetone. Metabolism of ketone body (also known as ketogenesis) contains several reactions. Acetoacetic acid (acetoacetate) will be catalyzed to form acetoacetyl-CoA irreversibly by 3-oxoacid CoA-transferase 1 that also coupled with interconversion of succinyl-CoA and succinic acid. Acetoacetic acid can also be catalyzed by mitochondrial D-beta-hydroxybutyrate dehydrogenase to form (R)-3-Hydroxybutyric acid with NADH. Ketogenesis occurs mostly during fasting and starvation.  Stored fatty acids will be broken down and mobilized to produce large amount of acetyl-CoA for ketogenesis in liver, which can reduce the demand of glucose for other tissues. Acetone cannot be converted back to acetyl-CoA; therefore, they are either breathed out through the lungs or excreted in urine. ",
      "type":"biolink:Pathway",
      "xrefs":{
         "smpdb":"SMPDB:SMP0000071",
         "kegg":"KEGG:map00072"
      }
   },
   "subject":{
      "id":"PUBCHEM.COMPOUND:92135",
      "name":"3-hydroxybutyric acid",
      "synonym":[
         "(r)-(-)-beta-hydroxybutyric acid",
         "(r)-3-hydroxybutanoic acid",
         "3-d-hydroxybutyric acid",
         "d-3-hydroxybutyric acid",
         "d-beta-hydroxybutyric acid",
         "(r)-(-)-b-hydroxybutyrate",
         "(r)-(-)-b-hydroxybutyric acid",
         "(r)-(-)-beta-hydroxybutyrate",
         "(r)-(-)-β-hydroxybutyrate",
         "(r)-(-)-β-hydroxybutyric acid",
         "(r)-3-hydroxybutanoate",
         "3-d-hydroxybutyrate",
         "d-3-hydroxybutyrate",
         "d-b-hydroxybutyrate",
         "d-b-hydroxybutyric acid",
         "d-beta-hydroxybutyrate",
         "d-β-hydroxybutyrate",
         "d-β-hydroxybutyric acid",
         "(r)-3-hydroxybutyrate",
         "3-delta-hydroxybutyrate",
         "3-delta-hydroxybutyric acid",
         "bhib",
         "d-(-)-3-hydroxybutyrate",
         "delta-(-)-3-hydroxybutyrate",
         "delta-3-hydroxybutyrate",
         "delta-3-hydroxybutyric acid",
         "delta-beta-hydroxybutyrate",
         "3r-hydroxy-butanoate",
         "(r)-3-hydroxybutyric acid",
         "(-)-3-hydroxy-n-butyric acid",
         "(-)-3-hydroxybutyric acid",
         "(3r)-3-hydroxybutanoic acid",
         "(3r)-3-hydroxybutyric acid",
         "(3r)-hydroxybutyrate",
         "(r)-(-)-3-hydroxybutyric acid",
         "(r)-beta-hydroxybutanoic acid",
         "(r)-beta-hydroxybutyric acid",
         "(r)-β-hydroxybutanoic acid",
         "(r)-β-hydroxybutyric acid",
         "3-hydroxy-n-butyric acid",
         "3-hydroxybutanoic acid",
         "3-hydroxybutyric acid",
         "3r-hydroxybutanoic acid",
         "d-(-)-3-hydroxybutanoic acid",
         "d-(-)-3-hydroxybutyric acid",
         "d-(-)-beta-hydroxybutyric acid",
         "d-(-)-β-hydroxybutyric acid",
         "beta-hydroxy-n-butyric acid",
         "beta-hydroxybutanoic acid",
         "beta-hydroxybutyric acid",
         "β-hydroxy-n-butyric acid",
         "β-hydroxybutanoic acid",
         "β-hydroxybutyric acid"
      ],
      "description":"3-Hydroxybutyric acid (CAS: 300-85-6), also known as beta-hydroxybutanoic acid, is a typical partial-degradation product of branched-chain amino acids (primarily valine) released from muscle for hepatic and renal gluconeogenesis. This acid is metabolized by 3-hydroxybutyrate dehydrogenase (catalyzes the oxidation of 3-hydroxybutyrate to form acetoacetate, using NAD+ as an electron acceptor). The enzyme functions in nervous tissues and muscles, enabling the use of circulating hydroxybutyrate as a fuel. In the liver mitochondrial matrix, the enzyme can also catalyze the reverse reaction, a step in ketogenesis. 3-Hydroxybutyric acid is a chiral compound having two enantiomers, D-3-hydroxybutyric acid and L-3-hydroxybutyric acid, and is a ketone body. Like the other ketone bodies (acetoacetate and acetone), levels of 3-hydroxybutyrate in blood and urine are raised in ketosis. In humans, 3-hydroxybutyrate is synthesized in the liver from acetyl-CoA and can be used as an energy source by the brain when blood glucose is low. Blood levels of 3-hydroxybutyric acid levels may be monitored in diabetic patients to look for diabetic ketoacidosis. Persistent mild hyperketonemia is a common finding in newborns. Ketone bodies serve as an indispensable source of energy for extrahepatic tissues, especially the brain and lung of developing mammals. Another important function of ketone bodies is to provide acetoacetyl-CoA and acetyl-CoA for the synthesis of cholesterol, fatty acids, and complex lipids. During the early postnatal period, acetoacetate (AcAc) and beta-hydroxybutyrate are preferred over glucose as substrates for the synthesis of phospholipids and sphingolipids in accord with requirements for brain growth and myelination. Thus, during the first two weeks of postnatal development, when the accumulation of cholesterol and phospholipids accelerates, the proportion of ketone bodies incorporated into these lipids increases. On the other hand, an increased proportion of ketone bodies is utilized for cerebroside synthesis during the period of active myelination. In the lung, AcAc serves better than glucose as a precursor for the synthesis of lung phospholipids. The synthesized lipids, particularly dipalmitoylphosphatidylcholine, are incorporated into surfactant, and thus have a potential role in supplying adequate surfactant lipids to maintain lung function during the early days of life (PMID: 3884391). 3-Hydroxybutyric acid is found to be associated with fumarase deficiency and medium-chain acyl-CoA dehydrogenase deficiency, which are inborn errors of metabolism. 3-Hydroxybutyric acid is a metabolite of Alcaligenes and can be produced from plastic metabolization or incorporated into polymers, depending on the species (PMID: 7646009, 18615882).",
      "chemical_formula":"C4H8O3",
      "molecular_weight":{
         "average_molecular_weight":104.1045,
         "monoisotopic_molecular_weight":104.047344122
      },
      "state":"solid",
      "melting_point":"49 - 50 °C",
      "type":"biolink:SmallMolecule",
      "xrefs":{
         "pubchem":"PUBCHEM.COMPOUND:92135",
         "inchikey":"INCHIKEY:WHBMMWSBFZVSSR-GSVOUGTGSA-N",
         "chebi":"CHEBI:17066",
         "hmdb":"HMDB:HMDB0000011",
         "cas":"CAS:625-72-9",
         "kegg":"KEGG.COMPOUND:C01089",
         "chemspider":"chemspider:83181",
         "foodb":"foodb.compound:FDB021869",
         "bigg":"BIGG.METABOLITE:36784",
         "vmh":"VMH:BHB"
      }
   }
}
```