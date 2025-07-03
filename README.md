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
- **Total records with pathways that do not have ids only names: 815803**
- **Total unique records: 815,611 (192 fewer)**
- **Total records with pathway ids: 813,854**
- Unique pathway identifiers: 48,585

## Example Output
```json
{
   "_id":"1738118_participates_in_SMP0000071",
   "association":{
      "id":"RO:0000056",
      "predicate":"biolink:ChemicalToPathwayAssociation",
      "type":"participates_in",
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
      "id":"PUBCHEM.COMPOUND:1738118",
      "name":"succinic acid",
      "synonym":[
         "1,2-ethanedicarboxylic acid",
         "acide butanedioique",
         "acide succinique",
         "acidum succinicum",
         "amber acid",
         "asuccin",
         "bernsteinsaeure",
         "butandisaeure",
         "butanedionic acid",
         "dihydrofumaric acid",
         "e363",
         "ethylenesuccinic acid",
         "hooc-ch2-ch2-cooh",
         "spirit of amber",
         "1,2-ethanedicarboxylate",
         "butanedionate",
         "dihydrofumarate",
         "ethylenesuccinate",
         "succinate",
         "2-acetamido-2-deoxy-d-glucose",
         "d-glcnac",
         "n-acetyl-d-glucosamine",
         "n-acetylchitosamine",
         "n acetyl d glucosamine",
         "2 acetamido 2 deoxy d glucose",
         "2 acetamido 2 deoxyglucose",
         "2-acetamido-2-deoxyglucose",
         "acetylglucosamine",
         "1,4-butanedioate",
         "1,4-butanedioic acid",
         "katasuccin",
         "wormwood acid",
         "1,2 ethanedicarboxylic acid",
         "1,4 butanedioic acid",
         "potassium succinate",
         "succinate, ammonium",
         "butanedioic acid",
         "succinate, potassium",
         "ammonium succinate"
      ],
      "description":"Succinic acid (succinate) is a dicarboxylic acid. It is an important component of the citric acid or TCA cycle and is capable of donating electrons to the electron transfer chain. Succinate is found in all living organisms ranging from bacteria to plants to mammals. In eukaryotes, succinate is generated in the mitochondria via the tricarboxylic acid cycle (TCA). Succinate can readily be imported into the mitochondrial matrix by the n-butylmalonate- (or phenylsuccinate-) sensitive dicarboxylate carrier in exchange with inorganic phosphate or another organic acid, e. g. malate (PMID 16143825). Succinate can exit the mitochondrial matrix and function in the cytoplasm as well as the extracellular space. Succinate has multiple biological roles including roles as a metabolic intermediate and roles as a cell signalling molecule. Succinate can alter gene expression patterns, thereby modulating the epigenetic landscape or it can exhibit hormone-like signaling functions (PMID: 26971832). As such, succinate links cellular metabolism, especially ATP formation, to the regulation of cellular function. Succinate can be broken down or metabolized into fumarate by the enzyme succinate dehydrogenase (SDH), which is part of the electron transport chain involved in making ATP. Dysregulation of succinate synthesis, and therefore ATP synthesis, can happen in a number of genetic mitochondrial diseases, such as Leigh syndrome, and Melas syndrome. Succinate has been found to be associated with D-2-hydroxyglutaric aciduria, which is an inborn error of metabolism. Succinic acid has recently been identified as an oncometabolite or an endogenous, cancer causing metabolite. High levels of this organic acid can be found in tumors or biofluids surrounding tumors. Its oncogenic action appears to due to its ability to inhibit prolyl hydroxylase-containing enzymes. In many tumours, oxygen availability becomes limited (hypoxia) very quickly due to rapid cell proliferation and limited blood vessel growth. The major regulator of the response to hypoxia is the HIF transcription factor (HIF-alpha). Under normal oxygen levels, protein levels of HIF-alpha are very low due to constant degradation, mediated by a series of post-translational modification events catalyzed by the prolyl hydroxylase domain-containing enzymes PHD1, 2 and 3, (also known as EglN2, 1 and 3) that hydroxylate HIF-alpha and lead to its degradation. All three of the PHD enzymes are inhibited by succinate. In humans, urinary succinic acid is produced by Escherichia coli, Pseudomonas aeruginosa, Klebsiella pneumonia, Enterobacter, Acinetobacter, Proteus mirabilis, Citrobacter frundii, Enterococcus faecalis (PMID: 22292465). Succinic acid is also found in Actinobacillus, Anaerobiospirillum, Mannheimia, Corynebacterium and Basfia (PMID: 22292465; PMID: 18191255; PMID: 26360870).",
      "chemical_formula":"C4H6O4",
      "molecular_weight":{
         "average_molecular_weight":118.088,
         "monoisotopic_molecular_weight":118.02660868
      },
      "state":"solid",
      "water_solubility":"83.2 mg/mL",
      "logp":"-0.59",
      "melting_point":"185 - 188 °C",
      "type":"biolink:SmallMolecule",
      "xrefs":{
         "pubchem":"PUBCHEM.COMPOUND:1738118",
         "inchikey":"INCHIKEY:KDYFGRWQOYBRFD-UHFFFAOYSA-N",
         "drugbank":"DRUGBANK:DB00139",
         "chebi":"CHEBI:15741",
         "hmdb":"HMDB:HMDB0000254",
         "cas":"CAS:110-15-6",
         "kegg":"KEGG.COMPOUND:C00042",
         "chemspider":"chemspider:1078",
         "foodb":"foodb.compound:FDB001931",
         "vmh":"VMH:SUCC"
      },
      "cellular_component":[
         "endoplasmic reticulum",
         "extracellular",
         "mitochondria",
         "peroxisome"
      ],
      "biosample":[
         "blood",
         "breast milk",
         "cerebrospinal fluid (csf)",
         "feces",
         "saliva",
         "sweat",
         "urine"
      ],
      "anatomical_entity":[
         "adipose tissue",
         "brain",
         "fibroblasts",
         "kidney",
         "liver",
         "pancreas",
         "placenta",
         "prostate",
         "skeletal muscle",
         "spleen"
      ]
   }
}
```

<br><br>
***
## Protein-Pathway
Last Modified: 07/03/2025

- **Total records: **

## Example Output
```json
{
   "_id":"O00180_participates_in_SMP0000375",
   "association":{
      "id":"RO:0000056",
      "predicate":"biolink:GeneToPathwayAssociation",
      "type":"participates_in",
      "infores":"hmdb_v5",
      "publication":{
         "id":"PMID:8605869",
         "pmid":[
            8605869,
            9362344,
            9462864,
            15489334,
            16710414
         ],
         "type":"biolink:Publication"
      }
   },
   "object":{
      "id":"SMPDB:SMP0000375",
      "name":"verapamil action pathway",
      "description":"Verapamil is a phenylalkylamine calcium channel blocker (CCB) or antagonist. There are at least five different types of calcium channels in Homo sapiens: L-, N-, P/Q-, R- and T-type. CCBs target L-type calcium channels, the major channel in muscle cells that mediates contraction. Verapamil, an organic cation, is thought to primarily block L-type calcium channels in their open state by interfering with the binding of calcium ions to the extracellular opening of the channel. It is one of only two clinically used CCBs that are cardioselective. Verapamil and diltiazem and, the other cardioselective CCB, shows greater activity against cardiac calcium channels than those of the peripheral vasculature. Other CCBs, such as nifedipine and amlodipine, have little to no effect on cardiac cells (cardiac myocytes and cells of the SA and AV nodes). Due to its cardioselective properties, verapamil may be used to treat arrhythmias (e.g. atrial fibrillation) as well as hypertension. \r\r\r\rThe first part of this pathway depicts the pharmacological action of verapamil on cardiac myocytes and peripheral arterioles and coronary arteries. Verapamil decreases cardiac myocyte contractility by inhibiting the influx of calcium ions. Calcium ions entering the cell through L-type calcium channels bind to calmodulin. Calcium-bound calmodulin then binds to and activates myosin light chain kinase (MLCK). Activated MLCK catalyzes the phosphorylation of the regulatory light chain subunit of myosin, a key step in muscle contraction. Signal amplification is achieved by calcium-induced calcium release from the sarcoplasmic reticulum through ryanodine receptors. Inhibition of the initial influx of calcium decreases the contractile activity of cardiac myocytes and results in an overall decreased force of contraction by the heart. Verapamil affects smooth muscle contraction and subsequent vasoconstriction in peripheral arterioles and coronary arteries by the same mechanism.  Decreased cardiac contractility and vasodilation lower blood pressure. \r\r\r\rThe second part of this pathway illustrates the effect of calcium channel antagonism on the cardiac action potentials. Contractile activity of cardiac myocytes is elicited via action potentials mediated by a number of ion channel proteins. During rest, or diastole, cells maintain a negative membrane potential; i.e. the inside of the cell is negatively charged relative to the cells\\x8a\\x97È extracellular environment. Membrane ion pumps, such as the sodium-potassium ATPase and sodium-calcium exchanger (NCX), maintain low intracellular sodium (5 mM) and calcium (100 nM) concentrations and high intracellular potassium (140 mM) concentrations. Conversely, extracellular concentrations of sodium (140 mM) and calcium (1.8 mM) are relatively high and extracellular potassium concentrations are low (5 mM). At rest, the cardiac cell membrane is impermeable to sodium and calcium ions, but is permeable to potassium ions via inward rectifier potassium channels (I-K1), which allow an outward flow of potassium ions down their concentration gradient. The positive outflow of potassium ions aids in maintaining the negative intracellular electric potential. When cells reach a critical threshold potential, voltage-gated sodium channels (I-Na) open and the rapid influx of positive sodium ions into the cell occurs as the ions travel down their electrochemical gradient. This is known as the rapid depolarization or upstroke phase of the cardiac action potential. Sodium channels then close and rapidly activated potassium channels such as the voltage-gated transient outward delayed rectifying potassium channel (I-Kto) and the voltage-gated ultra rapid delayed rectifying potassium channel (I-Kur) open. These events make up the early repolarization phase during which potassium ions flow out of the cell and sodium ions are continually pumped out. During the next phase, known as the plateau phase, calcium L-type channels (I-CaL) open and the resulting influx of calcium ions roughly balances the outward flow of potassium channels. During the final repolarization phase, the voltage-gated rapid (I-Kr) and slow (I-Ks) delayed rectifying potassium channels open increasing the outflow of potassium ions and repolarizing the cell. The extra sodium and calcium ions that entered the cell during the action potential are extruded via sodium-potassium ATPases and NCX and intra- and extracellular ion concentrations are restored. In specialized pacemaker cells, gradual depolarization to threshold occurs via funny channels (I-f). Blocking L-type calcium channels decreases conduction and increases the refractory period. Verapamil\\x8a\\x97Ès effects on pacemaker cells enable its use as a rate-controlling agent in atrial fibrillation. \r\r",
      "type":"biolink:Pathway",
      "xrefs":{
         "smpdb":"SMPDB:SMP0000375"
      }
   },
   "subject":{
      "id":"UniProtKB:O00180",
      "name":"KCNK1",
      "full_name":"potassium channel subfamily k member 1",
      "synonym":[
         "inward rectifying potassium channel protein twik-1",
         "potassium channel kcno1"
      ],
      "description":"Ion channel that contributes to passive transmembrane potassium transport and to the regulation of the resting membrane potential in brain astrocytes, but also in kidney and in other tissues (PubMed:15820677, PubMed:21653227). Forms dimeric channels through which potassium ions pass in accordance with their electrochemical gradient. The channel is selective for K(+) ions at physiological potassium concentrations and at neutral pH, but becomes permeable to Na(+) at subphysiological K(+) levels and upon acidification of the extracellular medium (PubMed:21653227, PubMed:22431633). The homodimer has very low potassium channel activity, when expressed in heterologous systems, and can function as weakly inward rectifying potassium channel (PubMed:15820677, PubMed:21653227, PubMed:22431633, PubMed:23169818, PubMed:25001086, PubMed:8605869, PubMed:8978667). Channel activity is modulated by activation of serotonin receptors (By similarity). Heterodimeric channels containing KCNK1 and KCNK2 have much higher activity, and may represent the predominant form in astrocytes (By similarity). Heterodimeric channels containing KCNK1 and KCNK3 or KCNK9 have much higher activity (PubMed:23169818). Heterodimeric channels formed by KCNK1 and KCNK9 may contribute to halothane-sensitive currents (PubMed:23169818). Mediates outward rectifying potassium currents in dentate gyrus granule cells and contributes to the regulation of their resting membrane potential (By similarity). Contributes to the regulation of action potential firing in dentate gyrus granule cells and down-regulates their intrinsic excitability (By similarity). In astrocytes, the heterodimer formed by KCNK1 and KCNK2 is required for rapid glutamate release in response to activation of G-protein coupled receptors, such as F2R and CNR1 (By similarity). Required for normal ion and water transport in the kidney (By similarity). Contributes to the regulation of the resting membrane potential of pancreatic beta cells (By similarity). The low channel activity of homodimeric KCNK1 may be due to sumoylation (PubMed:15820677, PubMed:20498050, PubMed:23169818). The low channel activity may be due to rapid internalization from the cell membrane and retention in recycling endosomes (PubMed:19959478). Permeable to monovalent cations with ion selectivity for K(+) > Rb(+) >> NH4(+) >> Cs(+) = Na(+) = Li(+)",
      "function":"Involved in potassium channel activity",
      "specific_function":"Weakly inward rectifying potassium channel",
      "residue_num":336,
      "molecular_weight":38142.8,
      "theoretical_pi":6.34,
      "transmembrane_region":[
         [
            21,
            41
         ],
         [
            133,
            153
         ],
         [
            178,
            198
         ],
         [
            247,
            267
         ]
      ],
      "protein_seq":">Potassium channel subfamily K member 1\nMLQSLAGSSCVRLVERHRSAWCFGFLVLGYLLYLVFGAVVFSSVELPYEDLLRQELRKLK\nRRFLEEHECLSEQQLEQFLGRVLEASNYGVSVLSNASGNWNWDFTSALFFASTVLSTTGY\nGHTVPLSDGGKAFCIIYSVIGIPFTLLFLTAVVQRITVHVTRRPVLYFHIRWGFSKQVVA\nIVHAVLLGFVTVSCFFFIPAAVFSVLEDDWNFLESFYFCFISLSTIGLGDYVPGEGYNQK\nFRELYKIGITCYLLLGLIAMLVVLETFCELHELKKFRKMFYVKKDKDEDQVHIIEHDQLS\nFSSITDQAAGMKEDQKQNEPFVATQSSACVDGPANH",
      "chromosomal_location":1,
      "locus":"1q42-q43",
      "gene_seq":">1011 bp\nATGCTGCAGTCCCTGGCCGGCAGCTCGTGCGTGCGCCTGGTGGAGCGGCACCGCTCGGCC\nTGGTGCTTCGGCTTCCTGGTGCTGGGCTACTTGCTCTACCTGGTCTTCGGCGCAGTGGTC\nTTCTCCTCGGTGGAGCTGCCCTATGAGGACCTGCTGCGCCAGGAGCTGCGCAAGCTGAAG\nCGACGCTTCTTGGAGGAGCACGAGTGCCTGTCTGAGCAGCAGCTGGAGCAGTTCCTGGGC\nCGGGTGCTGGAGGCCAGCAACTACGGCGTGTCGGTGCTCAGCAACGCCTCGGGCAACTGG\nAACTGGGACTTCACCTCCGCGCTCTTCTTCGCCAGCACCGTGCTCTCCACCACAGGTTAT\nGGCCACACCGTGCCCTTGTCAGATGGAGGTAAGGCCTTCTGCATCATCTACTCCGTCATT\nGGCATTCCCTTCACCCTCCTGTTCCTGACGGCTGTGGTCCAGCGCATCACCGTGCACGTC\nACCCGCAGGCCGGTCCTCTACTTCCACATCCGCTGGGGCTTCTCCAAGCAGGTGGTGGCC\nATCGTCCATGCCGTGCTCCTTGGGTTTGTCACTGTGTCCTGCTTCTTCTTCATCCCGGCC\nGCTGTCTTCTCAGTCCTGGAGGATGACTGGAACTTCCTGGAATCCTTTTATTTTTGTTTT\nATTTCCCTGAGCACCATTGGCCTGGGGGATTATGTGCCTGGGGAAGGCTACAATCAAAAA\nTTCAGAGAGCTCTATAAGATTGGGATCACGTGTTACCTGCTACTTGGCCTTATTGCCATG\nTTGGTAGTTCTGGAAACCTTCTGTGAACTCCATGAGCTGAAAAAATTCAGAAAAATGTTC\nTATGTGAAGAAGGACAAGGACGAGGATCAGGTGCACATCATAGAGCATGACCAACTGTCC\nTTCTCCTCGATCACAGACCAGGCAGCTGGCATGAAAGAGGACCAGAAGCAAAATGAGCCT\nTTTGTGGCCACCCAGTCATCTGCCTGCGTGGATGGCCCTGCAAACCATTGA",
      "gene_description":"This gene encodes one of the members of the superfamily of potassium channel proteins containing two pore-forming P domains. The product of this gene has not been shown to be a functional channel, however, it may require other non-pore-forming proteins for activity. [provided by RefSeq, Jul 2008]",
      "type":"biolink:Protein",
      "protein_type":"unknown",
      "cellular_component":[
         "membrane",
         "multi-pass membrane protein (potential)"
      ],
      "xrefs":{
         "uniprotkb":"UniProtKB:O00180",
         "hgnc":"HGNC:HGNC:6272",
         "genbank_gene":"GBG:U33632",
         "genecard":"GENECARD:KCNK1",
         "geneatlas":"GENEATLAS:KCNK1",
         "hmdbp":"HMDBP:HMDBP10851",
         "pfam":[
            {
               "id":"PFAM:PF07885",
               "name":"ion_trans_2"
            }
         ],
         "entrezgene":"NCBIGene:3775"
      }
   }
}
```

<br><br>
***

## Protein-Biological Process
Last Modified: 07/03/2025

- **Total records: **

## Example Output
```json
{
   "_id":"P43220_involved_in_0007189",
   "association":{
      "id":"RO:0002331",
      "predicate":"biolink:MacromolecularMachineToBiologicalProcessAssociation",
      "type":"involved_in",
      "infores":"hmdb_v5",
      "publication":{
         "id":"PMID:7517895",
         "pmid":[
            7517895,
            7843404,
            8216285,
            8404634,
            8405712,
            9213353,
            14574404,
            15489334,
            18287102,
            19861722,
            20869417,
            21901419,
            22412906,
            27196125,
            28514449
         ],
         "type":"biolink:Publication"
      }
   },
   "object":{
      "id":"GO:0007189",
      "name":"adenylate cyclase-activating g-protein coupled receptor signaling pathway",
      "description":"A G protein-coupled receptor signaling pathway in which the signal is transmitted via the activation of adenylyl cyclase activity which results in an increase in the intracellular concentration of cyclic AMP (cAMP). This pathway is negatively regulated by phosphodiesterase, which cleaves cAMP and terminates the signaling.",
      "type":"biolink:BiologicalProcess",
      "xrefs":{
         "go":"GO:0007189"
      }
   },
   "subject":{
      "id":"UniProtKB:P43220",
      "name":"GLP1R",
      "full_name":"glucagon-like peptide 1 receptor",
      "synonym":[
         "glp-1 receptor",
         "glp-1-r",
         "glp-1r"
      ],
      "description":"G-protein coupled receptor for glucagon-like peptide 1 (GLP-1) (PubMed:19861722, PubMed:26308095, PubMed:27196125, PubMed:28514449, PubMed:7517895, PubMed:8216285, PubMed:8405712). Ligand binding triggers activation of a signaling cascade that leads to the activation of adenylyl cyclase and increased intracellular cAMP levels (PubMed:19861722, PubMed:26308095, PubMed:27196125, PubMed:28514449, PubMed:7517895, PubMed:8216285, PubMed:8405712). Plays a role in regulating insulin secretion in response to GLP-1 (By similarity)",
      "specific_function":"G-protein coupled receptor for glucagon-like peptide 1 (GLP-1) (PubMed:8405712, PubMed:8216285, PubMed:7517895, PubMed:19861722, PubMed:26308095, PubMed:27196125, PubMed:28514449). Ligand binding triggers activation of a signaling cascade that leads to the activation of adenylyl cyclase and increased intracellular cAMP levels (PubMed:8405712, PubMed:8216285, PubMed:7517895, PubMed:19861722, PubMed:26308095, PubMed:27196125, PubMed:28514449). Plays a role in regulating insulin secretion in response to GLP-1 (By similarity).",
      "residue_num":463,
      "molecular_weight":53025.22,
      "theoretical_pi":8.137,
      "gene_description":"This gene encodes a 7-transmembrane protein that functions as a receptor for glucagon-like peptide 1 (GLP-1) hormone, which stimulates glucose-induced insulin secretion. This receptor, which functions at the cell surface, becomes internalized in response to GLP-1 and GLP-1 analogs, and it plays an important role in the signaling cascades leading to insulin secretion. It also displays neuroprotective effects in animal models. Polymorphisms in this gene are associated with diabetes. The protein is an important drug target for the treatment of type 2 diabetes and stroke. Alternative splicing of this gene results in multiple transcript variants. [provided by RefSeq, Apr 2016]",
      "type":"biolink:Protein",
      "protein_type":"unknown",
      "xrefs":{
         "uniprotkb":"UniProtKB:P43220",
         "hmdbp":"HMDBP:HMDBP14647",
         "pdb":[
            "3c59",
            " 3c5t",
            " 3iol",
            " 4zgm",
            " 5e94",
            " 5nx2",
            " 5ott",
            " 5otu",
            " 5otv",
            " 5otw",
            " 5otx",
            " 5vew",
            " 5vex",
            " 6b3j",
            " 6gb1",
            " 6orv",
            " 6vcb",
            " 6x18",
            " 6x19",
            " 6x1a",
            " 6xox",
            " 7c2e",
            " 7lci",
            " 7lcj",
            " 7lck"
         ],
         "pfam":[
            {
               "id":"PFAM:PF00002",
               "name":"7tm_2"
            },
            {
               "id":"PFAM:PF02793",
               "name":"hrm"
            }
         ],
         "entrezgene":"NCBIGene:2740"
      }
   }
}
``` 