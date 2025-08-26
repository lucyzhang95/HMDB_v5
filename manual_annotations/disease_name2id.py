"""Manual disease name to ID mappings for HMDB data."""

# manual disease name to ontology ID mappings
MANUAL_DISEASE_MAPPINGS = {
    "early preeclampsia": "MONDO:0005081",  # preeclampsia
    "late-onset preeclampsia": "MONDO:0005081",  # preeclampsia
    "perillyl alcohol administration for cancer treatment": "GO:0018457",  # perillyl-alcohol dehydrogenase (NAD+) activity
    "3-hydroxyisobutyric acid dehydrogenase deficiency": "MONDO:0009371",  # 3-hydroxyisobutyric aciduria
    "attachment loss": "UMLS:C0206114",  # periodontal attachment loss
    "periodontal probing depth": "UMLS:C1882338",  # periodontal probing
    "methylmalonic aciduria mitochondrial encephelopathy leigh-like": "UMLS:C1855119",  # methylmalonic aciduria
    "functional hypothalamic amenorrhea": "UMLS:C0341862",  # hypothalamic amenorrhea
    "peroxisomal disorders, new type, liver": "UMLS:C5568675",  # liver disease due to peroxisomal disease
    "dopamine-serotonin vesicular transport defect": "MONDO:0018130",  # brain dopamine-serotonin vesicular transport disease
    "prosthesis/missing teeth": "UMLS:C0080233",  # tooth Loss
    "small intestinal malabsorption": "UMLS:C1833057",  # malabsorption (small intestine)
    "refractory localization-related epilepsy": "UMLS:C0472349",  # localization-related symptomatic epilepsy
    "d-glyceric acidura": "MONDO:0009070",  # d-glyceric aciduria
    "methyl formate exposure": "ECTO:9001470",  # exposure to methyl formate
    "formic acid intoxication": "ECTO:9000376",  # exposure to formic acid
    "idiopathic oro-facial pain": "MONDO:0018362",  # persistent idiopathic facial pain
    "serine deficiency syndrome, infantile": "MONDO:0035004",  # serine biosynthesis pathway deficiency, infantile/juvenile form
    "hepatic and biliary malignancies": "MONDO:0002514",  # hepatobiliary neoplasm
    "acute seizures": "UMLS:C0036572",  # seizures
    "methamphetamine (map) psychosis": "MONDO:0005465",  # methamphetamine-induced psychosis
    "nicotinamide adenine dinucleotide deficiency": "UMLS:C1283629",  # deficiency of NAD+ nucleosidase
    "homozygous sickle cell disease": "MONDO:0011382",  # sickle cell disease
    "prepartum depression": "UMLS:C0011570",  # mental depressionm -> only postpartum depression exists
    "nucleotide depletion syndrome": "MONDO:0018158",  # mitochondrial DNA depletion syndrome
    "terminal aldosterone biosynthesis defects": "MONDO:0018541",  # familial hypoaldosteronism (Aldosterone synthase deficiency is a rare inherited defect of the final step of aldosterone biosynthesis)
    "glutaryl-coa dehydrogenase deficiency (gdhd)": "MONDO:0009281",  # glutaryl-CoA dehydrogenase deficiency
    "cancer with metastatic bone disease": "UMLS:C5444038",  # metastatic bone disease
    "neuroinfection": "UMLS:C0870953",  # neuroinfections
    "dermal fibroproliferative disorder": "UMLS:C1304434",  # dermal elastolysis
    "d-lactic acidosis and short bowel syndrome": "MONDO:0015183",  # short bowel syndrome
    "tert-amyl-methyl ether exposed": "ECTO:9000278",  # exposure to ether
    "vessel occlusion": "MONDO:0020673",  # arterial occlusion
    "cresol poisoning ibs": "ECTO:9001047",  # exposure to cresol
    "dimethyl sulfide poisoning": "UMLS:C2062726",  # poisoning by sulfides
    "quetiapine poisoning": "ECTO:9000356",  # exposure to quetiapine
}


def get_manual_disease_mappings() -> dict:
    """Get manual disease name to ontology ID mappings."""
    return MANUAL_DISEASE_MAPPINGS.copy()


def format_manual_disease_mappings(disease_names: list[str]) -> dict:
    """Format manual disease mappings for given disease names."""
    manual_mappings = get_manual_disease_mappings()
    return {
        name: {"id": manual_mappings[name], "mapping_tool": "manual"}
        for name in set(disease_names)
        if name in manual_mappings
    }
