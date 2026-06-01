import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
ARTICLE_JSON = WORK / "article_extraction.json"
INV_CSV = WORK / "pdf_inventory.csv"
OUT_JSON = WORK / "article_extraction_qc.json"
OUT_CSV = WORK / "article_extraction_qc.csv"


def norm(s):
    return re.sub(r"\s+", " ", str(s or "")).strip()


def clean_doi(doi):
    return str(doi or "").strip().rstrip(".,);]")


def title_from_filename(filename):
    stem = Path(str(filename).split("|")[0].strip()).stem.replace("_", " ")
    stem = re.sub(r"\s+Source\s+.*", "", stem, flags=re.I)
    stem = re.sub(r"\s+SO\s+\d{4}.*", "", stem, flags=re.I)
    return norm(stem)


def pubmed_search(term):
    params = {"db": "pubmed", "retmode": "json", "retmax": "3", "term": term, "sort": "relevance"}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=20) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    return data.get("esearchresult", {}).get("idlist", [])


def pubmed_fetch(pmid):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        + urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "xml"})
    )
    with urllib.request.urlopen(url, timeout=20) as resp:
        root = ET.fromstring(resp.read())
    art = root.find(".//PubmedArticle")
    if art is None:
        return {}
    article = art.find(".//Article")
    title = ""
    if article is not None and article.find("ArticleTitle") is not None:
        title = "".join(article.find("ArticleTitle").itertext())
    journal = article.findtext("Journal/Title", default="") if article is not None else ""
    year = ""
    if article is not None:
        year = article.findtext("Journal/JournalIssue/PubDate/Year", default="") or article.findtext(
            "Journal/JournalIssue/PubDate/MedlineDate", default=""
        )[:4]
    first = ""
    auth = article.find("AuthorList/Author") if article is not None else None
    if auth is not None:
        first = auth.findtext("LastName", default="") or auth.findtext("CollectiveName", default="")
    doi = ""
    for aid in art.findall(".//ArticleId"):
        if aid.attrib.get("IdType") == "doi":
            doi = aid.text or ""
    return {
        "PMID": str(pmid),
        "Title": norm(title),
        "First_author": norm(first),
        "Year": norm(year[:4]),
        "Journal": norm(journal),
        "DOI": clean_doi(doi),
    }


def token_set(s):
    stop = {"lung", "adenocarcinoma", "source", "prognostic", "signature", "model", "based", "with"}
    return {x for x in re.findall(r"[a-z0-9]{4,}", str(s).lower()) if x not in stop}


def sim(a, b):
    aa, bb = token_set(a), token_set(b)
    if not aa or not bb:
        return 0
    return len(aa & bb) / min(len(aa), len(bb))


def lookup_pubmed(title, doi=""):
    queries = []
    if doi:
        queries.append(f"{doi}[AID]")
    if title:
        queries.append(f'"{title[:180]}"[Title]')
        words = " ".join(re.findall(r"[A-Za-z0-9-]{4,}", title)[:10])
        if words:
            queries.append(words + "[Title]")
    for q in queries:
        try:
            for pmid in pubmed_search(q):
                meta = pubmed_fetch(pmid)
                if not title or sim(meta.get("Title", ""), title) >= 0.2:
                    return meta
        except Exception:
            continue
    return {}


def source_text(row, inv):
    names = [x.strip() for x in str(row.get("Source_files", "")).split("|")]
    chunks = []
    for name in names:
        hit = inv[inv["file_name"].eq(name)]
        if not hit.empty:
            chunks.append(Path(hit.iloc[0]["text_file"]).read_text(encoding="utf-8", errors="ignore"))
    return "\n\n".join(chunks)


def context(text, pattern, w=650):
    m = re.search(pattern, text, re.I)
    if not m:
        return ""
    return norm(text[max(0, m.start() - w) : min(len(text), m.end() + w)])


def set_formula(row, genes, coefs, formula):
    row["Gene_list"] = "; ".join(genes)
    row["Gene_number"] = len(genes)
    row["Coefficients"] = "; ".join(f"{g}:{c}" for g, c in zip(genes, coefs))
    row["Coefficients_available"] = "yes"
    row["Formula_available"] = "yes"
    row["Risk_score_formula"] = formula


def apply_manual_repairs(row, text):
    aid = row["Article_ID"]
    if aid == "A033":
        set_formula(
            row,
            ["GSEC", "FAM83A-AS1", "AL606489.1", "AL034397.3", "AC010980.2"],
            ["0.7919", "0.2422", "0.3857", "-0.9045", "0.4247"],
            "Risk Score = GSEC expression x 0.7919 + FAM83A-AS1 expression x 0.2422 + AL606489.1 expression x 0.3857 - AL034397.3 expression x 0.9045 + AC010980.2 expression x 0.4247",
        )
    elif aid == "A060":
        set_formula(
            row,
            ["AC245014.3", "LINC00941", "AL358115.1", "AC243960.3", "AC005180.1"],
            ["0.655730531649772", "0.462817768020813", "0.877653684548024", "-1.51062067434958", "-1.35918746869745"],
            "Risk score = 0.655730531649772 x AC245014.3 expression + 0.462817768020813 x LINC00941 expression + 0.877653684548024 x AL358115.1 expression - 1.51062067434958 x AC243960.3 expression - 1.35918746869745 x AC005180.1 expression",
        )
        row["Cutoff_value"] = "0.994682147"
        row["Cutoff_type"] = "median / training-derived"
    elif aid == "A062":
        set_formula(
            row,
            ["PLK1", "PLK2", "PRKDC", "NBN"],
            ["0.054", "0.199", "0.229", "0.320"],
            "Risk score = 0.054 x PLK1 expression + 0.199 x PLK2 expression + 0.229 x PRKDC expression + 0.320 x NBN expression",
        )
        row["GEO_accession"] = re.sub(r"\bGSE312\b;?\s*", "", str(row.get("GEO_accession", ""))).strip("; ")
        row["Training_cohort"] = re.sub(r"\bGSE312\b;?\s*", "", str(row.get("Training_cohort", ""))).strip("; ")
    elif aid == "A069":
        set_formula(
            row,
            ["IL1A", "P2RX1", "GSTM2", "ESYT3", "ZNF682", "KCNF1", "STK32A", "HHIPL2", "GDF10", "NDC80", "GSTA1", "BCL2L10", "CCR2"],
            ["0.1077", "-0.0146", "-0.0630", "-0.0344", "-0.1546", "0.0292", "-0.0316", "0.0186", "-0.0115", "0.0609", "-0.0298", "0.0571", "-0.2199"],
            "Risk Score = 0.1077 x IL1A - 0.0146 x P2RX1 - 0.0630 x GSTM2 - 0.0344 x ESYT3 - 0.1546 x ZNF682 + 0.0292 x KCNF1 - 0.0316 x STK32A + 0.0186 x HHIPL2 - 0.0115 x GDF10 + 0.0609 x NDC80 - 0.0298 x GSTA1 + 0.0571 x BCL2L10 - 0.2199 x CCR2",
        )
    elif aid == "A070":
        set_formula(
            row,
            ["RP11-83A24.2", "AC007128.1", "RP11-363N22.3", "HNRNPUP1", "RP11-15H20.6", "CTD-2017D11.1"],
            ["-0.018392", "0.339099", "-0.034758", "0.116564", "-0.421656", "-0.103568"],
            "Risk Score = -0.018392 x RP11-83A24.2 + 0.339099 x AC007128.1 - 0.034758 x RP11-363N22.3 + 0.116564 x HNRNPUP1 - 0.421656 x RP11-15H20.6 - 0.103568 x CTD-2017D11.1",
        )
    elif aid == "A084":
        row["PMID"] = "37197625"
        row["Title"] = "Development of a cuproptosis-related signature for prognosis prediction in lung adenocarcinoma based on WGCNA."
        row["First_author"] = "Ling"
        row["Year"] = "2023"
        row["Journal"] = "Translational lung cancer research"
        row["DOI"] = "10.21037/tlcr-23-157"
    elif aid == "A085":
        row["PMID"] = "37691861"
        row["Title"] = "Construction and validation of a novel ferroptosis-related prognostic signature for lung adenocarcinoma."
        row["First_author"] = "Yang"
        row["Year"] = "2023"
        row["Journal"] = "Translational lung cancer research"
        row["DOI"] = "10.21037/tlcr-23-351"
        set_formula(
            row,
            ["PRDX6", "ACSL3", "ARNTL", "ACSL4"],
            ["0.354", "0.237", "-0.547", "0.080"],
            "Risk score = 0.354 x PRDX6 expression + 0.237 x ACSL3 expression - 0.547 x ARNTL expression + 0.080 x ACSL4 expression",
        )
    elif aid == "A086":
        row["PMID"] = "38854940"
        row["Title"] = "Identification of a novel prognostic signature based on vitamin metabolism clustering-related genes in lung adenocarcinoma."
        row["First_author"] = "Chen"
        row["Year"] = "2024"
        row["Journal"] = "Translational lung cancer research"
        row["DOI"] = "10.21037/tlcr-24-245"
        row["GEO_accession"] = str(row.get("GEO_accession", "")).replace("GSE72904", "GSE72094")
        row["Validation_cohort"] = str(row.get("Validation_cohort", "")).replace("GSE72904", "GSE72094")
    elif aid == "A094":
        set_formula(
            row,
            ["LINC00115", "LINC00173", "LINC00968", "LINC01352"],
            ["-0.289", "-0.151", "-0.166", "-0.175"],
            "RiskScore = -0.289 x LINC00115 - 0.151 x LINC00173 - 0.166 x LINC00968 - 0.175 x LINC01352",
        )
    elif aid == "A098":
        row["Title"] = "Development and validation of a robust immune-related prognostic signature in lung adenocarcinoma"
        row["PMID"] = ""
    elif aid == "A097":
        row["PMID"] = "36971680"
        row["Title"] = "Construction and validation of a hypoxia-related risk signature identified EXO1 as a prognostic biomarker based on 12 genes in lung adenocarcinoma."
        row["First_author"] = "Chen"
        row["Year"] = "2023"
        row["Journal"] = "Aging"
        row["DOI"] = "10.18632/aging.204613"
        row["Source_files"] = "Construction and validation of a hyp Source Aging Albany NY SO 2023.pdf | Construction_and_validation_of_a_hyp_Source_Aging_Albany_NY_SO_2023.pdf"
        row["GEO_accession"] = "GSE19188; GSE10072"
        row["Training_cohort"] = "TCGA-LUAD"
        row["Validation_cohort"] = "GSE19188; GSE10072"
    elif aid == "A100":
        row["Gene_list"] = str(row.get("Gene_list", "")).replace("SERINE2", "SERPINE2")
        row["Coefficients"] = str(row.get("Coefficients", "")).replace("SERINE2", "SERPINE2")
        row["Risk_score_formula"] = str(row.get("Risk_score_formula", "")).replace("SERINE2", "SERPINE2")
    elif aid == "A110":
        set_formula(
            row,
            ["CDKN1A", "DDIT4", "IL33", "KRAS", "MTDH", "PANX1", "RRM2", "SLC2A1", "SLC7A5", "TLR4", "TXNRD1"],
            ["0.0049", "0.0037", "-0.0121", "0.0057", "0.0089", "0.0169", "0.0060", "0.0057", "0.0041", "-0.0380", "0.0003"],
            "Risk score = 0.0049 x CDKN1A expression + 0.0037 x DDIT4 expression - 0.0121 x IL33 expression + 0.0057 x KRAS expression + 0.0089 x MTDH expression + 0.0169 x PANX1 expression + 0.0060 x RRM2 expression + 0.0057 x SLC2A1 expression + 0.0041 x SLC7A5 expression - 0.0380 x TLR4 expression + 0.0003 x TXNRD1 expression",
        )
        row["Supplement_needed"] = "no"
        row["Cutoff_available"] = "yes"
        row["Cutoff_type"] = "median"
        row["Risk_group_definition"] = "High- and low-risk groups were defined by the median value of risk scores in the TCGA LUAD cohort; the median value was applied to the GEO validation cohort."
    elif aid == "A111":
        row["Gene_list"] = str(row.get("Gene_list", "")).replace("; AS", "; SERPINE1-AS1")
        row["Coefficients"] = str(row.get("Coefficients", "")).replace("; AS:0.104", "; SERPINE1-AS1:0.104")
        row["Risk_score_formula"] = str(row.get("Risk_score_formula", "")) + " [SERPINE1-AS1 coefficient from Supplementary Table context should be manually confirmed.]"
    elif aid == "A112":
        row["Cutoff_type"] = "median"
        row["Cutoff_available"] = "yes"
    elif aid == "A120":
        row["Gene_list"] = str(row.get("Gene_list", "")).replace("; FRONTIERS", "")
        row["Coefficients"] = re.sub(r";\s*FRONTIERS:[^;]+", "", str(row.get("Coefficients", "")))
        row["Gene_number"] = len([x for x in str(row["Gene_list"]).split(";") if x.strip()])
        row["Coefficients"] = str(row.get("Coefficients", "")).replace("AC027020.2:2022.951311", "AC027020.2:-0.918")
        row["Risk_score_formula"] = str(row.get("Risk_score_formula", "")).replace("−0.918 × Frontiers in Genetics frontiersin.org04 Fang et al. 10.3389/fgene.2022.951311 AC027020.2", "−0.918 × AC027020.2")
    elif aid == "A123":
        set_formula(
            row,
            ["HNRNPA2B1", "IGF2BP2", "NSUN4", "ALYREF"],
            ["0.383", "0.0548", "-0.595", "0.113"],
            "Risk score = 0.383 x HNRNPA2B1 expression + 0.0548 x IGF2BP2 expression - 0.595 x NSUN4 expression + 0.113 x ALYREF expression",
        )
    elif aid == "A131":
        row["Signature_name"] = "PTAAMG-Sig"
        set_formula(
            row,
            ["PPAT", "MIF", "GCLC", "PSPH", "KYNU", "ALDH2", "ACAD8", "HDC", "TYRP1"],
            ["0.0717", "0.0649", "0.0095", "0.1195", "0.1262", "-0.1394", "-0.0422", "-0.1782", "-0.1274"],
            "PTAAMG-Sig = 0.0717 x EV(PPAT) + 0.0649 x EV(MIF) + 0.0095 x EV(GCLC) + 0.1195 x EV(PSPH) + 0.1262 x EV(KYNU) - 0.1394 x EV(ALDH2) - 0.0422 x EV(ACAD8) - 0.1782 x EV(HDC) - 0.1274 x EV(TYRP1)",
        )
        row["Cutoff_available"] = "yes"
        row["Cutoff_value"] = "-0.366"
        row["Cutoff_type"] = "training-derived / optimal"
        row["Risk_group_definition"] = "High-risk group was defined by PTAAMG-Sig risk score above the -0.366 cutoff in TCGA-LUAD."
        row["Supplement_needed"] = "no"
        row["Reported_HR"] = "2.718 (95% CI 1.937-3.815)"
        row["Reported_C_index"] = "0.641"
        row["Reported_AUC"] = "1-year AUC 0.730"
        row["Exclusion_reason"] = ""
        row["Screening_level"] = "Level 1"
        row["Screening_label"] = "Included_reproducible"
        row["Reproducible"] = "yes"
    elif aid == "A148":
        row["Signature_name"] = "ERSS"
        set_formula(
            row,
            ["SERPINH1", "DSG2", "GPR37", "PCSK9", "TRPA1", "F2", "CDKN3", "DMD", "NR3C2", "CFTR", "CYP1A2", "MAPT", "SYT2", "CYP2D6", "SCN4A", "NUPR1", "PIK3CG", "DERL3"],
            ["0.13882291", "0.09216343", "0.07294183", "0.06978185", "0.06473978", "0.03490076", "0.02606376", "-0.01049731", "-0.01071941", "-0.01443196", "-0.01872414", "-0.01915136", "-0.03551274", "-0.03613574", "-0.0400695", "-0.0432351", "-0.07008271", "-0.0722234"],
            "ERSS = 0.13882291 x SERPINH1 + 0.09216343 x DSG2 + 0.07294183 x GPR37 + 0.06978185 x PCSK9 + 0.06473978 x TRPA1 + 0.03490076 x F2 + 0.02606376 x CDKN3 - 0.01049731 x DMD - 0.01071941 x NR3C2 - 0.01443196 x CFTR - 0.01872414 x CYP1A2 - 0.01915136 x MAPT - 0.03551274 x SYT2 - 0.03613574 x CYP2D6 - 0.0400695 x SCN4A - 0.0432351 x NUPR1 - 0.07008271 x PIK3CG - 0.0722234 x DERL3",
        )
        row["Cutoff_available"] = "yes"
        row["Cutoff_type"] = "median"
        row["Risk_group_definition"] = "High- and low-risk groups were defined by the median ERSS/risk score."
        row["Supplement_needed"] = "no"
        row["Exclusion_reason"] = ""
        row["Screening_level"] = "Level 1"
        row["Screening_label"] = "Included_reproducible"
        row["Reproducible"] = "yes"
    elif aid == "A160":
        row["Cancer_type"] = "LUAD"
        row["LUAD_specific"] = "yes"
        row["Omics_type"] = "methylation"
        row["Screening_level"] = "Level 5"
        row["Screening_label"] = "Exclude"
        row["Reproducible"] = "no"
        row["Exclusion_reason"] = "Primary model is a DNA methylation CpG signature, not a transcriptome/gene-expression score suitable for expression-matrix reproduction."
    elif aid == "A163":
        set_formula(
            row,
            ["MMADHC", "SOD2", "HTATIP2", "CYCS", "MPC1", "AK2", "MRPL44", "TBRG4", "TFAM", "TTC19", "CHCHD4", "SFXN1", "ABCD1", "NDUFAF7", "NOP2", "NME6", "XPNPEP3", "LIPT1", "METAP1D", "CA5B", "KMO", "ADHFE1"],
            ["0.04223596", "0.0910999", "0.08776809", "0.04116991", "-0.05948219", "-0.09171416", "0.05753366", "0.01824556", "0.09307276", "-0.03425415", "0.13106342", "0.09911277", "0.06724052", "-0.08541973", "-0.11770888", "-0.19223805", "0.18218013", "-0.09327797", "0.11149579", "-0.11764198", "-0.19248092", "-0.13196013"],
            "Risk score = 0.04223596 x MMADHC + 0.0910999 x SOD2 + 0.08776809 x HTATIP2 + 0.04116991 x CYCS - 0.05948219 x MPC1 - 0.09171416 x AK2 + 0.05753366 x MRPL44 + 0.01824556 x TBRG4 + 0.09307276 x TFAM - 0.03425415 x TTC19 + 0.13106342 x CHCHD4 + 0.09911277 x SFXN1 + 0.06724052 x ABCD1 - 0.08541973 x NDUFAF7 - 0.11770888 x NOP2 - 0.19223805 x NME6 + 0.18218013 x XPNPEP3 - 0.09327797 x LIPT1 + 0.11149579 x METAP1D - 0.11764198 x CA5B - 0.19248092 x KMO - 0.13196013 x ADHFE1",
        )
    elif aid == "A173":
        set_formula(
            row,
            ["ALOX15", "DDIT4", "HNF4A", "IL33", "GDF15"],
            ["-0.061165052", "0.169594473", "0.072356312", "-0.131023567", "-0.107882858"],
            "Risk score = -0.061165052 x ALOX15 expression + 0.169594473 x DDIT4 expression + 0.072356312 x HNF4A expression - 0.131023567 x IL33 expression - 0.107882858 x GDF15 expression",
        )

    # Better titles from filename when text title is clearly abstract/header debris.
    if re.match(r"^(background|abstract|lung cancer|prevalent form|/g1|signature to|lncRNA signature|m6A-related|Construction of an immunogenic cell$|A novel mitochondrial quality regula$)", str(row.get("Title", "")), re.I):
        row["Title"] = title_from_filename(row.get("Source_files", ""))

    # Try PubMed after title repair for missing/suspicious metadata.
    if not str(row.get("PMID", "")).strip() or str(row.get("PMID", "")).lower() == "nan":
        meta = lookup_pubmed(row.get("Title", ""), row.get("DOI", ""))
        if meta:
            for k in ["PMID", "Title", "First_author", "Year", "Journal", "DOI"]:
                if meta.get(k):
                    row[k] = meta[k]
    return row


def recompute_basic(row):
    if row.get("Screening_label") == "Exclude":
        return row
    if row.get("Screening_label") == "Pending_supplement":
        return row
    genes = [x.strip() for x in str(row.get("Gene_list", "")).split(";") if x.strip()]
    coefs = [x.strip() for x in str(row.get("Coefficients", "")).split(";") if x.strip()]
    coef_ok = row.get("Coefficients_available") == "yes" and len(coefs) >= 2
    formula_ok = row.get("Formula_available") == "yes" and bool(row.get("Risk_score_formula"))
    cutoff_ok = row.get("Cutoff_available") == "yes"
    cohort_ok = bool(row.get("Training_cohort")) and bool(row.get("Validation_cohort"))
    if len(genes) >= 2 and coef_ok and formula_ok and cutoff_ok and cohort_ok:
        row["Screening_level"] = "Level 1"
        row["Screening_label"] = "Included_reproducible"
        row["Reproducible"] = "yes"
        row["Exclusion_reason"] = ""
    return row


def main():
    rows = json.loads(ARTICLE_JSON.read_text(encoding="utf-8"))
    inv = pd.read_csv(INV_CSV)
    repaired = []
    for row in rows:
        text = source_text(row, inv)
        row = apply_manual_repairs(dict(row), text)
        row = recompute_basic(row)
        if row.get("Article_ID") in {
            "A033",
            "A060",
            "A062",
            "A069",
            "A070",
            "A084",
            "A085",
            "A086",
            "A094",
            "A097",
            "A098",
            "A100",
            "A110",
            "A111",
            "A112",
            "A120",
            "A123",
            "A131",
            "A148",
            "A160",
            "A163",
            "A173",
        }:
            note = "QC repaired/checked on 2026-05-30."
            row["Notes"] = norm((row.get("Notes", "") + " " + note).strip())
        repaired.append(row)

    if not any("Development and validation of a nomogram with an epigenetic signature" in r.get("Title", "") for r in repaired):
        repaired.append(
            {
                "Article_ID": "A178",
                "PMID": "33221751",
                "Title": "Development and validation of a nomogram with an epigenetic signature for predicting survival in patients with lung adenocarcinoma",
                "First_author": "Wang",
                "Year": "2020",
                "Journal": "Aging",
                "Cancer_type": "LUAD",
                "LUAD_specific": "yes",
                "Omics_type": "mRNA",
                "Endpoint": "OS",
                "Model_type": "LASSO-Cox / multivariate Cox / nomogram",
                "Training_cohort": "TCGA-LUAD",
                "Validation_cohort": "GSE31210",
                "Public_data_used": "TCGA / GEO",
                "GEO_accession": "GSE31210",
                "Signature_name": "11-gene epigenetic signature",
                "Gene_number": 11,
                "Gene_list": "DMAP1; ENY2; EPC1; GADD45A; HCFC2; PHC2; RCOR1; SMARCAL1; TLE2; TRIM28; ZNF516",
                "Coefficients_available": "yes",
                "Coefficients": "DMAP1:-0.070819821; ENY2:0.093606965; EPC1:-0.141271509; GADD45A:0.01034072; HCFC2:-0.356532015; PHC2:0.012487505; RCOR1:0.073312056; SMARCAL1:0.139640667; TLE2:-0.018668209; TRIM28:0.005275523; ZNF516:-0.088282786",
                "Formula_available": "yes",
                "Risk_score_formula": "risk score = -0.070819821 x DMAP1 level + 0.093606965 x ENY2 level - 0.141271509 x EPC1 level + 0.01034072 x GADD45A level - 0.356532015 x HCFC2 level + 0.012487505 x PHC2 level + 0.073312056 x RCOR1 level + 0.139640667 x SMARCAL1 level - 0.018668209 x TLE2 level + 0.005275523 x TRIM28 level - 0.088282786 x ZNF516 level",
                "Cutoff_available": "yes",
                "Cutoff_value": "",
                "Cutoff_type": "median",
                "Risk_group_definition": "Patients were clustered into high- and low-risk groups based on their median risk score.",
                "Reported_HR": "",
                "Reported_C_index": "",
                "Reported_AUC": "",
                "Independent_prognostic": "yes",
                "Nomogram_included": "yes",
                "External_validation": "yes",
                "Supplement_needed": "no",
                "Reproducible": "yes",
                "Exclusion_reason": "",
                "Notes": "QC added as separate article on 2026-05-30 after detecting erroneous merge with A097.",
                "Screening_level": "Level 1",
                "Screening_label": "Included_reproducible",
                "DOI": "10.18632/aging.104090",
                "Source_files": "Development and validation of a nomogram with an epigenetic .pdf",
                "Evidence_formula_context": "A risk score for each patient was then calculated as follows: risk score = (-0.070819821 x DMAP1 level) + ... + (-0.088282786 x ZNF516 level). Patients were clustered into high- and low-risk group based on their median risk score.",
                "Manual_check_needed": "no",
            }
        )

    OUT_JSON.write_text(json.dumps(repaired, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame(repaired).to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    print("wrote", OUT_JSON)
    print(pd.Series([r["Screening_label"] for r in repaired]).value_counts().to_string())


if __name__ == "__main__":
    main()
