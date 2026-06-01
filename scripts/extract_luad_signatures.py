import csv
import html
import json
import math
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
INV_CSV = WORK / "pdf_inventory.csv"
OUT_JSON = WORK / "article_extraction.json"
OUT_CSV = WORK / "article_extraction.csv"
PUBMED_CACHE = WORK / "pubmed_cache.json"


EXTRACTION_FIELDS = [
    "PMID",
    "Title",
    "First_author",
    "Year",
    "Journal",
    "Cancer_type",
    "LUAD_specific",
    "Omics_type",
    "Endpoint",
    "Model_type",
    "Training_cohort",
    "Validation_cohort",
    "Public_data_used",
    "GEO_accession",
    "Signature_name",
    "Gene_number",
    "Gene_list",
    "Coefficients_available",
    "Coefficients",
    "Formula_available",
    "Risk_score_formula",
    "Cutoff_available",
    "Cutoff_value",
    "Cutoff_type",
    "Risk_group_definition",
    "Reported_HR",
    "Reported_C_index",
    "Reported_AUC",
    "Independent_prognostic",
    "Nomogram_included",
    "External_validation",
    "Supplement_needed",
    "Reproducible",
    "Exclusion_reason",
    "Notes",
]


GENE_STOPWORDS = {
    "AUC",
    "ROC",
    "TCGA",
    "GEO",
    "GSE",
    "LUAD",
    "LUSC",
    "NSCLC",
    "OS",
    "DFS",
    "RFS",
    "PFS",
    "DSS",
    "HR",
    "CI",
    "KM",
    "LASSO",
    "Cox",
    "RNA",
    "DNA",
    "EMT",
    "TME",
    "ICGC",
    "FPKM",
    "TPM",
    "GO",
    "KEGG",
    "PCA",
    "GSVA",
    "GSEA",
    "NES",
    "FDR",
    "WHO",
    "USA",
    "AND",
    "FOR",
    "TO",
    "WITH",
    "BASED",
    "FIGURE",
    "TABLE",
    "FIG",
    "FIG.",
    "OF",
    "CUTOFF",
    "SUPPLE",
    "SUPPLEMENTARY",
    "EXPRESSION",
    "RESSION",
    "SION",
    "RISK",
    "SCORE",
    "OVERALL",
    "SENSITIVITY",
    "SPECIFICITY",
    "TIME",
    "TIME-DEPENDENT",
    "NEGATIVE",
    "POSITIVE",
    "BLOOD",
    "AGE",
    "STAGE",
    "MUTATION",
    "PERFORMANCEEVALUA",
    "PERFORMANCEEVALUATI",
    "PERFORMANCEEVALUATIONAND",
    "RESPECTIVELY",
    "GENDER",
    "HAZARD",
    "MALE",
    "FEMALE",
    "TABLE1",
    "AT",
    "AN",
    "MO",
    "RES",
}


def norm_space(s: str) -> str:
    return re.sub(r"\s+", " ", str(s or "")).strip()


def norm_key(s: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "", str(s or "").lower())
    return s[:180]


def title_similarity(a: str, b: str) -> float:
    stop = {
        "the",
        "and",
        "for",
        "with",
        "based",
        "source",
        "journal",
        "article",
        "lung",
        "adenocarcinoma",
    }
    ta = {x for x in re.findall(r"[a-z0-9]{4,}", str(a).lower()) if x not in stop}
    tb = {x for x in re.findall(r"[a-z0-9]{4,}", str(b).lower()) if x not in stop}
    if not ta or not tb:
        return 0.0
    return len(ta & tb) / max(1, min(len(ta), len(tb)))


def clean_doi(doi: str) -> str:
    doi = doi.strip().rstrip(".,);]")
    doi = doi.replace("%20", "")
    return doi


def find_doi(text: str) -> str:
    m = re.findall(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", text[:12000])
    if not m:
        m = re.findall(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", text)
    cleaned = [clean_doi(x) for x in m if len(clean_doi(x)) > 8]
    if not cleaned:
        return ""
    return Counter(cleaned).most_common(1)[0][0]


def contexts(text: str, patterns: list[str], window: int = 420, limit: int = 10) -> list[str]:
    out = []
    for pat in patterns:
        for m in re.finditer(pat, text, flags=re.I):
            start = max(0, m.start() - window)
            end = min(len(text), m.end() + window)
            snip = norm_space(text[start:end])
            if snip not in out:
                out.append(snip)
            if len(out) >= limit:
                return out
    return out


def first_match(text: str, pats: list[str], flags=re.I) -> str:
    for pat in pats:
        m = re.search(pat, text, flags)
        if m:
            return norm_space(m.group(1) if m.groups() else m.group(0))
    return ""


def find_pmids(text: str) -> str:
    hits = re.findall(r"\bPMID\s*:?\s*(\d{6,9})\b", text, flags=re.I)
    return hits[0] if hits else ""


def pubmed_fetch_by_id(pmid: str) -> dict:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        + urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "xml"})
    )
    with urllib.request.urlopen(url, timeout=20) as resp:
        xml = resp.read()
    root = ET.fromstring(xml)
    article = root.find(".//PubmedArticle")
    if article is None:
        return {}
    medline = article.find("MedlineCitation")
    art = medline.find("Article") if medline is not None else None
    title = "".join(art.findtext("ArticleTitle", default="").split()) if False else ""
    if art is not None and art.find("ArticleTitle") is not None:
        title = "".join(art.find("ArticleTitle").itertext())
    journal = ""
    year = ""
    if art is not None:
        journal = art.findtext("Journal/Title", default="") or art.findtext("Journal/ISOAbbreviation", default="")
        year = (
            art.findtext("Journal/JournalIssue/PubDate/Year", default="")
            or art.findtext("Journal/JournalIssue/PubDate/MedlineDate", default="")[:4]
        )
    first = ""
    auth = art.find("AuthorList/Author") if art is not None else None
    if auth is not None:
        first = auth.findtext("LastName", default="") or auth.findtext("CollectiveName", default="")
    doi = ""
    for aid in article.findall(".//ArticleId"):
        if aid.attrib.get("IdType") == "doi":
            doi = aid.text or ""
    return {
        "PMID": pmid,
        "Title": norm_space(title),
        "First_author": norm_space(first),
        "Year": norm_space(year[:4]),
        "Journal": norm_space(journal),
        "DOI": clean_doi(doi),
    }


def pubmed_search(term: str) -> str:
    params = {
        "db": "pubmed",
        "retmode": "json",
        "retmax": "1",
        "term": term,
        "sort": "relevance",
    }
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=20) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    ids = data.get("esearchresult", {}).get("idlist", [])
    return ids[0] if ids else ""


def load_pubmed_cache() -> dict:
    if PUBMED_CACHE.exists():
        return json.loads(PUBMED_CACHE.read_text(encoding="utf-8"))
    return {}


def save_pubmed_cache(cache: dict) -> None:
    PUBMED_CACHE.write_text(json.dumps(cache, ensure_ascii=False, indent=2), encoding="utf-8")


def pubmed_lookup(doi: str, title_guess: str, text: str, cache: dict) -> dict:
    direct_pmid = find_pmids(text)
    key = direct_pmid or doi or norm_key(title_guess)
    if key in cache:
        cached = cache[key]
        if not cached.get("Title") or title_similarity(cached.get("Title", ""), title_guess) >= 0.18:
            return cached
    meta = {}
    try:
        pmid = direct_pmid
        if not pmid and doi:
            pmid = pubmed_search(f"{doi}[AID]")
            time.sleep(0.12)
        if not pmid and title_guess:
            query_title = re.sub(r"\s+", " ", title_guess)
            query_title = re.sub(r"Source .*", "", query_title, flags=re.I)
            pmid = pubmed_search(f'"{query_title[:180]}"[Title]')
            time.sleep(0.12)
        if pmid:
            meta = pubmed_fetch_by_id(pmid)
            time.sleep(0.12)
        if meta.get("Title") and title_similarity(meta["Title"], title_guess) < 0.18:
            # The first DOI in a PDF can belong to a cited reference. Fall back to title search.
            meta = {}
            query_title = re.sub(r"\s+", " ", title_guess)
            query_title = re.sub(r"Source .*", "", query_title, flags=re.I)
            pmid = pubmed_search(f'"{query_title[:180]}"[Title]')
            time.sleep(0.12)
            if pmid:
                meta = pubmed_fetch_by_id(pmid)
                time.sleep(0.12)
            if meta.get("Title") and title_similarity(meta["Title"], title_guess) < 0.18:
                meta = {}
    except Exception as exc:
        meta = {"lookup_error": repr(exc)}
    cache[key] = meta
    save_pubmed_cache(cache)
    return meta


def infer_title_from_text(text: str, fallback: str) -> str:
    lines = [norm_space(x) for x in text[:7000].splitlines()]
    lines = [x for x in lines if 20 <= len(x) <= 240]
    bad = re.compile(r"^(abstract|introduction|keywords|original article|research article|frontiers|www\.|doi:)", re.I)
    for line in lines[:80]:
        if bad.search(line):
            continue
        if re.search(r"lung|adenocarcinoma|nsclc|progn|signature|model|risk", line, re.I):
            return line
    return fallback


def title_from_filename(filename: str) -> str:
    stem = Path(filename).stem.replace("_", " ")
    stem = re.sub(r"\s+Source\s+.*", "", stem, flags=re.I)
    stem = re.sub(r"\s+SO\s+\d{4}.*", "", stem, flags=re.I)
    return norm_space(stem)


def infer_cancer_type(text: str) -> tuple[str, str]:
    head = text[:50000]
    luad = bool(re.search(r"\bLUAD\b|lung adenocarcinoma|lung adeno-carcinoma", head, re.I))
    nsclc = bool(re.search(r"\bNSCLC\b|non-small cell lung cancer", head, re.I))
    pan = bool(re.search(r"pan-cancer|multiple cancers|across cancers", head, re.I))
    if luad and not nsclc:
        return "LUAD", "yes"
    if luad and nsclc:
        return "NSCLC/LUAD", "yes"
    if nsclc:
        return "NSCLC", "unclear"
    if pan:
        return "mixed", "unclear"
    return "unclear", "unclear"


def infer_omics_type(text: str, title: str) -> str:
    sample = f"{title} {text[:70000]}".lower()
    types = []
    if re.search(r"\blncrna|long non.?coding rna", sample):
        types.append("lncRNA")
    if re.search(r"\bmirna|micro.?rna|miR-", sample, re.I):
        types.append("miRNA")
    if re.search(r"\bcircrna|circular rna", sample):
        types.append("circRNA")
    if re.search(r"immune|immun", sample):
        types.append("immune")
    if re.search(r"metabolism|metabolic|glycolysis|fatty acid|lactate|propionate|glutathione", sample):
        types.append("metabolism")
    if re.search(r"methylation|epigenetic", sample) and not re.search(r"gene expression|transcriptome|rna-seq|mrna", sample):
        types.append("methylation")
    if re.search(r"multi.?omics|multi omics|proteom", sample):
        types.append("multi-omics")
    if not any(t in types for t in ["lncRNA", "miRNA", "circRNA"]) and re.search(
        r"\bmrna\b|gene expression|transcriptom|rna-seq|rna seq|expression profile", sample
    ):
        types.insert(0, "mRNA")
    return " / ".join(dict.fromkeys(types)) if types else "unclear"


def infer_endpoint(text: str) -> str:
    candidates = []
    for label, pat in [
        ("OS", r"\boverall survival\b|\bOS\b"),
        ("DFS", r"\bdisease[- ]free survival\b|\bDFS\b"),
        ("RFS", r"\brecurrence[- ]free survival\b|\brelapse[- ]free survival\b|\bRFS\b"),
        ("PFS", r"\bprogression[- ]free survival\b|\bPFS\b"),
        ("DSS", r"\bdisease[- ]specific survival\b|\bDSS\b"),
    ]:
        if re.search(pat, text[:120000], re.I):
            candidates.append(label)
    return " / ".join(candidates) if candidates else "unclear"


def infer_model_type(text: str) -> str:
    s = text[:140000]
    types = []
    if re.search(r"LASSO", s, re.I) and re.search(r"Cox", s, re.I):
        types.append("LASSO-Cox")
    if re.search(r"multivariate Cox|multivariable Cox", s, re.I):
        types.append("multivariate Cox")
    elif re.search(r"\bCox\b", s, re.I):
        types.append("Cox")
    if re.search(r"random survival forest|\bRSF\b", s, re.I):
        types.append("RSF")
    if re.search(r"machine learning|support vector|xgboost|random forest|elastic net|neural network", s, re.I):
        types.append("ML")
    if re.search(r"nomogram", s, re.I):
        types.append("nomogram")
    return " / ".join(dict.fromkeys(types)) if types else "other"


def infer_cohorts(text: str) -> tuple[str, str, str, str]:
    geo = sorted(set(re.findall(r"\bGSE\d{3,8}\b", text, flags=re.I)), key=lambda x: (len(x), x))
    geo = [x.upper() for x in geo]
    public = []
    if re.search(r"\bTCGA\b|The Cancer Genome Atlas", text, re.I):
        public.append("TCGA")
    if geo or re.search(r"\bGEO\b|Gene Expression Omnibus", text, re.I):
        public.append("GEO")
    if re.search(r"\bICGC\b|International Cancer Genome Consortium", text, re.I):
        public.append("ICGC")
    if re.search(r"\bCPTAC\b", text, re.I):
        public.append("CPTAC")

    train_ctx = contexts(text, [r"training cohort", r"training set", r"train(?:ed|ing) data", r"TCGA[- ]LUAD"], 260, 8)
    valid_ctx = contexts(text, [r"validation cohort", r"validation set", r"external validation", r"validated in", r"test set"], 260, 8)

    def summarize(ctxs):
        found = []
        for c in ctxs:
            if re.search(r"\bTCGA\b|TCGA[- ]LUAD", c, re.I):
                found.append("TCGA-LUAD")
            found.extend(re.findall(r"\bGSE\d{3,8}\b", c, flags=re.I))
            if re.search(r"\bICGC\b", c, re.I):
                found.append("ICGC")
            if re.search(r"internal|clinical|hospital|institution", c, re.I):
                found.append("internal cohort")
        return "; ".join(dict.fromkeys(x.upper() if x.lower().startswith("gse") else x for x in found))

    training = summarize(train_ctx)
    validation = summarize(valid_ctx)

    if not training:
        if re.search(r"TCGA[- ]LUAD|TCGA LUAD|LUAD.*TCGA", text[:100000], re.I):
            training = "TCGA-LUAD"
    if not validation:
        # If multiple public cohorts occur and validation language is present, list GEO cohorts as likely validation.
        if geo and re.search(r"validat|test set|external", text[:160000], re.I):
            validation = "; ".join(dict.fromkeys(geo[:8]))
    return training, validation, " / ".join(public), "; ".join(dict.fromkeys(geo))


def formula_snippets(text: str) -> list[str]:
    pats = [
        r"risk score",
        r"riskscore",
        r"risk-score",
        r"prognostic score",
        r"score\s*=",
        r"formula",
        r"coefficient",
    ]
    snippets = contexts(text, pats, 700, 14)
    ranked = sorted(
        snippets,
        key=lambda s: (
            bool(re.search(r"risk score\s*=|score\s*=|coef|coefficient|\*|\+|\(|\)", s, re.I)),
            len(re.findall(r"[-+]?\d+\.\d+", s)),
        ),
        reverse=True,
    )
    return ranked[:8]


def likely_formula(snippets: list[str]) -> str:
    for snip in snippets:
        if not re.search(r"risk\s*score|prognostic\s*score|riskscore", snip, re.I):
            continue
        if not (re.search(r"[-+]?\d+\.\d+", snip) and re.search(r"\+|\*|×|¡Á|coefficient|coef", snip, re.I)):
            continue
        m = re.search(r"(risk\s*score|riskscore|prognostic\s*score)\s*(?:=|/equals|:)", snip, flags=re.I)
        if not m:
            continue
        start = m.start() if m else max(0, re.search(r"risk\s*score|riskscore|prognostic\s*score", snip, re.I).start())
        seg = snip[start : start + 2200]
        if len(re.findall(r"\+|×|\*|¡Á", seg)) < 2:
            continue
        # Stop after the equation but do not treat decimal points as sentence endings.
        stop = re.search(
            r"\s(?:Patients were|The TCGA|Kaplan|Risk scores were|Based on the expression|According to the median|We then|Subsequently|Thereafter|Figure\s+\d|Fig\.|Performance\s*evaluation|Performanceevaluation)\b|\d+\.\d+\s*Performance",
            seg,
            flags=re.I,
        )
        if stop and stop.start() > 80:
            seg = seg[: stop.start()]
        return norm_space(seg)
    return ""


def parse_gene_coeffs(formula: str) -> tuple[list[str], list[str]]:
    if not formula:
        return [], []
    f = formula.replace("×", "*").replace("¡Á", "*").replace("−", "-").replace("–", "-")
    # pypdf sometimes inserts spaces inside symbols (FER MT1, AC AT2).
    f = re.sub(r"\b([A-Z]{1,5})\s+([A-Z]{0,3}\d[A-Z0-9-]{0,8})\b", r"\1\2", f)
    f = re.sub(r"\b([A-Z]{1,5}\d)\s+([A-Z0-9-]{1,8})\b", r"\1\2", f)
    # Coeff near gene, both common orders.
    gene_re = r"[A-Z][A-Z0-9.]*(?:-[A-Z][A-Z0-9.]*)?"
    expr_re = r"(?:Exp(?:ression)?|expres\s*sion|expression)(?:\s+(?:level|value))?(?:\s+of)?"
    patterns = [
        rf"\(?\s*([-+]?\d+\.\d+(?:e[-+]?\d+)?)\s*\)?\s*[\*\u00d7]?\s*(?:[xX]\s*)?(?:{expr_re}\s*)?[\(\[]?({gene_re})[\)\]]?",
        rf"({gene_re})\s*(?:exp|expression)?\.?\s*[\)\]]?\s*[\*\u00d7]\s*\(?\s*([-+]?\d+\.\d+(?:e[-+]?\d+)?)\s*\)?",
        rf"({gene_re})\s*\(([-+]?\d+\.\d+(?:e[-+]?\d+)?)\)",
    ]
    genes = []
    coeffs = []
    for pat in patterns:
        for m in re.finditer(pat, f, re.I):
            a, b = m.group(1), m.group(2)
            if re.match(r"[-+]?\d", a):
                coef, gene = a, b
            else:
                gene, coef = a, b
            gene = gene.strip().upper()
            gene = re.sub(r"[-+]\d.*$", "", gene)
            gene = gene.strip(" .;,:+-_()[]")
            gene = re.sub(r"(EXP|EXPRESSION)$", "", gene)
            gene = gene.strip(" .;,:+-_()[]")
            if gene in GENE_STOPWORDS or len(gene) < 2:
                continue
            if re.fullmatch(r"E-?\d+", gene):
                continue
            if not re.search(r"\d|\.", coef):
                continue
            item = f"{gene}:{coef}"
            if gene not in genes:
                genes.append(gene)
                coeffs.append(item)
    return genes, coeffs


def gene_list_from_text(text: str, formula_genes: list[str]) -> list[str]:
    if formula_genes:
        return formula_genes
    pats = [
        r"(?:signature|model|risk model|prognostic model).*?(?:consisted of|composed of|included|contains|comprised)\s+(?:\d+\s+)?(?:genes|lncRNAs|miRNAs|circRNAs)?[:：]?\s*([A-Za-z0-9,; /\-()]+)",
        r"(?:identified|screened)\s+(?:\d+\s+)?(?:genes|lncRNAs|miRNAs|circRNAs).*?[:：]\s*([A-Za-z0-9,; /\-()]+)",
    ]
    for pat in pats:
        m = re.search(pat, text[:160000], re.I | re.S)
        if not m:
            continue
        chunk = norm_space(m.group(1)[:900])
        hits = re.findall(r"\b[A-Z][A-Z0-9-]{1,14}\b", chunk)
        hits = [h.upper() for h in hits if h.upper() not in GENE_STOPWORDS]
        if len(hits) >= 2:
            return list(dict.fromkeys(hits))[:80]
    return []


def infer_signature_name(title: str, text: str, omics: str) -> str:
    m = re.search(r"((?:\d+|[A-Za-z]+)[- ](?:gene|lncRNA|miRNA|circRNA)[- ](?:prognostic )?signature)", title, re.I)
    if m:
        return norm_space(m.group(1))
    m = re.search(r"(\d+[- ](?:gene|lncRNA|miRNA|circRNA)[- ](?:signature|model))", text[:50000], re.I)
    if m:
        return norm_space(m.group(1))
    if re.search(r"signature|risk model|prognostic model", title, re.I):
        return title[:180]
    return ""


def infer_cutoff(text: str) -> tuple[str, str, str, str]:
    cut_ctxs = contexts(text, [r"cut-?off", r"median risk", r"median value", r"X-tile", r"Youden", r"surv_cutpoint", r"high-risk", r"low-risk"], 360, 12)
    combined = " ".join(cut_ctxs)
    ctype = ""
    if re.search(r"median", combined, re.I):
        ctype = "median"
    if re.search(r"optimal|optimum|best cutoff|best cut-off", combined, re.I):
        ctype = "optimal" if not ctype else ctype + " / optimal"
    if re.search(r"Youden", combined, re.I):
        ctype = "Youden" if not ctype else ctype + " / Youden"
    if re.search(r"surv_cutpoint|survminer", combined, re.I):
        ctype = "surv_cutpoint" if not ctype else ctype + " / surv_cutpoint"
    if re.search(r"X-tile", combined, re.I):
        ctype = "X-tile" if not ctype else ctype + " / X-tile"
    if re.search(r"training set|training cohort|derived from training", combined, re.I):
        ctype = "training-derived" if not ctype else ctype + " / training-derived"
    value = first_match(
        combined,
        [
            r"cut-?off(?: value)?(?: was| is| of| =|:)?\s*([-+]?\d+(?:\.\d+)?)",
            r"threshold(?: value)?(?: was| is| of| =|:)?\s*([-+]?\d+(?:\.\d+)?)",
            r"median(?: risk score| value)?(?: was| is| of| =|:)?\s*([-+]?\d+(?:\.\d+)?)",
        ],
    )
    definition = ""
    for c in cut_ctxs:
        if re.search(r"high-risk|low-risk", c, re.I):
            definition = c[:900]
            break
    available = "yes" if ctype or value or definition else "no"
    return available, value, ctype or "unclear", definition


def infer_metrics(text: str) -> tuple[str, str, str]:
    metric_ctx = " ".join(contexts(text, [r"hazard ratio", r"\bHR\b", r"C-index", r"\bAUC\b"], 180, 16))
    hr = "; ".join(dict.fromkeys(re.findall(r"(?:HR|hazard ratio)[^\d]{0,20}(\d+\.\d+(?:\s*\([^)]*\))?)", metric_ctx, re.I)))[:500]
    cindex = "; ".join(dict.fromkeys(re.findall(r"C[- ]?index[^\d]{0,20}(\d+\.\d+)", metric_ctx, re.I)))[:500]
    auc = "; ".join(
        dict.fromkeys(
            re.findall(
                r"(?:(?:1|3|5)[- ]year[^.;]{0,40})?AUC[^0-9]{0,12}(\d+\.\d+)",
                metric_ctx,
                re.I,
            )
        )
    )[:700]
    return hr, cindex, auc


def exclusion_and_level(row: dict, text: str) -> tuple[str, str, str, str]:
    title = row["Title"]
    omics = row["Omics_type"]
    luad_ok = row["LUAD_specific"] == "yes" or row["Cancer_type"].startswith("LUAD")
    expression_ok = bool(re.search(r"mRNA|lncRNA|miRNA|circRNA|immune|metabolism|multi-omics", omics))
    prognosis_ok = row["Endpoint"] != "unclear" or re.search(r"prognos|survival", title + " " + text[:40000], re.I)
    try:
        gene_n = int(row["Gene_number"]) if row["Gene_number"] not in ("", None) and not (isinstance(row["Gene_number"], float) and math.isnan(row["Gene_number"])) else 0
    except Exception:
        gene_n = len([x for x in str(row.get("Gene_list", "")).split(";") if x.strip()])
    gene_ok = bool(row["Gene_list"]) and gene_n >= 2
    coef_ok = row["Coefficients_available"] == "yes" or row["Formula_available"] == "yes"
    cutoff_ok = row["Cutoff_available"] == "yes"
    cohort_ok = bool(row["Training_cohort"]) and bool(row["Validation_cohort"])
    supp_hint = row["Supplement_needed"] == "yes"
    cluster_only = bool(re.search(r"subtype|cluster", title, re.I)) and not re.search(r"risk score|signature|risk model", text[:160000], re.I)
    methyl_only = "methylation" in omics and not expression_ok
    diagnosis_only = bool(re.search(r"diagnos", title, re.I)) and not prognosis_ok
    review = bool(re.search(r"\breview\b|meta-analysis|systematic review", title, re.I))

    if review:
        return "Level 5", "Exclude", "no", "Review/meta-analysis, not an original reproducible LUAD expression prognostic model."
    if not luad_ok:
        return "Level 5", "Exclude", "no", "LUAD is not clearly the target cancer type or LUAD-specific results cannot be separated."
    if methyl_only:
        return "Level 5", "Exclude", "no", "Primary model is methylation/epigenetic rather than transcriptome/gene-expression score."
    if not prognosis_ok:
        return "Level 5", "Exclude", "no", "No clear survival prognosis endpoint."
    if not expression_ok:
        return "Level 5", "Exclude", "no", "Not clearly based on transcriptomic/gene-expression data."
    if cluster_only:
        return "Level 4", "Relevant_not_reproducible", "no", "Subtype/cluster analysis without a clear risk-score formula."
    if gene_ok and coef_ok and cutoff_ok and cohort_ok:
        return "Level 1", "Included_reproducible", "yes", ""
    if gene_ok and supp_hint and cohort_ok:
        return "Level 2", "Pending_supplement", "partial", "Formula/coefficient or cutoff appears to require supplementary material."
    if gene_ok and (not coef_ok or not cutoff_ok):
        missing = []
        if not coef_ok:
            missing.append("complete coefficient/formula")
        if not cutoff_ok:
            missing.append("cutoff/grouping rule")
        if not cohort_ok:
            missing.append("clear training/validation cohort")
        return "Level 3", "Partially_reproducible", "partial", "Missing " + ", ".join(missing) + "; alternative gene-set style reproduction may be possible."
    return "Level 4", "Relevant_not_reproducible", "no", "Relevant LUAD expression prognosis article, but gene list and/or risk-score parameters are not sufficiently extractable from text."


def build_articles() -> list[dict]:
    inv = pd.read_csv(INV_CSV)
    texts = []
    for _, r in inv.iterrows():
        text = Path(r["text_file"]).read_text(encoding="utf-8", errors="ignore")
        doi = find_doi(text)
        title_guess = norm_space(r.get("title_guess", "")) or infer_title_from_text(text, r["file_name"])
        key = f"doi:{doi.lower()}" if doi else f"title:{norm_key(title_guess)}"
        texts.append({**r.to_dict(), "text": text, "doi": doi, "article_key": key, "title_guess2": title_guess})

    by_key: dict[str, list[dict]] = {}
    for item in texts:
        # Exact duplicate text with different file name will share DOI/title; keep grouped.
        by_key.setdefault(item["article_key"], []).append(item)

    cache = load_pubmed_cache()
    rows = []
    for idx, (key, group) in enumerate(sorted(by_key.items()), 1):
        group_sorted = sorted(group, key=lambda x: int(x.get("text_chars", 0)), reverse=True)
        combined_text = "\n\n".join(g["text"] for g in group_sorted)
        primary = group_sorted[0]
        doi = primary["doi"] or find_doi(combined_text)
        title_guess = primary["title_guess2"]
        meta = pubmed_lookup(doi, title_guess, combined_text, cache)
        title = meta.get("Title") or infer_title_from_text(combined_text, title_guess)
        if not meta.get("Title") and re.search(r"^(background|lung cancer remains|abstract\b)", title, re.I):
            title = title_from_filename(primary["file_name"]) or title
        first_author = meta.get("First_author", "")
        year = meta.get("Year", "") or first_match(primary["file_name"], [r"\b(20\d{2}|19\d{2})\b"])
        journal = meta.get("Journal", "") or first_match(primary["file_name"], [r"Source\s+(.+?)\s+SO\s+\d{4}", r"Source\s+(.+?)\s+\d{4}"])

        cancer, luad_specific = infer_cancer_type(combined_text)
        omics = infer_omics_type(combined_text, title)
        endpoint = infer_endpoint(combined_text)
        model = infer_model_type(combined_text)
        training, validation, public_data, geo = infer_cohorts(combined_text)
        snippets = formula_snippets(combined_text)
        formula = likely_formula(snippets)
        formula_genes, coeffs = parse_gene_coeffs(formula)
        if len(coeffs) < 2:
            formula_for_repro = ""
            coeffs_for_repro = []
            formula_genes_for_repro = []
        else:
            formula_for_repro = formula
            coeffs_for_repro = coeffs
            formula_genes_for_repro = formula_genes
        genes = gene_list_from_text(combined_text, formula_genes_for_repro)
        cutoff_avail, cutoff_value, cutoff_type, risk_def = infer_cutoff(combined_text)
        hr, cindex, auc = infer_metrics(combined_text)

        supplement_needed = "yes" if re.search(r"supplementary table|supplemental table|Table S\d|Additional file|Supplementary Data", combined_text, re.I) and not coeffs else "no"
        independent = "yes" if re.search(r"independent prognostic|independent predictor|multivariate Cox.*independent", combined_text[:180000], re.I | re.S) else "no"
        nomogram = "yes" if re.search(r"nomogram", combined_text[:180000], re.I) else "no"
        external = "yes" if validation and re.search(r"external|GSE|GEO|ICGC|independent cohort|validation cohort", validation + " " + combined_text[:140000], re.I) else "no"

        row = {
            "Article_ID": f"A{idx:03d}",
            "PMID": meta.get("PMID", ""),
            "Title": title,
            "First_author": first_author,
            "Year": year,
            "Journal": journal,
            "Cancer_type": cancer,
            "LUAD_specific": luad_specific,
            "Omics_type": omics,
            "Endpoint": endpoint,
            "Model_type": model,
            "Training_cohort": training,
            "Validation_cohort": validation,
            "Public_data_used": public_data,
            "GEO_accession": geo,
            "Signature_name": infer_signature_name(title, combined_text, omics),
            "Gene_number": len(genes) if genes else "",
            "Gene_list": "; ".join(genes),
            "Coefficients_available": "yes" if coeffs_for_repro else "no",
            "Coefficients": "; ".join(coeffs_for_repro),
            "Formula_available": "yes" if formula_for_repro else "no",
            "Risk_score_formula": formula_for_repro or formula,
            "Cutoff_available": cutoff_avail,
            "Cutoff_value": cutoff_value,
            "Cutoff_type": cutoff_type if cutoff_avail == "yes" else "",
            "Risk_group_definition": risk_def,
            "Reported_HR": hr,
            "Reported_C_index": cindex,
            "Reported_AUC": auc,
            "Independent_prognostic": independent,
            "Nomogram_included": nomogram,
            "External_validation": external,
            "Supplement_needed": supplement_needed,
            "Reproducible": "",
            "Exclusion_reason": "",
            "Notes": "",
            "Screening_level": "",
            "Screening_label": "",
            "DOI": doi or meta.get("DOI", ""),
            "Source_files": " | ".join(g["file_name"] for g in group_sorted),
            "Evidence_formula_context": " || ".join(snippets[:4])[:4000],
            "Manual_check_needed": "yes" if (not formula and re.search(r"formula|coefficient|risk score", combined_text, re.I)) else "no",
        }
        level, label, reproducible, reason = exclusion_and_level(row, combined_text)
        row["Screening_level"] = level
        row["Screening_label"] = label
        row["Reproducible"] = reproducible
        row["Exclusion_reason"] = reason
        notes = []
        if meta.get("lookup_error"):
            notes.append("PubMed lookup failed: " + meta["lookup_error"])
        if len(group_sorted) > 1:
            notes.append(f"Merged {len(group_sorted)} PDFs by DOI/title.")
        if row["Manual_check_needed"] == "yes" and label in {"Pending_supplement", "Partially_reproducible", "Relevant_not_reproducible"}:
            notes.append("Formula/coefficient terms detected but complete machine-readable equation not found; manual PDF/image/supplement check recommended.")
        row["Notes"] = " ".join(notes)
        rows.append(row)
    return rows


def main() -> None:
    rows = build_articles()
    OUT_JSON.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    fieldnames = ["Article_ID", *EXTRACTION_FIELDS, "Screening_level", "Screening_label", "DOI", "Source_files", "Manual_check_needed", "Evidence_formula_context"]
    with OUT_CSV.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUT_CSV} with {len(rows)} article-level rows")
    print(Counter(r["Screening_label"] for r in rows))


if __name__ == "__main__":
    main()
