import csv
import hashlib
import re
from pathlib import Path

from pypdf import PdfReader


ROOT = Path(__file__).resolve().parents[1]
PDF_DIR = ROOT / "source_pdfs"
WORK_DIR = ROOT / "work"
TEXT_DIR = WORK_DIR / "texts"
OUT_CSV = WORK_DIR / "pdf_inventory.csv"


def safe_stem(path: Path, idx: int) -> str:
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem).strip("_")
    return f"{idx:03d}_{stem[:140]}"


def clean_text(text: str) -> str:
    text = text.replace("\x00", " ")
    text = re.sub(r"-\n(?=[a-z])", "", text)
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def title_guess(text: str, metadata_title: str | None, filename: str) -> str:
    if metadata_title and len(metadata_title.strip()) > 10:
        return re.sub(r"\s+", " ", metadata_title.strip())
    lines = [re.sub(r"\s+", " ", ln).strip() for ln in text[:6000].splitlines()]
    lines = [ln for ln in lines if 20 <= len(ln) <= 220]
    skip = re.compile(
        r"^(abstract|introduction|original article|research article|article|www\.|doi:|open access|frontiers|scientific reports)",
        re.I,
    )
    for ln in lines[:60]:
        if skip.search(ln):
            continue
        if re.search(r"lung|adenocarcinoma|nsclc|progn|signature|model|risk", ln, re.I):
            return ln
    return Path(filename).stem


def main() -> None:
    TEXT_DIR.mkdir(parents=True, exist_ok=True)
    pdfs = sorted(PDF_DIR.rglob("*.pdf")) + sorted(PDF_DIR.rglob("*.PDF"))
    seen = set()
    rows = []
    for raw_idx, pdf in enumerate(pdfs, 1):
        if pdf in seen:
            continue
        seen.add(pdf)
        out_stem = safe_stem(pdf, len(rows) + 1)
        out_txt = TEXT_DIR / f"{out_stem}.txt"
        status = "ok"
        err = ""
        page_count = None
        metadata_title = ""
        text = ""
        try:
            reader = PdfReader(str(pdf))
            page_count = len(reader.pages)
            metadata_title = str(reader.metadata.title or "") if reader.metadata else ""
            chunks = []
            for page in reader.pages:
                try:
                    chunks.append(page.extract_text() or "")
                except Exception as exc:  # keep extracting other pages
                    chunks.append(f"\n[PAGE_EXTRACT_ERROR: {exc}]\n")
            text = clean_text("\n".join(chunks))
            out_txt.write_text(text, encoding="utf-8", errors="replace")
        except Exception as exc:
            status = "error"
            err = repr(exc)
            out_txt.write_text("", encoding="utf-8")

        text_hash = hashlib.sha256(text.encode("utf-8", errors="ignore")).hexdigest()
        first_2k_hash = hashlib.sha256(text[:2000].encode("utf-8", errors="ignore")).hexdigest()
        lower_name = pdf.name.lower()
        lower_text_head = text[:10000].lower()
        is_supp = bool(
            re.search(r"supplement|supplementary|supporting information|appendix", lower_name)
            or re.search(r"\bsupplementary (table|figure|material|methods)\b", lower_text_head)
        )
        rows.append(
            {
                "pdf_id": out_stem,
                "file_name": pdf.name,
                "file_path": str(pdf),
                "text_file": str(out_txt),
                "bytes": pdf.stat().st_size,
                "pages": page_count,
                "text_chars": len(text),
                "status": status,
                "error": err,
                "metadata_title": metadata_title,
                "title_guess": title_guess(text, metadata_title, pdf.name),
                "is_supplement_like": "yes" if is_supp else "no",
                "sha256_text": text_hash,
                "sha256_first_2k": first_2k_hash,
            }
        )

    WORK_DIR.mkdir(exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUT_CSV} with {len(rows)} rows")
    print(f"text dir: {TEXT_DIR}")


if __name__ == "__main__":
    main()
