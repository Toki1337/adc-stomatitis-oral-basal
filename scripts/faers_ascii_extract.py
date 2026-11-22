"""Stream-download and parse official FDA FAERS Quarterly ASCII extracts (2012Q4-2026Q3)
into compact per-quarter ADC case extracts and whole-FAERS deduplication arrays under
dl/faers/ascii_extract/ (Track 2 raw regeneration; Track 1 uses results/faers_counts.tsv).
"""

import concurrent.futures
import gzip
import io
import json
import re
import sys
import time
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "dl" / "faers" / "ascii_extract"
OUT_DIR.mkdir(parents=True, exist_ok=True)
INDEX_URL = "https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html"
UA = "Mozilla/5.0 (compatible; AcademicPharmacovigilanceResearch/1.0)"

OUT_GROUPS = [
    ("oral", {
        "STOMATITIS", "MUCOSAL INFLAMMATION", "MOUTH ULCERATION",
        "APHTHOUS ULCER", "GLOSSITIS", "ORAL MUCOSAL ERYTHEMA",
        "ORAL MUCOSAL BLISTERING",
    }),
    ("neutropenia", {
        "NEUTROPENIA", "FEBRILE NEUTROPENIA", "NEUTROPHIL COUNT DECREASED",
    }),
    ("neuropathy", {
        "PERIPHERAL SENSORY NEUROPATHY", "NEUROPATHY PERIPHERAL",
        "PERIPHERAL MOTOR NEUROPATHY", "POLYNEUROPATHY",
    }),
    ("skin", {
        "RASH", "RASH MACULO-PAPULAR", "RASH PRURITIC", "RASH ERYTHEMATOUS",
        "DERMATITIS", "SKIN EXFOLIATION", "STEVENS-JOHNSON SYNDROME",
        "TOXIC EPIDERMAL NECROLYSIS",
    }),
    ("ocular", {
        "KERATOPATHY", "KERATITIS", "VISION BLURRED", "DRY EYE",
        "CONJUNCTIVITIS", "VISUAL IMPAIRMENT", "CORNEAL DISORDER",
    }),
    ("ild", {
        "INTERSTITIAL LUNG DISEASE", "PNEUMONITIS",
    }),
    ("diarrhea", {
        "DIARRHOEA", "DIARRHEA",
    }),
]
PT_TO_BIT = {}
for bit_idx, (_, pts) in enumerate(OUT_GROUPS):
    for pt in pts:
        PT_TO_BIT[pt] = PT_TO_BIT.get(pt, 0) | (1 << bit_idx)

ADC_REGEX = {
    "trastuzumab emtansine": re.compile(
        r"\b(?:KADCYLA|ADO[-\s]?TRASTUZUMAB\s+EMTANSINE|TRASTUZUMAB[-\s]+EMTANSINE|T-?DM1)\b"
    ),
    "trastuzumab deruxtecan": re.compile(
        r"\b(?:ENHERTU|FAM[-\s]?TRASTUZUMAB\s+DERUXTECAN(?:-[A-Z0-9]+)?|TRASTUZUMAB[-\s]+DERUXTECAN(?:-[A-Z0-9]+)?|T-?DXD)\b"
    ),
    "datopotamab deruxtecan": re.compile(
        r"\b(?:DATROWAY|DATOPOTAMAB(?:\s+DERUXTECAN(?:-[A-Z0-9]+)?)?|DS-?1062A?|DATO-?DXD)\b"
    ),
    "sacituzumab govitecan": re.compile(
        r"\b(?:TRODELVY|SACITUZUMAB\s+GOVITECAN(?:-[A-Z0-9]+)?|IMMU-?132)\b"
    ),
    "enfortumab vedotin": re.compile(
        r"\b(?:PADCEV|ENFORTUMAB(?:\s+VEDOTIN(?:-[A-Z0-9]+)?)?|ASG-?22ME)\b"
    ),
    "tisotumab vedotin": re.compile(
        r"\b(?:TIVDAK|TISOTUMAB(?:\s+VEDOTIN(?:-[A-Z0-9]+)?)?)\b"
    ),
    "mirvetuximab soravtansine": re.compile(
        r"\b(?:ELAHERE|MIRVETUXIMAB(?:\s+SORAVTANSINE(?:-[A-Z0-9]+)?)?|IMGN-?853)\b"
    ),
    "brentuximab vedotin": re.compile(
        r"\b(?:ADCETRIS|BRENTUXIMAB(?:\s+VEDOTIN(?:-[A-Z0-9]+)?)?|SGN-?35)\b"
    ),
    "polatuzumab vedotin": re.compile(
        r"\b(?:POLIVY|POLATUZUMAB(?:\s+VEDOTIN(?:-[A-Z0-9]+)?)?)\b"
    ),
    "loncastuximab tesirine": re.compile(
        r"\b(?:ZYNLONTA|LONCASTUXIMAB(?:\s+TESIRINE(?:-[A-Z0-9]+)?)?|ADCT-?402)\b"
    ),
    "inotuzumab ozogamicin": re.compile(
        r"\b(?:BESPONSA|INOTUZUMAB(?:\s+OZOGAMICIN)?)\b"
    ),
    "gemtuzumab ozogamicin": re.compile(
        r"\b(?:MYLOTARG|GEMTUZUMAB(?:\s+OZOGAMICIN)?)\b"
    ),
    "belantamab mafodotin": re.compile(
        r"\b(?:BLENREP|BELANTAMAB(?:\s+MAFODOTIN(?:-[A-Z0-9]+)?)?)\b"
    ),
    "telisotuzumab vedotin": re.compile(
        r"\b(?:EMRELIS|TELISOTUZUMAB\s+VEDOTIN(?:-[A-Z0-9]+)?|ABBV-?399)\b"
    ),
    "disitamab vedotin": re.compile(
        r"\b(?:AIDIXI|DISITAMAB(?:\s+VEDOTIN)?|RC-?48)\b"
    ),
}
FAST_ADC_RE = re.compile(
    r"(?:KADCYLA|EMTANSINE|T-?DM1|ENHERTU|DERUXTECAN|T-?DXD|DATROWAY|DATOPOTAMAB|DS-?1062|DATO-?DXD|"
    r"TRODELVY|SACITUZUMAB|IMMU-?132|PADCEV|ENFORTUMAB|ASG-?22ME|TIVDAK|TISOTUMAB|"
    r"ELAHERE|MIRVETUXIMAB|IMGN-?853|ADCETRIS|BRENTUXIMAB|SGN-?35|POLIVY|POLATUZUMAB|"
    r"ZYNLONTA|LONCASTUXIMAB|ADCT-?402|BESPONSA|INOTUZUMAB|MYLOTARG|GEMTUZUMAB|"
    r"BLENREP|BELANTAMAB|EMRELIS|TELISOTUZUMAB|ABBV-?399|AIDIXI|DISITAMAB|RC-?48)"
)


def fetch_bytes(url, retries=5):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:
            wait = 4 + attempt * 6
            print(f"  RETRY {attempt+1} {url.split('/')[-1]}: {e} (sleep {wait}s)", file=sys.stderr, flush=True)
            time.sleep(wait)
    raise RuntimeError(f"Failed downloading {url}")


def get_quarter_urls(max_quarter="2026Q3"):
    html = fetch_bytes(INDEX_URL).decode("utf-8", "replace")
    full = re.findall(r"https://fis\.fda\.gov/content/Exports/faers_ascii_\d{4}[qQ]\d\.zip", html)
    by_q = {}
    for u in full:
        m = re.search(r"faers_ascii_(\d{4})[qQ](\d)\.zip$", u)
        if m:
            key = f"{m.group(1)}Q{m.group(2)}"
            if key <= max_quarter:
                by_q[key] = u
    return sorted(by_q.items())


def find_member(zf, prefix):
    cands = []
    for name in zf.namelist():
        base = name.rsplit("/", 1)[-1].lower()
        if base.startswith(prefix) and base.endswith(".txt") and "stat" not in base:
            cands.append(name)
    if not cands:
        return None
    cands.sort(key=lambda s: (0 if "ascii" in s.lower() or "delet" in s.lower() else 1, len(s)))
    return cands[0]


def parse_header(line):
    if isinstance(line, bytes):
        if line.startswith(b"\xef\xbb\xbf"):
            line = line[3:]
        line = line.decode("latin-1", "replace")
    line = line.lstrip("\ufeff")
    if line.startswith("ï»¿"):
        line = line[3:]
    cols = [c.strip().lstrip("\ufeff").lower() for c in line.rstrip("\r\n").split("$")]
    return {c: i for i, c in enumerate(cols)}


def safe_int(s):
    s = s.strip()
    if not s:
        return 0
    try:
        return int(s)
    except ValueError:
        digits = re.sub(r"\D", "", s)
        return int(digits) if digits else 0


def match_adcs(text_upper):
    if not text_upper or not FAST_ADC_RE.search(text_upper):
        return []
    return [drug for drug, rx in ADC_REGEX.items() if rx.search(text_upper)]


def process_quarter(qkey, url):
    done_flag = OUT_DIR / f"{qkey}.done"
    if done_flag.exists():
        return qkey, "cached", 0, 0

    t0 = time.time()
    raw = fetch_bytes(url)
    mb = len(raw) / (1024 * 1024)

    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        del_mem = find_member(zf, "delete")
        del_ids = []
        if del_mem:
            with zf.open(del_mem) as f:
                for bline in f:
                    if bline.startswith(b"\xef\xbb\xbf"):
                        bline = bline[3:]
                    s = bline.decode("latin-1", "replace").strip().lstrip("\ufeff")
                    if s.isdigit():
                        del_ids.append(s)

        drug_mem = find_member(zf, "drug")
        if not drug_mem:
            raise RuntimeError(f"{qkey}: missing DRUG table")
        adc_pids = {}
        with zf.open(drug_mem) as f:
            hdr = parse_header(f.readline())
            i_pid = hdr["primaryid"]
            i_dn = hdr.get("drugname")
            i_ai = hdr.get("prod_ai")
            for bline in f:
                parts = bline.decode("latin-1", "replace").rstrip("\r\n").split("$")
                if len(parts) <= i_pid:
                    continue
                dn = parts[i_dn].upper() if (i_dn is not None and i_dn < len(parts)) else ""
                ai = parts[i_ai].upper() if (i_ai is not None and i_ai < len(parts)) else ""
                combo = f"{dn} | {ai}" if ai else dn
                hits = match_adcs(combo)
                if hits:
                    pid = parts[i_pid].strip()
                    s = adc_pids.get(pid)
                    if s is None:
                        adc_pids[pid] = set(hits)
                    else:
                        s.update(hits)

        adc_drugs = defaultdict(list)
        if adc_pids:
            with zf.open(drug_mem) as f:
                hdr = parse_header(f.readline())
                i_pid = hdr["primaryid"]
                i_seq = hdr.get("drug_seq")
                i_role = hdr.get("role_cod")
                i_dn = hdr.get("drugname")
                i_ai = hdr.get("prod_ai")
                for bline in f:
                    parts = bline.decode("latin-1", "replace").rstrip("\r\n").split("$")
                    if len(parts) <= i_pid:
                        continue
                    pid = parts[i_pid].strip()
                    if pid in adc_pids:
                        adc_drugs[pid].append({
                            "seq": parts[i_seq].strip() if (i_seq is not None and i_seq < len(parts)) else "",
                            "role": parts[i_role].strip().upper() if (i_role is not None and i_role < len(parts)) else "",
                            "dn": parts[i_dn].strip().upper() if (i_dn is not None and i_dn < len(parts)) else "",
                            "ai": parts[i_ai].strip().upper() if (i_ai is not None and i_ai < len(parts)) else "",
                        })

        reac_mem = find_member(zf, "reac")
        pid_mask = {}
        adc_pts = defaultdict(list)
        if reac_mem:
            with zf.open(reac_mem) as f:
                hdr = parse_header(f.readline())
                i_pid = hdr["primaryid"]
                i_pt = hdr["pt"]
                for bline in f:
                    parts = bline.decode("latin-1", "replace").rstrip("\r\n").split("$")
                    if len(parts) <= max(i_pid, i_pt):
                        continue
                    pid = parts[i_pid].strip()
                    pt = parts[i_pt].strip().upper()
                    bit = PT_TO_BIT.get(pt, 0)
                    if bit:
                        pid_i = safe_int(pid)
                        pid_mask[pid_i] = pid_mask.get(pid_i, 0) | bit
                    if pid in adc_pids:
                        adc_pts[pid].append(pt)

        outc_mem = find_member(zf, "outc")
        adc_outc = defaultdict(list)
        if outc_mem and adc_pids:
            with zf.open(outc_mem) as f:
                hdr = parse_header(f.readline())
                i_pid = hdr["primaryid"]
                i_oc = hdr.get("outc_cod", hdr.get("outc_code"))
                for bline in f:
                    parts = bline.decode("latin-1", "replace").rstrip("\r\n").split("$")
                    if len(parts) <= i_pid:
                        continue
                    pid = parts[i_pid].strip()
                    if pid in adc_pids and i_oc is not None and i_oc < len(parts):
                        oc = parts[i_oc].strip().upper()
                        if oc:
                            adc_outc[pid].append(oc)

        indi_mem = find_member(zf, "indi")
        adc_indi = defaultdict(list)
        if indi_mem and adc_pids:
            with zf.open(indi_mem) as f:
                hdr = parse_header(f.readline())
                i_pid = hdr["primaryid"]
                i_seq = hdr.get("indi_drug_seq")
                i_pt = hdr.get("indi_pt")
                for bline in f:
                    parts = bline.decode("latin-1", "replace").rstrip("\r\n").split("$")
                    if len(parts) <= i_pid:
                        continue
                    pid = parts[i_pid].strip()
                    if pid in adc_pids and i_pt is not None and i_pt < len(parts):
                        ipt = parts[i_pt].strip().upper()
                        seq = parts[i_seq].strip() if (i_seq is not None and i_seq < len(parts)) else ""
                        if ipt:
                            adc_indi[pid].append([seq, ipt])

        demo_mem = find_member(zf, "demo")
        if not demo_mem:
            raise RuntimeError(f"{qkey}: missing DEMO table")
        bg_cid, bg_fda, bg_ver, bg_pid, bg_msk = [], [], [], [], []
        adc_cases = []
        with zf.open(demo_mem) as f:
            hdr = parse_header(f.readline())
            i_pid = hdr["primaryid"]
            i_cid = hdr["caseid"]
            i_ver = hdr.get("caseversion")
            i_fda = hdr.get("fda_dt")
            i_ev = hdr.get("event_dt")
            i_age = hdr.get("age")
            i_agec = hdr.get("age_cod")
            i_sex = hdr.get("sex", hdr.get("gndr_cod"))
            i_occp = hdr.get("occp_cod")
            i_ctry = hdr.get("occr_country")
            i_rctry = hdr.get("reporter_country")
            i_rept = hdr.get("rept_cod")
            for bline in f:
                parts = bline.decode("latin-1", "replace").rstrip("\r\n").split("$")
                if len(parts) <= max(i_pid, i_cid):
                    continue
                pid = parts[i_pid].strip()
                cid = parts[i_cid].strip()
                pid_i = safe_int(pid)
                cid_i = safe_int(cid)
                ver_i = safe_int(parts[i_ver]) if (i_ver is not None and i_ver < len(parts)) else 0
                fda_i = safe_int(parts[i_fda]) if (i_fda is not None and i_fda < len(parts)) else 0
                msk = pid_mask.get(pid_i, 0)
                bg_cid.append(cid_i)
                bg_fda.append(fda_i)
                bg_ver.append(ver_i)
                bg_pid.append(pid_i)
                bg_msk.append(msk)
                if pid in adc_pids:
                    getf = lambda idx: parts[idx].strip() if (idx is not None and idx < len(parts)) else ""
                    adc_cases.append({
                        "q": qkey,
                        "pid": pid,
                        "cid": cid,
                        "ver": ver_i,
                        "fda_dt": fda_i,
                        "ev_dt": getf(i_ev),
                        "age": getf(i_age),
                        "age_cod": getf(i_agec).upper(),
                        "sex": getf(i_sex).upper(),
                        "occp": getf(i_occp).upper(),
                        "ctry": getf(i_ctry).upper() or getf(i_rctry).upper(),
                        "rept": getf(i_rept).upper(),
                        "matched": sorted(adc_pids[pid]),
                        "drugs": adc_drugs.get(pid, []),
                        "pts": adc_pts.get(pid, []),
                        "outc": adc_outc.get(pid, []),
                        "indi": adc_indi.get(pid, []),
                    })

    np.savez_compressed(
        OUT_DIR / f"{qkey}_bg.npz",
        cid=np.asarray(bg_cid, dtype=np.int64),
        fda=np.asarray(bg_fda, dtype=np.int32),
        ver=np.asarray(bg_ver, dtype=np.int32),
        pid=np.asarray(bg_pid, dtype=np.int64),
        msk=np.asarray(bg_msk, dtype=np.uint8),
    )
    with gzip.open(OUT_DIR / f"{qkey}_adcs.jsonl.gz", "wt", encoding="utf-8") as f:
        for rec in adc_cases:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    with open(OUT_DIR / f"{qkey}_del.txt", "w") as f:
        for d_id in del_ids:
            f.write(d_id + "\n")
    done_flag.write_text(f"{len(bg_cid)}\t{len(adc_cases)}\t{len(del_ids)}\n")
    dt = time.time() - t0
    return qkey, f"{mb:.1f}MB/{dt:.1f}s", len(bg_cid), len(adc_cases)


def main():
    quarters = get_quarter_urls()
    print(f"Found {len(quarters)} FAERS ASCII quarters: {quarters[0][0]} .. {quarters[-1][0]}", file=sys.stderr, flush=True)
    t_start = time.time()
    done, errors = 0, []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(process_quarter, qk, u): qk for qk, u in quarters}
        for fut in concurrent.futures.as_completed(futs):
            qk = futs[fut]
            try:
                qkey, info, n_bg, n_adc = fut.result()
                done += 1
                print(f"[{done:2d}/{len(quarters)}] {qkey}: {info} | demo={n_bg:,} | adc_cases={n_adc:,}", file=sys.stderr, flush=True)
            except Exception as e:
                print(f"ERROR in {qk}: {repr(e)}", file=sys.stderr, flush=True)
                errors.append((qk, repr(e)))
    if errors:
        sys.exit(1)
    print(f"ALL {done} QUARTERS COMPLETED in {time.time() - t_start:.1f}s", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
