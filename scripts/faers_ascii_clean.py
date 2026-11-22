"""Deduplicate extracted FDA FAERS Quarterly ASCII records (2012Q4-2026Q3) per FDA
and READUS-PV guidelines (official DELETE removal + CASEID latest-version retention),
apply exact drug-role and antineoplastic exclusion filters, and write processed summary
tables to results/faers_counts.tsv, results/faers_background.tsv, and
results/faers_readus_baseline.tsv.
"""

import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXT_DIR = ROOT / "dl" / "faers" / "ascii_extract"
RES_DIR = ROOT / "results"

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

CORE_15 = [
    "trastuzumab emtansine",
    "trastuzumab deruxtecan",
    "datopotamab deruxtecan",
    "sacituzumab govitecan",
    "enfortumab vedotin",
    "tisotumab vedotin",
    "mirvetuximab soravtansine",
    "brentuximab vedotin",
    "polatuzumab vedotin",
    "loncastuximab tesirine",
    "inotuzumab ozogamicin",
    "gemtuzumab ozogamicin",
    "belantamab mafodotin",
    "telisotuzumab vedotin",
    "disitamab vedotin",
]

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

# Comprehensive ATC L01 antineoplastic, anti-HER2, and checkpoint-inhibitor exclusion list
ANTINEOPLASTIC_TERMS = [
    "CYCLOPHOSPHAMIDE", "ENDOXAN", "CYTOXAN",
    "DOXORUBICIN", "ADRIAMYCIN", "CAELYX", "DOXIL", "MYOCET",
    "EPIRUBICIN", "DAUNORUBICIN", "IDARUBICIN", "MITOXANTRONE", "LIPOSOMAL DOXORUBICIN",
    "VINCRISTINE", "ONCOVIN", "VINBLASTINE", "VINORELBINE", "NAVELBINE", "VINDESINE",
    "RITUXIMAB", "MABTHERA", "RITUXAN", "TRUXIMA", "RUXIENCE", "RIABNI",
    "OBINUTUZUMAB", "GAZYVA", "GAZYVARO", "OFATUMUMAB",
    "BENDAMUSTINE", "TREANDA", "BENDEKA", "LEVACT",
    "CYTARABINE", "ARA-C", "ARAC", "FLUDARABINE", "CLADRIBINE", "CLOFARABINE", "NELARABINE",
    "DECITABINE", "AZACITIDINE", "VIDAZA", "DACOGEN", "INQOVI",
    "METHOTREXATE", "PEMETREXED", "ALIMTA", "PRALATREXATE", "RALTITREXED",
    "ETOPOSIDE", "VEPESID", "TENIPOSIDE",
    "CISPLATIN", "CARBOPLATIN", "OXALIPLATIN", "NEDAPLATIN", "LOBAPLATIN",
    "FLUOROURACIL", "5-FU", "5FU", "CAPECITABINE", "XELODA", "TEGAFUR", "TRIFLURIDINE", "LONSURF",
    "GEMCITABINE", "GEMZAR",
    "PACLITAXEL", "TAXOL", "ABRAXANE", "DOCETAXEL", "TAXOTERE", "CABAZITAXEL", "JEVTANA",
    "ERIBULIN", "HALAVEN", "IXABEPILONE", "TRABECTEDIN", "LURBINECTEDIN",
    "PEMBROLIZUMAB", "KEYTRUDA", "NIVOLUMAB", "OPDIVO", "ATEZOLIZUMAB", "TECENTRIQ",
    "DURVALUMAB", "IMFINZI", "AVELUMAB", "BAVENCIO", "IPILIMUMAB", "YERVOY",
    "TREMELIMUMAB", "CEMIPLIMAB", "DOSTARLIMAB", "TORIPALIMAB", "SINTILIMAB", "CAMRELIZUMAB", "TISLELIZUMAB", "CADONILIMAB",
    "DACARBAZINE", "TEMOZOLOMIDE", "TEMODAR", "BLEOMYCIN", "MELPHALAN", "BUSULFAN", "THIOTEPA",
    "CARMUSTINE", "LOMUSTINE", "PROCARBAZINE", "IFOSFAMIDE", "CHLORAMBUCIL", "MITOMYCIN",
    "HYDROXYUREA", "HYDROXYCARBAMIDE", "ASPARAGINASE", "PEGASPARGASE",
    "LENALIDOMIDE", "REVLIMID", "POMALIDOMIDE", "POMALYST", "IMNOVID", "THALIDOMIDE",
    "BORTEZOMIB", "VELCADE", "CARFILZOMIB", "KYPROLIS", "IXAZOMIB", "NINLARO",
    "DARATUMUMAB", "DARZALEX", "ISATUXIMAB", "SARCLISA", "ELOTUZUMAB",
    "VENETOCLAX", "VENCLEXTA", "VENCLYXTO",
    "BLINATUMOMAB", "BLINCYTO", "TAFASITAMAB", "GLOFITAMAB", "MOSUNETUZUMAB",
    "EPCORITAMAB", "TECLISTAMAB", "ELRANATAMAB", "TALQUETAMAB",
    "IBRUTINIB", "IMBRUVICA", "ACALABRUTINIB", "CALQUENCE", "ZANUBRUTINIB", "BRUKINSA",
    "PIRTOBRUTINIB", "IDELALISIB", "COPANLISIB", "DUVELISIB",
    "IRINOTECAN", "CAMPTOSAR", "ONIVYDE", "TOPOTECAN",
    "EVEROLIMUS", "AFINITOR", "TEMSIROLIMUS",
    "BEVACIZUMAB", "AVASTIN", "MVASI", "ZIRABEV", "RAMUCIRUMAB", "CYRAMZA",
    "CETUXIMAB", "ERBITUX", "PANITUMUMAB", "VECTIBIX",
    "TRASTUZUMAB", "HERCEPTIN", "KANJINTI", "OGIVRI", "TRAZIMERA", "HERZUMA", "ONTRUZANT", "PHESGO",
    "PERTUZUMAB", "PERJETA", "MARGETUXIMAB",
    "OLAPARIB", "LYNPARZA", "NIRAPARIB", "ZEJULA", "RUCAPARIB", "TALAZOPARIB",
    "PALBOCICLIB", "IBRANCE", "ABEMACICLIB", "VERZENIO", "VERZENIOS", "RIBOCICLIB", "KISQALI",
    "TUCATINIB", "TUKYSA", "LAPATINIB", "TYKERB", "TYVERB", "NERATINIB", "NERLYNX", "PYROTINIB",
    "ALPELISIB", "PIQRAY", "CAPIVASERTIB", "TRUQAP",
    "OSIMERTINIB", "TAGRISSO", "ERLOTINIB", "GEFITINIB", "AFATINIB", "DACOMITINIB", "AMIVANTAMAB", "LAZERTINIB",
    "IMATINIB", "GLEEVEC", "GLIVEC", "DASATINIB", "SPRYCEL", "NILOTINIB", "TASIGNA", "BOSUTINIB", "PONATINIB", "ASCIMINIB",
    "MIDOSTAURIN", "GILTERITINIB", "XOSPATA", "QUIZARTINIB", "IVOSIDENIB", "ENASIDENIB", "GLASDEGIB",
    "SORAFENIB", "LENVATINIB", "CABOZANTINIB", "SUNITINIB", "PAZOPANIB", "AXITINIB", "REGORAFENIB",
    "TIVOZANIB", "ANLOTINIB", "APATINIB", "FRUQUINTINIB", "ERDAFITINIB",
    "CRIZOTINIB", "ALECTINIB", "BRIGATINIB", "LORLATINIB", "ENTRECTINIB", "LAROTRECTINIB",
    "SELPERCATINIB", "PRALSETINIB", "CAPMATINIB", "TEPOTINIB", "SOTORASIB", "ADAGRASIB",
    "DABRAFENIB", "TRAMETINIB", "VEMURAFENIB", "COBIMETINIB", "ENCORAFENIB", "BINIMETINIB",
    "RUXOLITINIB", "PANOBINOSTAT", "VORINOSTAT", "ROMIDEPSIN", "BELINOSTAT", "SELINEXOR",
    "ARSENIC TRIOXIDE", "TRETINOIN", "ATRA", "MERCAPTOPURINE", "THIOGUANINE",
    "CAR-T", "TISAGENLECLEUCEL", "AXICABTAGENE", "BREXUCABTAGENE", "LISOCABTAGENE", "IDECABTAGENE", "CILTACABTAGENE",
    "MOXETUMOMAB", "PATRITUMAB", "IFINATAMAB", "SIGVOTATUG", "SACITUZUMAB TIRUMOTECAN",
    "TELISOTUZUMAB ADIZUTECAN",
]
ANTINEO_RE = re.compile(
    r"\b(?:" + "|".join(re.escape(t) for t in sorted(set(ANTINEOPLASTIC_TERMS), key=len, reverse=True)) + r")\b"
)
HNSCC_INDI_RE = re.compile(
    r"(?:HEAD AND NECK|ORAL CAVITY|ORAL CANCER|ORAL SQUAMOUS|TONGUE|OROPHARYNG|HYPOPHARYNG|"
    r"NASOPHARYNG|LARYNG|PHARYNG|BUCCAL|GINGIVA|PALATE|LIP CARCINOMA|SALIVARY|ESOPHAG|OESOPHAG|"
    r"RADIOTHERAPY|RADIATION)"
)


def parse_age_years(age_str, cod):
    if not age_str:
        return np.nan
    try:
        v = float(age_str)
    except ValueError:
        return np.nan
    if v <= 0 or v > 12000:
        return np.nan
    cod = (cod or "YR").upper()
    scale = {"YR": 1.0, "DEC": 10.0, "MON": 1.0 / 12.0, "WK": 1.0 / 52.18, "DY": 1.0 / 365.25, "HR": 1.0 / 8766.0}
    if cod not in scale:
        return np.nan
    yrs = v * scale[cod]
    return yrs if 0 <= yrs <= 115 else np.nan


def match_drug_row(dn, ai):
    if "PLACEBO" in dn:
        return []
    hits = []
    for d, rx in ADC_REGEX.items():
        if ai:
            if rx.search(ai):
                hits.append(d)
            elif rx.search(dn) and not ANTINEO_RE.search(ai):
                hits.append(d)
        else:
            if rx.search(dn):
                hits.append(d)
    return hits


def is_other_antineoplastic(dn, ai, target_adc):
    combo = f"{dn} | {ai}" if ai else dn
    cleaned = ADC_REGEX[target_adc].sub(" ", combo)
    for d, rx in ADC_REGEX.items():
        if d != target_adc and rx.search(cleaned):
            return True
    return bool(ANTINEO_RE.search(cleaned))


def run_clean():
    if not EXT_DIR.exists() or not list(EXT_DIR.glob("*_bg.npz")):
        raise FileNotFoundError(f"Missing extracted quarterly FAERS files in {EXT_DIR}; run faers_ascii_extract.py first.")

    # 1. Deleted caseids
    del_files = sorted(EXT_DIR.glob("*_del.txt"))
    deleted_cids = set()
    for df in del_files:
        for line in df.read_text().splitlines():
            s = line.strip()
            if s.isdigit():
                deleted_cids.add(int(s))

    # 2. Whole-FAERS background deduplication
    bg_files = sorted(EXT_DIR.glob("*_bg.npz"))
    quarters = sorted(f.name.replace("_bg.npz", "") for f in bg_files)
    cids_l, fdas_l, vers_l, pids_l, msks_l = [], [], [], [], []
    raw_bg_total = 0
    for bf in bg_files:
        z = np.load(bf)
        cids_l.append(z["cid"])
        fdas_l.append(z["fda"])
        vers_l.append(z["ver"])
        pids_l.append(z["pid"])
        msks_l.append(z["msk"])
        raw_bg_total += len(z["cid"])
    cid_all = np.concatenate(cids_l)
    fda_all = np.concatenate(fdas_l)
    ver_all = np.concatenate(vers_l)
    pid_all = np.concatenate(pids_l)
    msk_all = np.concatenate(msks_l)
    del cids_l, fdas_l, vers_l, pids_l, msks_l

    if deleted_cids:
        del_arr = np.fromiter(deleted_cids, dtype=np.int64)
        keep = ~np.isin(cid_all, del_arr)
        cid_all, fda_all, ver_all, pid_all, msk_all = cid_all[keep], fda_all[keep], ver_all[keep], pid_all[keep], msk_all[keep]
    after_del_bg = int(len(cid_all))

    order = np.lexsort((pid_all, ver_all, fda_all, cid_all))
    cid_s = cid_all[order]
    msk_s = msk_all[order]
    del cid_all, fda_all, ver_all, pid_all, msk_all, order
    is_last = np.empty(len(cid_s), dtype=bool)
    is_last[:-1] = cid_s[:-1] != cid_s[1:]
    is_last[-1] = True
    dedup_msk = msk_s[is_last]
    bg_tot = int(len(dedup_msk))

    bg_rows = [
        {"metric": "quarter_start", "value": quarters[0]},
        {"metric": "quarter_end", "value": quarters[-1]},
        {"metric": "n_quarters", "value": len(quarters)},
        {"metric": "n_deleted_caseids", "value": len(deleted_cids)},
        {"metric": "raw_faers_rows", "value": raw_bg_total},
        {"metric": "after_delete_rows", "value": after_del_bg},
        {"metric": "dedup_total_cases", "value": bg_tot},
    ]
    for bit_idx, (out, _) in enumerate(OUT_GROUPS):
        bg_rows.append({"metric": f"dedup_{out}", "value": int(((dedup_msk & (1 << bit_idx)) != 0).sum())})
    pd.DataFrame(bg_rows).to_csv(RES_DIR / "faers_background.tsv", sep="\t", index=False)
    del cid_s, msk_s, is_last, dedup_msk

    # 3. ADC case deduplication (strictly across CORE_15)
    adc_files = sorted(EXT_DIR.glob("*_adcs.jsonl.gz"))
    raw_adc_rows = 0
    del_adc_rows = 0
    best_by_cid = {}
    for af in adc_files:
        with gzip.open(af, "rt", encoding="utf-8") as f:
            for line in f:
                rec = json.loads(line)
                # Ensure record actually matches at least one CORE_15 ADC cleanly
                matched_core = set()
                for dr in rec["drugs"]:
                    for h in match_drug_row(dr["dn"], dr["ai"]):
                        if h in CORE_15:
                            matched_core.add(h)
                if not matched_core:
                    continue
                raw_adc_rows += 1
                cid_i = int(rec["cid"])
                if cid_i in deleted_cids:
                    del_adc_rows += 1
                    continue
                key = (int(rec["fda_dt"]), int(rec["ver"]), int(rec["pid"]))
                prev = best_by_cid.get(cid_i)
                if prev is None or key > prev[0]:
                    best_by_cid[cid_i] = (key, rec)

    dedup_cases = [v[1] for v in best_by_cid.values()]

    # Secondary exact demographic/event/drug/PT duplicate check
    seen_sig = set()
    sec_dup_cids = set()
    for rec in dedup_cases:
        if rec["ev_dt"] and rec["age"] and rec["sex"] in ("M", "F"):
            sig = (
                rec["ev_dt"], rec["age"], rec["age_cod"], rec["sex"], rec["ctry"],
                tuple(sorted(set(rec["pts"]))),
                tuple(sorted(set(dr["dn"] for dr in rec["drugs"]))),
            )
            if sig in seen_sig:
                sec_dup_cids.add(rec["cid"])
            else:
                seen_sig.add(sig)

    # Cohorts: primary = sole_suspect_noco (exact ANALYSIS_PLAN.md) + 6 sensitivity cohorts
    cohorts = [
        "sole_suspect_noco",          # Primary: dedup + sole suspect (PS, no SS) + no co-antineoplastic
        "ps_noco",                    # Sensitivity 1: dedup + PS (SS allowed) + no co-antineoplastic
        "sole_suspect_noco_hcp",      # Sensitivity 2: primary restricted to HCP reporters (MD/PH/OT)
        "sole_suspect_noco_serious",  # Sensitivity 3: primary restricted to serious cases (OUTC present)
        "sole_suspect_noco_nohnscc",  # Sensitivity 4: primary excluding HNSCC/oral-cancer/radiation indications
        "sole_suspect_noco_secdedup", # Sensitivity 5: primary + secondary demographic/event/drug/PT dedup
        "ps_all",                     # Sensitivity 6: all dedup PS cases without co-antineoplastic exclusion
    ]
    cnt = {c: {d: Counter() for d in CORE_15} for c in cohorts}

    # READUS-PV baseline tracking for primary cohort (sole_suspect_noco) and attrition funnel
    n_ps_any_unique = 0
    n_ps_noco_unique = 0
    n_sole_noco_unique = 0
    ds = {"ages": [], "sex": Counter(), "occp": Counter(), "ser": 0, "death": 0, "hosp": 0, "ctry": Counter(), "q": Counter(), "n": 0}

    for rec in dedup_cases:
        pt_set = set(rec["pts"])
        out_hits = [out for out, pts in OUT_GROUPS if pt_set & pts]
        is_hcp = rec["occp"] in ("MD", "PH", "OT")
        is_ser = len(rec["outc"]) > 0
        is_death = "DE" in rec["outc"]
        is_hosp = "HO" in rec["outc"]
        indi_text = " | ".join(x[1] for x in rec["indi"])
        is_hnscc = bool(HNSCC_INDI_RE.search(indi_text))
        is_sec_dup = rec["cid"] in sec_dup_cids

        drug_roles_by_adc = defaultdict(set)
        for dr in rec["drugs"]:
            for h in match_drug_row(dr["dn"], dr["ai"]):
                drug_roles_by_adc[h].add(dr["role"])

        ps_adcs = [adc for adc, roles in drug_roles_by_adc.items() if "PS" in roles]
        if ps_adcs:
            n_ps_any_unique += 1
        any_ps_noco = False
        any_sole_noco = False

        for adc in ps_adcs:
            has_other_antineo = False
            has_other_suspect = False
            for dr in rec["drugs"]:
                row_adcs = match_drug_row(dr["dn"], dr["ai"])
                if is_other_antineoplastic(dr["dn"], dr["ai"], adc):
                    has_other_antineo = True
                if (adc not in row_adcs) and dr["role"] in ("PS", "SS"):
                    has_other_suspect = True

            is_ps_noco = not has_other_antineo
            is_sole_noco = (not has_other_suspect) and (not has_other_antineo)
            if is_ps_noco:
                any_ps_noco = True
            if is_sole_noco:
                any_sole_noco = True

            flags = {
                "sole_suspect_noco": is_sole_noco,
                "ps_noco": is_ps_noco,
                "sole_suspect_noco_hcp": is_sole_noco and is_hcp,
                "sole_suspect_noco_serious": is_sole_noco and is_ser,
                "sole_suspect_noco_nohnscc": is_sole_noco and (not is_hnscc),
                "sole_suspect_noco_secdedup": is_sole_noco and (not is_sec_dup),
                "ps_all": True,
            }
            for cname, ok in flags.items():
                if ok:
                    cnt[cname][adc]["n_reports"] += 1
                    for out in out_hits:
                        cnt[cname][adc]["n_" + out] += 1

        if any_ps_noco:
            n_ps_noco_unique += 1
        if any_sole_noco:
            n_sole_noco_unique += 1
            ds["n"] += 1
            ay = parse_age_years(rec["age"], rec["age_cod"])
            if np.isfinite(ay):
                ds["ages"].append(ay)
            ds["sex"][rec["sex"] or "UNK"] += 1
            ds["occp"][rec["occp"] or "UNK"] += 1
            if is_ser:
                ds["ser"] += 1
            if is_death:
                ds["death"] += 1
            if is_hosp:
                ds["hosp"] += 1
            ds["ctry"][rec["ctry"] or "UNK"] += 1
            ds["q"][rec["q"]] += 1

    # Save results/faers_counts.tsv
    crows = []
    for cname in cohorts:
        for d in CORE_15:
            row = {"drug": d, "set": cname, "n_reports": cnt[cname][d]["n_reports"]}
            for out, _ in OUT_GROUPS:
                row["n_" + out] = cnt[cname][d]["n_" + out]
            crows.append(row)
    pd.DataFrame(crows).to_csv(RES_DIR / "faers_counts.tsv", sep="\t", index=False)

    # Save results/faers_readus_baseline.tsv
    ages = np.asarray(ds["ages"], dtype=float)
    n = ds["n"]
    hcp_n = ds["occp"]["MD"] + ds["occp"]["PH"] + ds["occp"]["OT"]
    qs = sorted(ds["q"].keys())
    top_ctry = "; ".join(f"{k}: {v:,} ({100*v/n:.1f}%)" for k, v in ds["ctry"].most_common(6))
    brows = [
        ("Database & Period", "FDA FAERS Quarterly ASCII Extracts (56 quarters)", f"{quarters[0]}–{quarters[-1]}"),
        ("Whole-FAERS Deduplication", "Raw DEMO rows -> after official DELETE -> unique CASEID", f"{raw_bg_total:,} -> {after_del_bg:,} -> {bg_tot:,} (removed {len(deleted_cids):,} deleted case IDs)"),
        ("ADC Case Attrition (Step 1)", "Raw FAERS records mentioning >=1 of 15 marketed ADCs", f"{raw_adc_rows:,}"),
        ("ADC Case Attrition (Step 2)", "After removing FDA-deleted case IDs & CASEID deduplication", f"{len(dedup_cases):,} unique cases ({del_adc_rows:,} deleted rows removed)"),
        ("ADC Case Attrition (Step 3)", "Primary Suspect (ROLE_COD = PS) ADC cases", f"{n_ps_any_unique:,} unique cases"),
        ("ADC Case Attrition (Step 4)", "Primary Suspect + no co-reported antineoplastic/anti-HER2 drug", f"{n_ps_noco_unique:,} unique cases"),
        ("ADC Case Attrition (Step 5, Primary)", "Sole Suspect (PS with zero SS) + no co-reported antineoplastic (14 FDA-approved ADCs >=50 reports; disitamab vedotin 0 PS)", f"{n_sole_noco_unique:,} unique cases"),
        ("Secondary Duplicate Check", "Exact demographic/date/drug/PT duplicates flagged in sensitivity", f"{len(sec_dup_cids):,} cases among all deduplicated ADC records"),
        ("Primary Cohort Age", f"Median (IQR), years (n = {len(ages):,} with reported age)", f"{np.median(ages):.1f} ({np.percentile(ages, 25):.1f}–{np.percentile(ages, 75):.1f})"),
        ("Primary Cohort Sex", "Female / Male / Unknown", f"{ds['sex']['F']:,} ({100*ds['sex']['F']/n:.1f}%) / {ds['sex']['M']:,} ({100*ds['sex']['M']/n:.1f}%) / {n - ds['sex']['F'] - ds['sex']['M']:,} ({100*(n - ds['sex']['F'] - ds['sex']['M'])/n:.1f}%)"),
        ("Primary Cohort Reporter", "Healthcare professional (MD/PH/OT) / Consumer (CN) / Unknown", f"{hcp_n:,} ({100*hcp_n/n:.1f}%) / {ds['occp']['CN']:,} ({100*ds['occp']['CN']/n:.1f}%) / {n - hcp_n - ds['occp']['CN']:,} ({100*(n - hcp_n - ds['occp']['CN'])/n:.1f}%)"),
        ("Primary Cohort Seriousness", "Any serious outcome (OUTC) / Hospitalization / Death", f"{ds['ser']:,} ({100*ds['ser']/n:.1f}%) / {ds['hosp']:,} ({100*ds['hosp']/n:.1f}%) / {ds['death']:,} ({100*ds['death']/n:.1f}%)"),
        ("Primary Cohort Top Countries", "Top 6 reporting countries", top_ctry),
        ("Primary Cohort Quarter Span", "First to last quarter with primary cases", f"{qs[0]}–{qs[-1]}"),
    ]
    pd.DataFrame(brows, columns=["category", "item", "value"]).to_csv(RES_DIR / "faers_readus_baseline.tsv", sep="\t", index=False)
    print(f"Wrote faers_counts.tsv ({len(crows)} rows), faers_background.tsv, faers_readus_baseline.tsv (primary N={n:,}).")


if __name__ == "__main__":
    run_clean()
