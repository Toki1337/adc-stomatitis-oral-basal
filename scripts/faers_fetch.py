import json, time, urllib.request, urllib.parse, sys, os, gzip
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
BASE="https://api.fda.gov/drug/event.json"
DRUGS={
 "trastuzumab emtansine":["ADO-TRASTUZUMAB EMTANSINE","TRASTUZUMAB EMTANSINE","KADCYLA"],
 "trastuzumab deruxtecan":["FAM-TRASTUZUMAB DERUXTECAN-NXKI","TRASTUZUMAB DERUXTECAN","ENHERTU"],
 "datopotamab deruxtecan":["DATOPOTAMAB DERUXTECAN","DATOPOTAMAB DERUXTECAN-DLNK","DATROWAY"],
 "sacituzumab govitecan":["SACITUZUMAB GOVITECAN","SACITUZUMAB GOVITECAN-HZIY","TRODELVY"],
 "enfortumab vedotin":["ENFORTUMAB VEDOTIN","ENFORTUMAB VEDOTIN-EJFV","PADCEV"],
 "tisotumab vedotin":["TISOTUMAB VEDOTIN","TISOTUMAB VEDOTIN-TFTV","TIVDAK"],
 "mirvetuximab soravtansine":["MIRVETUXIMAB SORAVTANSINE","MIRVETUXIMAB SORAVTANSINE-GYNX","ELAHERE"],
 "telisotuzumab vedotin":["TELISOTUZUMAB VEDOTIN","TELISOTUZUMAB VEDOTIN-TLLV","EMRELIS"],
 "brentuximab vedotin":["BRENTUXIMAB VEDOTIN","ADCETRIS"],
 "polatuzumab vedotin":["POLATUZUMAB VEDOTIN","POLATUZUMAB VEDOTIN-PIIQ","POLIVY"],
 "loncastuximab tesirine":["LONCASTUXIMAB TESIRINE","LONCASTUXIMAB TESIRINE-LPYL","ZYNLONTA"],
 "inotuzumab ozogamicin":["INOTUZUMAB OZOGAMICIN","BESPONSA"],
 "gemtuzumab ozogamicin":["GEMTUZUMAB OZOGAMICIN","MYLOTARG"],
 "belantamab mafodotin":["BELANTAMAB MAFODOTIN","BELANTAMAB MAFODOTIN-BLMF","BLENREP"],
 "disitamab vedotin":["DISITAMAB VEDOTIN","AIDIXI"],
 "sacituzumab tirumotecan":["SACITUZUMAB TIRUMOTECAN"],
}
def drugq(names):
    parts=[]
    for n in names:
        for f in ["patient.drug.medicinalproduct","patient.drug.openfda.generic_name","patient.drug.openfda.brand_name","patient.drug.activesubstance.activesubstancename"]:
            parts.append('%s:"%s"'%(f,n))
    return "("+" OR ".join(parts)+")"
def get(params):
    url=BASE+"?"+urllib.parse.urlencode(params)
    for k in range(6):
        try:
            with urllib.request.urlopen(url,timeout=120) as r: return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code==404: return {"results":[]}
            time.sleep(3+k*5)
        except Exception: time.sleep(3+k*5)
    raise RuntimeError(url)
os.makedirs(f'{ROOT}/dl/faers/reports',exist_ok=True)
for d,names in DRUGS.items():
    fn=f"{ROOT}/dl/faers/reports/{d.replace(' ','_')}.jsonl.gz"
    if os.path.exists(fn): continue
    q=drugq(names); recs=[]; skip=0
    while True:
        js=get({"search":q,"limit":1000,"skip":skip})
        res=js.get("results",[])
        for r in res:
            recs.append({"id":r.get("safetyreportid"),"date":r.get("receivedate"),"serious":r.get("serious"),
              "country":r.get("occurcountry") or r.get("primarysourcecountry"),
              "qual":(r.get("primarysource") or {}).get("qualification"),
              "drugs":[{"name":x.get("medicinalproduct"),"char":x.get("drugcharacterization"),
                        "gen":(x.get("openfda") or {}).get("generic_name"),
                        "ind":x.get("drugindication")} for x in r["patient"].get("drug",[])],
              "pts":[y.get("reactionmeddrapt") for y in r["patient"].get("reaction",[])]})
        skip+=len(res)
        if len(res)<1000 or skip>=25000: break
        time.sleep(0.3)
    with gzip.open(fn,"wt") as f:
        for r in recs: f.write(json.dumps(r)+"\n")
    print(d,len(recs),file=sys.stderr,flush=True)
