import json, time, urllib.request, urllib.parse, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
BASE="https://api.fda.gov/drug/event.json"
ORAL=["STOMATITIS","MUCOSAL INFLAMMATION","MOUTH ULCERATION","APHTHOUS ULCER","GLOSSITIS","ORAL MUCOSAL ERYTHEMA","ORAL MUCOSAL BLISTERING"]
def total(search):
    url=BASE+"?"+urllib.parse.urlencode({"search":search,"limit":1})
    for k in range(5):
        try:
            with urllib.request.urlopen(url,timeout=60) as r:
                return json.load(r)["meta"]["results"]["total"]
        except urllib.error.HTTPError as e:
            if e.code==404: return 0
            time.sleep(2+k*3)
        except Exception:
            time.sleep(2+k*3)
    return None
ev='patient.reaction.reactionmeddrapt.exact:('+"+OR+".join('"%s"'%p for p in ORAL)+')'
ev=ev.replace("+"," ")
# drug name variants: (field, value)
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
 "patritumab deruxtecan":["PATRITUMAB DERUXTECAN"],
 "sacituzumab tirumotecan":["SACITUZUMAB TIRUMOTECAN"],
 "ifinatamab deruxtecan":["IFINATAMAB DERUXTECAN"],
 "tusamitamab ravtansine":["TUSAMITAMAB RAVTANSINE"],
 "rovalpituzumab tesirine":["ROVALPITUZUMAB TESIRINE"],
 "depatuxizumab mafodotin":["DEPATUXIZUMAB MAFODOTIN"],
 "trastuzumab duocarmazine":["TRASTUZUMAB DUOCARMAZINE"],
 "upifitamab rilsodotin":["UPIFITAMAB RILSODOTIN"],
 "anetumab ravtansine":["ANETUMAB RAVTANSINE"],
}
def drugq(names, suspect=True):
    parts=[]
    for n in names:
        parts.append('patient.drug.medicinalproduct:"%s"'%n)
        parts.append('patient.drug.openfda.generic_name:"%s"'%n)
        parts.append('patient.drug.openfda.brand_name:"%s"'%n)
        parts.append('patient.drug.activesubstance.activesubstancename:"%s"'%n)
    return "("+" OR ".join(parts)+")"
N_all=total('receivedate:[19000101 TO 20991231]')
N_ev=total(ev)
print("ALL",N_all,"EV",N_ev,file=sys.stderr)
out=open(f"{ROOT}/dl/faers/faers_oral_counts.tsv","w")
out.write("drug\tn_reports\tn_oral\tN_all\tN_oral_all\n")
for d,names in DRUGS.items():
    q=drugq(names)
    a_b=total(q); a=total(q+" AND "+ev)
    print(d,a_b,a,file=sys.stderr)
    out.write(f"{d}\t{a_b}\t{a}\t{N_all}\t{N_ev}\n"); out.flush()
    time.sleep(0.5)
