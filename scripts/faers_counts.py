import json, time, urllib.request, urllib.parse, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
BASE="https://api.fda.gov/drug/event.json"
def total(search):
    url=BASE+"?"+urllib.parse.urlencode({"search":search,"limit":1})
    for k in range(6):
        try:
            with urllib.request.urlopen(url,timeout=60) as r: return json.load(r)["meta"]["results"]["total"]
        except urllib.error.HTTPError as e:
            if e.code==404: return 0
            print('HTTP',e.code,file=sys.stderr); time.sleep(3+k*5)
        except Exception as e: time.sleep(3+k*5)
    return None
def orq(field,vals): return "("+" OR ".join('%s:"%s"'%(field,v) for v in vals)+")"
OUT={
 "oral":["STOMATITIS","MUCOSAL INFLAMMATION","MOUTH ULCERATION","APHTHOUS ULCER","GLOSSITIS","ORAL MUCOSAL ERYTHEMA","ORAL MUCOSAL BLISTERING"],
 "neuropathy":["PERIPHERAL SENSORY NEUROPATHY","NEUROPATHY PERIPHERAL","PERIPHERAL MOTOR NEUROPATHY","POLYNEUROPATHY"],
 "neutropenia":["NEUTROPENIA","FEBRILE NEUTROPENIA","NEUTROPHIL COUNT DECREASED"],
 "skin":["RASH","RASH MACULO-PAPULAR","RASH PRURITIC","RASH ERYTHEMATOUS","DERMATITIS","SKIN EXFOLIATION","STEVENS-JOHNSON SYNDROME","TOXIC EPIDERMAL NECROLYSIS"],
 "ocular":["KERATOPATHY","KERATITIS","VISION BLURRED","DRY EYE","CONJUNCTIVITIS","VISUAL IMPAIRMENT","CORNEAL DISORDER"],
 "ild":["INTERSTITIAL LUNG DISEASE","PNEUMONITIS"],
 "diarrhea":["DIARRHOEA"],
}
EXCL_GEN=["CYCLOPHOSPHAMIDE","DOXORUBICIN","VINCRISTINE","RITUXIMAB","BENDAMUSTINE","CYTARABINE","DAUNORUBICIN","IDARUBICIN","FLUDARABINE","METHOTREXATE","ETOPOSIDE","CISPLATIN","CARBOPLATIN","OXALIPLATIN","FLUOROURACIL","CAPECITABINE","GEMCITABINE","PACLITAXEL","DOCETAXEL","PEMBROLIZUMAB","NIVOLUMAB","ATEZOLIZUMAB","DURVALUMAB","DACARBAZINE","VINBLASTINE","BLEOMYCIN","MELPHALAN","BUSULFAN","THIOTEPA","LENALIDOMIDE","POMALIDOMIDE","BORTEZOMIB","DARATUMUMAB","AZACITIDINE","VENETOCLAX","IFOSFAMIDE","BLINATUMOMAB","OBINUTUZUMAB","IBRUTINIB","IRINOTECAN","TOPOTECAN","EVEROLIMUS","BEVACIZUMAB","CETUXIMAB","OLAPARIB","PALBOCICLIB","ABEMACICLIB","RIBOCICLIB","TUCATINIB","LAPATINIB","NERATINIB","EPIRUBICIN","ERIBULIN","VINORELBINE","CARFILZOMIB","ISATUXIMAB","TAFASITAMAB","GLOFITAMAB","MOSUNETUZUMAB","EPCORITAMAB","CLOFARABINE","MITOXANTRONE","PEMETREXED","IPILIMUMAB","CYTARABINE AND DAUNORUBICIN"]
EXCL_EXACT=["TRASTUZUMAB","PERTUZUMAB","TRASTUZUMAB AND HYALURONIDASE-OYSK","PERTUZUMAB, TRASTUZUMAB, AND HYALURONIDASE-ZZXF"]
excl="("+orq("patient.drug.openfda.generic_name",EXCL_GEN)[1:-1]+" OR "+orq("patient.drug.openfda.generic_name.exact",EXCL_EXACT)[1:-1]+")"
DRUGS={
 "trastuzumab emtansine":["ADO-TRASTUZUMAB EMTANSINE","KADCYLA"],
 "trastuzumab deruxtecan":["FAM-TRASTUZUMAB DERUXTECAN-NXKI","TRASTUZUMAB DERUXTECAN","ENHERTU"],
 "datopotamab deruxtecan":["DATOPOTAMAB DERUXTECAN","DATOPOTAMAB DERUXTECAN-DLNK","DATROWAY"],
 "sacituzumab govitecan":["SACITUZUMAB GOVITECAN","SACITUZUMAB GOVITECAN-HZIY","TRODELVY"],
 "enfortumab vedotin":["ENFORTUMAB VEDOTIN","ENFORTUMAB VEDOTIN-EJFV","PADCEV"],
 "tisotumab vedotin":["TISOTUMAB VEDOTIN","TISOTUMAB VEDOTIN-TFTV","TIVDAK"],
 "mirvetuximab soravtansine":["MIRVETUXIMAB SORAVTANSINE","MIRVETUXIMAB SORAVTANSINE-GYNX","ELAHERE"],
 "brentuximab vedotin":["BRENTUXIMAB VEDOTIN","ADCETRIS"],
 "polatuzumab vedotin":["POLATUZUMAB VEDOTIN","POLATUZUMAB VEDOTIN-PIIQ","POLIVY"],
 "loncastuximab tesirine":["LONCASTUXIMAB TESIRINE","LONCASTUXIMAB TESIRINE-LPYL","ZYNLONTA"],
 "inotuzumab ozogamicin":["INOTUZUMAB OZOGAMICIN","BESPONSA"],
 "gemtuzumab ozogamicin":["GEMTUZUMAB OZOGAMICIN","MYLOTARG"],
 "belantamab mafodotin":["BELANTAMAB MAFODOTIN","BELANTAMAB MAFODOTIN-BLMF","BLENREP"],
 "telisotuzumab vedotin":["TELISOTUZUMAB VEDOTIN","TELISOTUZUMAB VEDOTIN-TLLV","EMRELIS"],
 "disitamab vedotin":["DISITAMAB VEDOTIN","AIDIXI"],
}
def drugq(names):
    parts=[]
    for n in names:
        for f in ["patient.drug.medicinalproduct","patient.drug.openfda.generic_name.exact","patient.drug.openfda.brand_name.exact","patient.drug.activesubstance.activesubstancename"]:
            parts.append('%s:"%s"'%(f,n))
    return "("+" OR ".join(parts)+")"
out=open(f"{ROOT}/results/faers_counts.tsv","w")
out.write("drug\tset\tn_reports\t"+"\t".join("n_"+k for k in OUT)+"\n")
for d,names in DRUGS.items():
    base=drugq(names)
    for sname,q in [("all",base),("noco",base+" AND NOT "+excl)]:
        row=[total(q)]
        for k,pts in OUT.items():
            row.append(total(q+" AND "+orq("patient.reaction.reactionmeddrapt.exact",pts))); time.sleep(0.25)
        out.write(f"{d}\t{sname}\t"+"\t".join(map(str,row))+"\n"); out.flush()
        print(d,sname,row,file=sys.stderr,flush=True)
