# -*- coding: utf-8 -*-
"""36_en_refs_renumber.py — 英文手稿引文按首现顺序重排 + 生成 References 节
输入: results/paper2_manuscript_EN_v1.0.md  ([#n] / [#a-#b] / [#a,#b] 令牌)
输出: 原文件替换令牌为 [n,m,...]，References 占位替换为正式编号列表
"""
import re

BASE = r"C:\Users\Administrator\WorkBuddy\2026-09-27-09-39-22\iua-exosome-bioinfo"
MD = BASE + r"\results\paper2_manuscript_EN_v1.0.md"

REFS = {
 1:"Dixon SJ, Lemberg KM, Lamprecht MR, et al. Ferroptosis: an iron-dependent form of nonapoptotic cell death. Cell. 2012;149(5):1060-1072. PMID 22632970",
 2:"Stockwell BR, Friedmann Angeli JP, Bayir H, et al. Ferroptosis: a regulated cell death nexus linking metabolism, redox biology, and disease. Cell. 2017;171(2):273-285. PMID 28985560",
 3:"Jiang X, Stockwell BR, Conrad M. Ferroptosis: mechanisms, biology and role in disease. Nat Rev Mol Cell Biol. 2021;22(4):266-282. PMID 33495651",
 4:"Friedmann Angeli JP, Schneider M, Proneth B, et al. Inactivation of the ferroptosis regulator Gpx4 triggers acute renal failure in mice. Nat Cell Biol. 2014;16(12):1180-1191. PMID 25402683",
 5:"Doll S, Freitas FP, Shah R, et al. FSP1 is a glutathione-independent ferroptosis suppressor. Nature. 2019;575(7784):693-698. PMID 31634899",
 6:"Bersuker K, Hendricks JM, Li Z, et al. The CoQ oxidoreductase FSP1 acts parallel to GPX4 to inhibit ferroptosis. Nature. 2019;575(7784):688-692. PMID 31634900",
 7:"Mao C, Liu X, Zhang Y, et al. DHODH-mediated ferroptosis defence is a targetable vulnerability in cancer. Nature. 2021;593(7860):586-590. PMID 33981038",
 8:"Ursini F, Maiorino M. Lipid peroxidation and ferroptosis: the role of GSH and GPx4. Free Radic Biol Med. 2020;152:175-185. PMID 32165281",
 9:"Tang D, Kang R, Coyne CB, Zeh HJ, Lotze MT. PAMPs and DAMPs: signal 0s that spur autophagy and immunity. Immunol Rev. 2012;249(1):158-175. PMID 22889221",
 10:"Wen Q, Liu J, Kang R, Zhou B, Tang D. The release and activity of HMGB1 in ferroptosis. Biochem Biophys Res Commun. 2019;510(2):278-283. PMID 30686534",
 11:"Chen X, Kang R, Kroemer G, Tang D. Ferroptosis in infection, inflammation, and immunity. J Exp Med. 2021;218(6):e20210518. PMID 33978684",
 12:"Sun J, Xie T, Chen Y, et al. Interplay between natural killer cells and ferroptosis: novel insights in tumor immunity and therapeutic potential. Cell Commun Signal. 2026;24(1). PMID 41620780",
 13:"Asherman JG. Amenorrhoea traumatica (atretica). J Obstet Gynaecol Br Emp. 1948;55(1):23-30. PMID 18902559",
 14:"March CM. Asherman's syndrome. Semin Reprod Med. 2011;29(2):83-94. PMID 21437822",
 16:"Yu D, Wong YM, Cheong Y, Xia E, Li TC. Asherman syndrome — one century later. Fertil Steril. 2008;89(4):759-779. PMID 18406834",
 17:"Zhu Q, Chen F, Luo F, et al. Ferroptosis contributes to endometrial fibrosis in intrauterine adhesions. Free Radic Biol Med. 2023;205:151-162. PMID 37302615",
 18:"Santamaria X, Izquierdo A, Gomez-Gil A, Simon C. Decoding the endometrial niche of Asherman's syndrome at single-cell resolution. Nat Commun. 2023;14:5890. PMID 37735465",
 20:"Gene Expression Omnibus accession GSE311899: paired pre/post-therapy endometrial biopsies from a phase 1/2 trial of autologous CD133+ bone-marrow-derived stem cell therapy for Asherman syndrome (EudraCT 2016-003975-23; archived 2025/12)",
 21:"Cousins FL, Gargett CE, Nguyen HPH. Endometrial stem/progenitor cells and their role in endometrial repair and regeneration. Front Reprod Health. 2021;3:811537. PMID 36304009",
 22:"Xiong D, Wang H, Chen L, et al. Platelet-rich plasma as an adjuvant therapy for intrauterine adhesions: a narrative review of mechanisms, clinical efficacy and combination strategies. Ann Med. 2026;58(1):2731715. PMID 42751776",
 23:"Xu X, Zhang Y, Liu J, et al. Endometrial stromal Menin supports endometrial receptivity by maintaining homeostasis of WNT signaling pathway through H3K4me3 during the window of implantation. Commun Biol. 2025;8:995. PMID 40610728",
 25:"Santamaria X, Cabanillas S, Cervello I, et al. Autologous cell therapy with CD133+ bone marrow-derived stem cells for refractory Asherman's syndrome and endometrial atrophy: a pilot cohort study. Hum Reprod. 2016;31(5):1087-1096. PMID 27005892",
 26:"Xiong R, Liu Y, Zhang S, et al. Regulated cell death at the maternal-fetal interface in preeclampsia: apoptosis, necroptosis, pyroptosis, ferroptosis, autophagic cell death, and cuproptosis — a narrative review. Int J Gen Med. 2026;19:637273. PMID 42733653",
 27:"Bashah A, Almohawes Z, Alrouji M, et al. Iron homeostasis and macrophage polarization in oral squamous cell carcinoma: mechanisms and therapeutic perspectives. Front Immunol. 2026;17:1863727. PMID 42591785",
 29:"Foroutan M, Bhuva DD, Cursons J, Davis MJ. Single sample scoring of molecular phenotypes. BMC Bioinformatics. 2018;19:404. PMID 30400809",
 30:"Zhou N, Bao J. FerrDb V2: update of the manually curated database of ferroptosis regulators and ferroptosis-disease associations. Nucleic Acids Res. 2023;51(D1):D571-D582. PMID 36305834",
}

text = open(MD, encoding="utf-8").read()

# 1) collect all bracket groups containing '#', expand, record first-appearance order
order = []  # old numbers in first-appearance order
groups = re.findall(r"\[#[0-9][^\]]*\]", text)
for g in groups:
    inner = g[1:-1]
    for part in re.split(r"[,;]", inner):
        part = part.strip().replace("–", "-").replace("—", "-")
        if not part.startswith("#"):
            continue
        m = re.match(r"#(\d+)(?:-(\d+))?$", part)
        if not m:
            continue
        a, b = int(m.group(1)), int(m.group(2) or m.group(1))
        for k in range(a, b + 1):
            if k not in order:
                order.append(k)
missing = [k for k in order if k not in REFS]
assert not missing, f"tokens without citation text: {missing}"

mapping = {old: i + 1 for i, old in enumerate(order)}

# 2) replace each bracket group
def repl_group(m):
    inner = m.group(0)[1:-1]
    nums = []
    for part in re.split(r"[,;]", inner):
        part = part.strip().replace("–", "-").replace("—", "-")
        if not part.startswith("#"):
            continue
        mm = re.match(r"#(\d+)(?:-(\d+))?$", part)
        if not mm:
            continue
        a, b = int(mm.group(1)), int(mm.group(2) or mm.group(1))
        nums.extend(mapping[k] for k in range(a, b + 1))
    nums = sorted(set(nums))
    return "[" + ",".join(str(n) for n in nums) + "]"

text = re.sub(r"\[#[0-9][^\]]*\]", repl_group, text)

# 3) build references section
lines = ["## References", ""]
for old in order:
    lines.append(f"{mapping[old]}. {REFS[old]}.")
lines.append("")
text = text.replace("## References\n\n[TO BE AUTO-NUMBERED BY SCRIPT]\n", "\n".join(lines) + "\n")

open(MD, "w", encoding="utf-8").write(text)
print(f"renumbered {len(order)} refs; order(old) = {order}")
print("new numbering:", {mapping[o]: f'#{o}' for o in order})
