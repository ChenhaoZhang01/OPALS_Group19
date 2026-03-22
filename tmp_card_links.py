import urllib.request
import re

u = "https://card.mcmaster.ca/download"
with urllib.request.urlopen(u, timeout=30) as r:
    h = r.read().decode("utf-8", "ignore")

print("len", len(h))
for token in ["latest/data", "readme/data", ".fasta", ".fa", "card.json", "aro_categories", "download"]:
    print(token, token in h)

for m in re.findall(r"https?://[^\"'<>\s]+", h):
    if "card.mcmaster.ca" in m:
        print(m)
