import urllib.request
import re

u = "https://card.mcmaster.ca/latest/data"
with urllib.request.urlopen(u, timeout=30) as r:
    final = r.geturl()
    body = r.read()

print("final", final)
print("bytes", len(body))
text = body.decode("utf-8", "ignore")
print(text[:1000])
print("--- urls ---")
for m in sorted(set(re.findall(r"https?://[^\"'<>\s]+", text))):
    print(m)
print("--- hrefs ---")
for m in sorted(set(re.findall(r"href=\"([^\"]+)\"", text))):
    print(m)
