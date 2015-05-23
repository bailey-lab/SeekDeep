#!/usr/bin/env python

import os, sys, base64

try:
    from bs4 import BeautifulSoup
except:
    print "please pip install beautifulsoup4"
    sys.exit(-1)

if len(sys.argv) != 2:
    print "please specify one HTML file"
    sys.exit(-1)

fn = sys.argv[1]

with open(fn) as f:
    soup = BeautifulSoup(f.read())

mimeForPng = "data:image/png;base64,"
num = 1
for a in soup.find_all('img'):
    data = a.get("src")
    if not data.startswith(mimeForPng):
        print "unknown mime type; skipping...."
        continue
    data = data.replace(mimeForPng, '')
    imgData = base64.b64decode(data)
    imgFn = fn + str(num) + ".png"
    num += 1
    with open(imgFn, 'wb') as f:
        f.write(imgData)
    print "wrote", imgFn
