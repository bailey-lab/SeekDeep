#!/usr/bin/env python

import os

path = os.path.join(os.path.dirname(__file__), "../")
path = os.path.abspath(path)

regex = '-regex ".*\.[cChH]\(pp\)?"'
exclude = '-not -path "*/external/*" -not -name "*#*"'
cmd = 'find  {p} {r} {e} -print | xargs etags '.format(p=path, e=exclude, r=regex)

print cmd

os.system(cmd)
