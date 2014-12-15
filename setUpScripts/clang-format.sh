#!/bin/bash

DIRS="src"

find $DIRS -regex ".*\.[cChH]\(pp\)?" | xargs clang-format-3.4 -i --style=Google
