#!/bin/sh

COMMIT_HASH=2f14aa748f200290400307b40116787a0dbf5d43
URL=https://raw.githubusercontent.com/kbarbary/bsplines

wget -nc ${URL}/${COMMIT_HASH}/src/bs.c
wget -nc ${URL}/${COMMIT_HASH}/src/bs.h
wget -nc ${URL}/${COMMIT_HASH}/bsplines.pxi
