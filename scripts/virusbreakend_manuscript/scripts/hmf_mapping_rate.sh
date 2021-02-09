#!/bin/bash
grep "mapped (" ../protecteddata/hmf/*.flagstat | cut -f 1 -d ' ' | tr ':/.' '\t' | cut -f 6,7,9 > flagstat.tsv
