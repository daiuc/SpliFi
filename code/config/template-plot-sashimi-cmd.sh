#!/bin/bash
pyGenomeTracks --tracks {{ iniFile }} \
  --region {{ plotRange }} \
  --title {{ plotTitle}} \
  --height 10 \
  --width 20 \
  --trackLabelFraction 0.1 \
  --fontSize 7 \
  --dpi 150 \
  -out {{ plotFile }}


