#!/usr/bin/env bash
#
# Written by Dave Tang
# Year 2023
#
set -euo pipefail

usage(){
>&2 cat << EOF
   Usage

            $ ffs aspera

   Examples

            Read from a file

            $ cat ffq.json | ffs aspera

            Read from STDIN directly

            $ ffq --ftp SRR22891572 | ffs aspera

EOF
exit 1
}

# Read all input from stdin into a variable
json_input=$(cat)

# Check if the input is effectively empty
if [[ -z "$json_input" || "$json_input" == "[]" ]]; then
    echo "Error: JSON input is empty or '[]'." >&2
    exit 1
fi

echo "$json_input" \
   | grep '"url"' \
   | sed 's/ftp:\/\/ftp.sra.ebi.ac.uk\//era-fasp@fasp.sra.ebi.ac.uk:/' \
   | perl -lane '
      $l = $F[1];
      $l =~ s/"//g;
      ($f = $l) =~ s/.*\/(.*)$/$1/;
      print "if [[ ! -f $f ]]; then";
      print "   ascp -QT -l 300m -P33001 -i asperaweb_id_dsa.openssh $l .";
      print "fi";
   '

exit 0
