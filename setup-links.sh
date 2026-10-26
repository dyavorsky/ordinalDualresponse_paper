#!/bin/sh
# Recreate the Dropbox symlinks this project expects. Safe to re-run.
#
# Only shared, paper-related links belong here. Personal links (draft talk
# material and the like) are kept in local notes outside the repo, since this
# repository is public.
set -e
cd "$(dirname "$0")"

link() {
  target=$1; name=$2
  [ -d "$target" ] || { echo "missing: $target"; echo "  (adjust the path in setup-links.sh)"; return 1; }
  [ -L "$name" ] && rm "$name"
  [ -e "$name" ] && { echo "refusing to replace real directory: $name"; return 1; }
  mkdir -p "$(dirname "$name")"
  ln -s "$target" "$name"
  echo "linked $name -> $target"
}

link "$HOME/Dropbox/research_projects/2026_ordinal_dualresponse/lit" "lit"
