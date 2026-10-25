#!/usr/bin/env bash
# Starts the Materilens server setup from Oracle Cloud Shell. Run:
#   bash <(curl -fsSL https://raw.githubusercontent.com/drxavier369-spec/materilens-xrd-library/main/deploy/oracle.sh)
# It asks for a GitHub token that can read the private materials-app repository, checks it,
# then runs materials-app/deploy/oracle/cloudshell.sh.
set -uo pipefail
API=https://api.github.com/repos/drxavier369-spec/materials-app
while true; do
  echo
  echo "Paste your GitHub token here (right-click > Paste), then press Enter:"
  read -r GH_TOKEN
  GH_TOKEN=$(printf '%s' "$GH_TOKEN" | tr -d '[:space:]')
  [ -n "$GH_TOKEN" ] || { echo "Nothing was pasted. Try again."; continue; }
  code=$(curl -s -o /dev/null -w '%{http_code}' -H "Authorization: Bearer $GH_TOKEN" https://api.github.com/user)
  if [ "$code" != 200 ]; then
    echo "GitHub does not accept this token (it may be deleted or expired). Make a new one and paste it."
    continue
  fi
  code=$(curl -s -o /dev/null -w '%{http_code}' -H "Authorization: Bearer $GH_TOKEN" "$API")
  if [ "$code" != 200 ]; then
    echo "The token works but cannot see materials-app. On GitHub, edit the token: Repository access > materials-app, Contents: Read-only."
    continue
  fi
  echo "Token OK."
  break
done
export GH_TOKEN
clear
# saved to a file first: piped into bash, the script's own ssh/oci calls could swallow the rest of it from stdin
curl -fsSL -H "Authorization: Bearer $GH_TOKEN" -H "Accept: application/vnd.github.raw" \
  "$API/contents/deploy/oracle/cloudshell.sh?ref=prototype" -o ~/materilens-setup.sh
bash ~/materilens-setup.sh
