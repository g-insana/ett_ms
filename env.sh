#!/usr/bin/env bash
set -u
set -e

BIN_DIR="./bin"
REQUIRED_CMDS=("tfastx36" "fasta36" "ssearch36")
declare -A CMD_URLS=(
  ["tfastx36"]="https://github.com/wrpearson/fasta36/releases/download/v36.3.8i_14-Nov-2020/fasta-36.3.8i-linux64.tar.gz"
  ["fasta36"]="https://github.com/wrpearson/fasta36/releases/download/v36.3.8i_14-Nov-2020/fasta-36.3.8i-linux64.tar.gz"
  ["ssearch36"]="https://github.com/wrpearson/fasta36/releases/download/v36.3.8i_14-Nov-2020/fasta-36.3.8i-linux64.tar.gz"
)

ensure_in_path() {
  case ":$PATH:" in
    *":$BIN_DIR:"*) ;;
    *) export PATH="$BIN_DIR:$PATH" ;;
  esac
}

download_pkg() {
  local name="$1"
  local url="$2"
  local dest="$BIN_DIR/"

  echo "Downloading $name from $url ..."
  curl -fsSL "$url" -o "$dest/package.tgz"
  (cd $dest && tar xfz package.tgz --strip-components=3 '*/bin/*')
  rm -f "$dest/package.tgz"
}

check_cmd() {
  local name="$1"
  local url="$2"

  if command -v "$name" >/dev/null 2>&1; then
    echo "$name: found"
  else
    echo "$name: not found, installing..."
    download_pkg "$name" "$url"
  fi
}

main() {
  mkdir -p $BIN_DIR
  mkdir -p outdir
  echo "Adding $BIN_DIR to PATH"
  ensure_in_path

  echo "Checking requirements"
  for cmd in "${REQUIRED_CMDS[@]}"; do
    check_cmd "$cmd" "${CMD_URLS[$cmd]}"
  done

  echo "Setting up environment variables"
  export ETT_BIN=$(pwd)/analysis_code/
  export DATA_DIR=$(pwd)/data
  export BS_ETT_DIR=$(pwd)/outdir
  export DATA_DATE=2026-05-18 #replace this with the date of your analysis
}
main "$@"
