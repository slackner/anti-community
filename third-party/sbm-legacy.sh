#!/bin/bash
set -e -x

DIR="$(cd "$(dirname "$0")" && pwd)"
"$DIR/../src/gml" "$1" | "$DIR/estimate"
