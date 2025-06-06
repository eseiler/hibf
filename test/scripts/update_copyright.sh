#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

usage="\
SYNOPSIS
    update_copyright.sh <old_year> <new_year> [<files>...]

DESCRIPTION
    Updates the copyright year of files that are formatted in a certain way. Prints the
    copyright years that it ignores.

EXAMPLES
    ./test/scripts/update_copyright.sh 2023 2024 \$(find . -type f -not -path './\.**' -and -not -path './build/**' -and -not -path './include/hibf/contrib/std/**')
        Updates all copyright entries from 2023 to 2024. Only scans non hidden directories. Does not scan build
        directory and git subtree.
"

if [ $# -eq 0 ]; then
    echo -e "$usage"
    exit 1
fi

# New update year
oldyear=$1
year=$2
shift 2 # "Consumes" 2 positionial arguments. After this, $1 will be the first file.

echo "Setting copyright dates from ${oldyear} to ${year}."

for file in "$@"; do
    perl -i -pe 's/^(.*SPDX-FileCopyrightText: [0-9]{4}-)'${oldyear}'(,? Knut Reinert.*$)/${1}'${year}'${2}/' $file
    perl -ne 'print "'$file':$.: $_" if (/^.*SPDX-FileCopyrightText.*'${oldyear}'.*$/);' $file
done
