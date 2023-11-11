#!/bin/zsh

#quarto convert "$1"
# Print the current working directory
echo "Current working directory: $(pwd)"

base=$(basename "$1" .qmd)
dir=$(dirname "$1")

# Print the message with the name of the file provided by the user
echo "Rendering file: $base.qmd"

mkdir -p _working
#cp "$1" "_working/$base.qmd"

#quarto render "_working/$base.qmd" --profile md --output "_${base}.md"
quarto render $1 --profile md 

mv _working/$dir/$base.md _working/$base.md

# if a folder that starts with the same characters as "$base" exists inside _working/$dir/, get its full name and
# move it to $_working 
# check if any folders exist in _working/$dir/ that start with the same characters as "$base"

if [[ -n $(find _working/$dir -maxdepth 1 -type d -name "$base*") ]]; then
  echo "Moving folder"
  mv _working/$dir/$base* _working/
fi

rm -r _working/$dir
rm -r _working/index.html

# Replace instances of "``` r" with "```{r}"
gsed -i 's/^``` r/```{r}/g' "_working/${base}.md"

# convert _working/$base.md to _working/$base.qmd, just by changing the extension
mv "_working/${base}.md" "_working/${base}.qmd"


# Here is the above command in a single line without inline comments
# replace code chunks with embed short codes by setting a range, deleting all but two lines, joining the lines, substituting, and deleting anything that does not match the shortcode pattern
# cat "${1%.*}.qmd"  | gsed '/^```{/,/^```$/ { /^```{.*}$\|#| tags:/!d; N;s/\n//g; s/^```{\(.*\)}#| tags: \[\(.*\)\]/{{< embed _\1\.ipynb#\2 >}}/g; /{{< embed .*\.ipynb#.* >}}/!d }' | sed "s/{{< embed \(.*\)\.ipynb#\(.*\) >}}/{{< embed ${1%.*}\1\.ipynb#\2 echo=true >}}/g" > "${1%.*}_embed.qmd"

# Ranges example: This deletes all code chunks
# cat mulang.qmd | sed '/^```{/,/^```$/d'\n

# This inserts a tabset in front of every level-four (####) JavaScript header in the embed qmd
# cat "${1%.*}"_embed.qmd | gsed '/#### JavaScript/i ::: {.panel-tabset group="language"}\n'> "${1%.*}"_tabset.qmd

# # recreate html output file with embedded notebooks
# quarto render "${1%.*}_tabset.qmd" --cache-refresh --profile knitr --metadata engine:knitr