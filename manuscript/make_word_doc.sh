for f in figure_*.pdf; do
  pdftoppm -png -singlefile -r 220 "$f" "${f%.pdf}"
done

pandoc main.tex --from latex --to docx --citeproc --bibliography references.bib --resource-path . --output main.docx

rm *.png
