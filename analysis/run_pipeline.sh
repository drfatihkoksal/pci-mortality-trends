#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
source_workbook="${project_root}/KoronerArterMudahale_birlesik_mortalite.xlsx"

if [[ ! -f "${source_workbook}" ]]; then
  echo "Source workbook not found: ${source_workbook}" >&2
  exit 1
fi

cd "${project_root}"

python3 analysis/01_extract_deidentify.py
Rscript analysis/02_analyze.R
python3 analysis/03_build_manuscript.py
python3 -m unittest discover -s tests -p 'test_*.py' -v

if command -v quarto >/dev/null 2>&1; then
  (
    cd manuscript
    quarto render paper.qmd --to docx
    quarto render paper.qmd --to html
  )
else
  (
    cd manuscript
    pandoc paper.qmd --from=markdown --citeproc --resource-path=.:../outputs/figures \
      --metadata link-citations=true -o paper.docx
    pandoc paper.qmd --from=markdown --citeproc --resource-path=.:../outputs/figures \
      --metadata link-citations=true --standalone -o paper.html
  )
fi

(
  cd manuscript
  pandoc supplement.md --resource-path=.:../outputs/figures -o supplement.docx
)

(
  cd submission
  pandoc cover_letter.md -o cover_letter.docx
  pandoc title_page.md -o title_page.docx
  pandoc pre_submission_checklist.md -o pre_submission_checklist.docx
)

if command -v libreoffice >/dev/null 2>&1; then
  render_dir="$(mktemp -d)"
  trap 'rm -rf "${render_dir}"' EXIT
  libreoffice --headless --convert-to pdf --outdir "${render_dir}" \
    "${project_root}/manuscript/paper.docx" >/dev/null
  libreoffice --headless --convert-to pdf --outdir "${render_dir}" \
    "${project_root}/manuscript/supplement.docx" >/dev/null
  mv -f "${render_dir}/paper.pdf" "${project_root}/manuscript/paper.pdf"
  mv -f "${render_dir}/supplement.pdf" "${project_root}/manuscript/supplement.pdf"
fi

echo "Pipeline completed successfully."
echo "Primary manuscript: ${project_root}/manuscript/paper.docx"
echo "Pre-submission status: INTERNAL WORKING DRAFT until checklist gates are closed."
