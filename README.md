## About

I am a Research Software Engineer and Data Scientist working at the intersection of machine learning, statistical modelling, and natural language processing.

My work spans marine and ecological research, medicine, and the social sciences. I build reproducible data pipelines, predictive models, and research software that combine domain expertise with robust, interpretable methods.

My publications range from methodological research—such as approaches to class imbalance and spatial zero-inflation—to applied studies of long-term species trends, gelatinous zooplankton, microplastics, and medical data.

I value practical, well-documented code and reproducible workflows that make complex analyses transparent, reliable, and reusable across projects.

## How this repository works

This repository is designed to keep publications (BibTeX entries and PDFs) together and to regenerate the publications list automatically in CI.

Key points:

- Place BibTeX files in the `bib/` directory and PDFs in the `pdf/` directory.
- Use the same key/name for a publication's BibTeX and PDF files. For example, for a publication with key `lyashevska2016a` provide `bib/lyashevska2016a.bib` and (optionally) `pdf/lyashevska2016a.pdf`.
- The generator script `bib_to_md.py` reads the files in `bib/`, normalizes author names, adds available metadata (volume/number/pages), and writes a generated `publications.md` (used as an intermediate artifact).

What CI does:

- The GitHub Actions workflow (`.github/workflows/update-readme.yml`) checks out the repository and installs Python dependencies (from `requirements.txt` when present).
- Before committing the updated README, CI runs `scripts/verify_readme_bibs.py` which ensures that every BibTeX key found under `bib/` is referenced in the publications block; the workflow fails if keys are missing.
- The workflow commits only the updated `README.md` (the generated `publications.md` is an intermediate artifact and does not need to be tracked).


### Testing locally

Before pushing changes, run `./test.sh` to verify everything works:

```bash
./test.sh
```

This script mirrors the CI workflow and will:
1. Check Python environment
2. Install dependencies
3. Generate publications.md
4. Update README
5. Verify all bib files are referenced

If it passes locally, CI will pass on GitHub.

**First-time setup:**
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

**Tips:**
- Keep file keys consistent: matching names between `bib/{key}.bib` and `pdf/{key}.pdf`
- See `DEVNOTES.md` for additional developer notes

<!-- PUBLICATIONS START -->
## 📚 Publications
- **Klimesova, Bela; O'Dwyer, Katie; D'Arcy, Jack; Talbot, Anita; Lyashevska, Olga; Rodger, Hamish; McManus, Catherine; Ruane, Neil M.** (2026). *Assessing Sea lice infection levels in Irish Atlantic salmon farms: Metrics and treatment trigger levels.* _Aquaculture._ 614: 743568. [PDF](pdf/klimesova2026a.pdf) | [BibTeX](bib/klimesova2026a.bib)
- **Klimesova, B.; Ruane, N. M.; Domingo-Bretón, R.; Moroni, F.; Naya-Català, F.; Pérez-Sánchez, J.; O'Dwyer, K.; Lyashevska, O.; Rodger, H.; Talbot, A.** (2026). *Sea Lice (Lepeophtheirus salmonis) Harbour Putative Fish Pathogens: Insights From Illumina and Nanopore Sequencing.* _Journal of Fish Diseases._ 49(9): e70182. [PDF](pdf/klimesova2026b.pdf) | [BibTeX](bib/klimesova2026b.bib)
- **Farmani, Vahid; Kniep, Helge; Maros, Mate E.; Lyashevska, Olga; Malone, Fiona; Fiehler, Jens; Morris, Liam** (2025). *Estimating Individualized Effectiveness of Receiving Successful Recanalization for Ischemic Stroke Cases Using Machine Learning Techniques.* _Journal of Stroke and Cerebrovascular Diseases._ [PDF](pdf/farmani2025.pdf) | [BibTeX](bib/farmani2025.bib)
- **Klimesova B., Lyashevska O., Ruane N., D’Arcy, J., Talbot A., Rodger H., O’Dwyer K.** (2025). *Effects of temporal, geographical and environmental factors on salmon lice (Lepeophtheirus
salmonis) levels of Atlantic salmon (Salmo salar) in Ireland.* _Scientific Reports._ 15: 34614. [PDF](pdf/klimesova2025.pdf) | [BibTeX](bib/klimesova2025.bib)
- **Long, A.P.; Bastian, T.; Haberlin, D.; Stokes, D.; Lyashevska, O.; Brophy, D.; Lawton, C.; Doyle, T.K.** (2024). *Regular widespread aggregations of the oceanic jellyfish Pelagia noctiluca in the northeast Atlantic over 11 years.* _Estuarine, Coastal and Shelf Science._ 303: 108805. [PDF](pdf/long2024.pdf) | [BibTeX](bib/long2024.bib)
- **Williams, Rosie S.; Brownlow, Andrew; Baillie, Andrew; Barber, Jonathan L.; Barnett, James; Davison, Nicholas J.; Deaville, Robert; Doeschate, Mariel ten; Penrose, Rod; Perkins, Matthew; Williams, Ruth; Jepson, Paul D.; Lyashevska, Olga; Murphy, Sinéad** (2023). *Evaluation of a marine mammal status and trends contaminants indicator for European waters.* _Science of The Total Environment._ 866: 161301. [PDF](pdf/williams2023.pdf) | [BibTeX](bib/williams2023.bib)
- **Lyashevska, Olga; Brophy, Deirdre; Wing, Steve; Johns, David G.; Haberlin, Damien; Doyle, Thomas K.** (2022). *Evidence of a range expansion in sunfish from 47 years of coastal sightings.* _Marine Biology._ 169(2): 20. [PDF](pdf/lyashevska2022-sunfish.pdf) | [BibTeX](bib/lyashevska2022-sunfish.bib)
- **Long, A P; Haberlin, D; Lyashevska, O; Brophy, D; O’ Hea, Brendan; O’Donnell, C; Scarrott, R G; Lawton, C; Doyle, T K** (2021). *Interannual variability of gelatinous mesozooplankton in a temperate shelf sea: greater abundance coincides with cooler sea surface temperatures.* _ICES Journal of Marine Science._ [PDF](pdf/long_etal_2021.pdf)
- **Deschepper, Inge; Lyons, Kieran; Lyashevska, Olga; Brophy, Deirdre** (2020). *Biophysical models reveal the role of tides, wind, and larval behaviour in early transport and retention of Atlantic herring (Clupea harengus) in the Celtic Sea.* _Canadian Journal of Fisheries and Aquatic Sciences._ 77(1): 90--107. [PDF](pdf/deschepper2019.pdf) | [BibTeX](bib/deschepper2019.bib)
- **Frias, João P.G.L.; Lyashevska, Olga; Joyce, Haleigh; Pagter, Elena; Nash, Róisín** (2020). *Floating microplastics in a coastal embayment: A multifaceted issue.* _Marine Pollution Bulletin._ 158: 111361. [PDF](pdf/frias2020.pdf) | [BibTeX](bib/frias2020.bib)
- **Lyashevska, Olga; Malone, Fiona; MacCarthy, Eugene; Fiehler, Jens; Buhk, Jan-Hendrik; Morris, Liam** (2020). *Class imbalance in gradient boosting classification algorithms: Application to experimental stroke data.* _Statistical Methods in Medical Research._ 0(0): 0962280220980484. [PDF](pdf/lyashevska2020-class-imbalance.pdf) | [BibTeX](bib/lyashevska2020-class-imbalance.bib)
- **Lyashevska, Olga; Harma, Clementine; Minto, Cóilín; Clarke, Maurice; Brophy, Deirdre** (2020). *Long-term trends in herring growth primarily linked to temperature by gradient boosting regression trees.* _Ecological Informatics._ 60: 101154. [PDF](pdf/lyashevska2020-herring.pdf) | [BibTeX](bib/lyashevska2020-herring.bib)
- **Zhang, Zhongheng; Zhao, Yiming; Canes, Aran; Steinberg, Dan; Lyashevska, Olga** (2019). *Predictive analytics with gradient boosting in clinical medicine.* _Annals of Translational Medicine._ 7(7): 152. [PDF](pdf/zhang2019.pdf) | [BibTeX](bib/zhang2019.bib)
- **Acampora, Heidi; White, Philip; Lyashevska, Olga; O'Connor, Ian** (2018). *Contrasting congener profiles for persistent organic pollutants and PAH monitoring in European storm petrels (Hydrobates pelagicus) breeding in Ireland: a preen oil versus feathers approach.* _Environmental Science and Pollution Research._ 25(17): 16933-16944. [PDF](pdf/acampora2018.pdf) | [BibTeX](bib/acampora2018.bib)
- **Kanhai, La Daana K.; Gårdfeldt, Katarina; Lyashevska, Olga; Hassellöv, Martin; Thompson, Richard C.; O'Connor, Ian** (2018). *Microplastics in sub-surface waters of the Arctic Central Basin.* _Marine Pollution Bulletin ._ 130: 8 - 18. [PDF](pdf/kanhai2018.pdf) | [BibTeX](bib/kanhai2018.bib)
- **Acampora, Heidi; White, Philip; Lyashevska, Olga; O'Connor, Ian** (2017). *Presence of persistent organic pollutants in a breeding common tern (Sterna hirundo) population in Ireland.* _Environmental Science and Pollution Research._ 1--11. [PDF](pdf/acampora2017.pdf) | [BibTeX](bib/acampora2017.bib)
- **Lyashevska, Olga; Brus, Dick J.; Meer, Jaap van der** (2016). *Grid-spacing and the quality of abundance maps for species that show spatial autocorrelation and zero-inflation.* _Spatial Statistics ._ 18, Part B: 386 - 395. [PDF](pdf/lyashevska2016b.pdf) | [BibTeX](bib/lyashevska2016b.bib)
- **Acampora, Heidi; Lyashevska, Olga; Franeker, Jan Andries Van; O'Connor, Ian** (2016). *The use of beached bird surveys for marine plastic litter monitoring in Ireland.* _Marine Environmental Research ._ 120: 122 - 129. [PDF](pdf/acampora2016.pdf) | [BibTeX](bib/acampora2016.bib)
- **Lyashevska, Olga; Brus, Dick J.; van der Meer, Jaap** (2016). *Mapping species abundance by a spatial zero-inflated Poisson model: a case study in the Wadden Sea, the Netherlands.* _Ecology and Evolution._ [PDF](pdf/lyashevska2016a.pdf) | [BibTeX](bib/lyashevska2016a.bib)
- **Kanhai, La Daana K.; Officer, Rick; Lyashevska, Olga; Thompson, Richard C.; O'Connor, Ian** (2016). *Microplastic abundance, distribution and composition along a latitudinal gradient in the Atlantic Ocean.* _Marine Pollution Bulletin ._ -. [PDF](pdf/kanhai2017.pdf) | [BibTeX](bib/kanhai2017.bib)
- **Lyashevska, Olga; Farnsworth, Keith D.** (2012). *How many dimensions of biodiversity do we need?.* _Ecological Indicators ._ 18(0): 485 - 492. [PDF](pdf/lyashevska2012.pdf) | [BibTeX](bib/lyashevska2012.bib)
- **Farnsworth, Keith D.; Lyashevska, Olga; Fung, Tak** (2012). *Functional complexity: The source of value in biodiversity.* _Ecological Complexity ._ 11(0): 46 - 52. [PDF](pdf/farnsworth2012.pdf) | [BibTeX](bib/farnsworth2012.bib)<!-- PUBLICATIONS END -->
