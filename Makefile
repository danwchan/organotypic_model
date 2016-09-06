#this the the rendered and knit R script to html for browsing intermediate figures and looking at analysis decisions
analysis_html_documentation: crystal_violet_stats_vis.html

crystal_violet_stats_vis.html: Data/crystal_violet.Rdata
	Rscript -e "rmarkdown::render(knitr::spin('crystal_violet_stats_vis.R', knit = FALSE), output_dir = 'analysis_html')";
	rm "crystal_violet_stats_vis.Rmd"

Data/crystal_violet.RData: 
	Rscript "crystal_violet_datawrangle.R";

rafts_mutants_stats_vis.html: Data/hla_tidy.RData
	Rscript -e "rmarkdown::render(knitr::spin('rafts_mutants_stats_vis.R', knit = FALSE), output_dir = 'analysis_html')";
	rm "crystal_violet_stats_vis.Rmd"

Data/hla_tidy.RData: 
	Rscript "hla_tidy.R"
	rm "Rplots.pdf"

#these were made later with moe understanding
22251_effsize.pdf:
	Rscript -e "rmarkdown::render(knitr::spin('./Figures/160606_clinical_strain_panel/160608_22251_effsize.R', knit = FALSE))"

22251_agrA_comp.pdf:
	Rscript -e "rmarkdown::render(knitr::spin('./Figures/160606_clinical_strain_panel/160608_22251_agrA_comp.R', knit = FALSE))"

#export svg inkscape files as pdfs for latex import
# why am I exporting my figure as .tiff files then deleting them?:
# because inkscape imports vector based things as a collection of objects that do not update dynamically in the inkscape file
# since I want dynamiclly updating figures and I want high quality I must export graphs as
allfigs: mainfigs suppfigs

mainfigs: fig1_characterize_model.pdf fig2_SEM.pdf fig3_growth.pdf fig4_biofilm.pdf

suppfigs: figS1_experimental_setup.pdf figS2_sonicated.pdf figS3_22251.pdf figS4_newman.pdf figS5_PBS.pdf

test.tiff:
	Rscript "crystal_violet_stats_vis.R";
	echo "$@ created, it is a large file..."

test_fig1.pdf: test.tiff
	inkscape "Figures/test_inkscape.svg" -A "Figures/for_publication/test_fig1.pdf";
	rm "Figures/test.tiff"

#the IF images are a blackspot on my ability to keep good records of my process... the original images are unclear and I am (at the moment) unwilling to dive through and find them, and also rederive the exact manipulations to the image
fig1_characterize_model.pdf: lucifer_yellow.pdf
	inkscape "./Figures/Fig1_forpub_model.svg" -A "./Figures/fig1_characterize_model.pdf"

fig2_SEM.pdf: 
	inkscape "./Figures/Fig2_forpub_SEM.svg" -A "./Figures/fig2_SEM.pdf"

fig3_growth.pdf:
	inkscape "./Figures/Fig3_forpub_mutant_growth.svg" -A "./Figures/fig3_mutant_growth.pdf"

fig4_biofilm.pdf:
	inkscape "./Figures/Fig4_forpub_biofilm.svg" -A "./Figures/fig4_biofilm.pdf"

figS1_experimental_setup.pdf:
	inkscape "./Figures/figS1_experimental_setup.svg" -A "./Figures/figS1_experimental_setup.pdf"

figS2_effsize.pdf:
	inkscape "./Figures/figS2_effsize.svg" -A "./Figures/figS2_effsize.pdf"

figS2_sonicated.pdf:
	inkscape "./Figures/figS3_sonicated.svg" -A "./Figures/figS2_sonicated.pdf"

figS3_22251.pdf:
	inkscape "./Figures/Fig4_forpub_22251.svg" -A "./Figures/figS3_22251.pdf"

figS4_newman.pdf:
	inkscape "./Figures/figS4_newman.svg" -A "./Figures/figS4_newman.pdf"

figS5_PBS.pdf:
	inkscape "./Figures/figS5_PBS.svg" -A "./Figures/figS5_PBS.pdf"

#having some trouble here... when the lesson is learned...
figures_proofing.pdf: allfigs
	cd "./Figures"
	pdflatex "./Figures/figures_proofing.tex"

#make references for pasting into plos latex template since they don't support .bib file submission
make_ref_plos:
	latex first_draft.tex
	bibtex first_draft.aux