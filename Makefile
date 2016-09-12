TEXFILE = second_draft

#this the the rendered and knit R script to html for browsing intermediate figures and looking at analysis decisions
analysis_html_documentation: *.html

Data/crystal_violet.RData: 
	Rscript crystal_violet_datawrangle.R;
	rm Rplots.pdf

crystal_violet_datawrangle.html:
	Rscript -e "rmarkdown::render(knitr::spin('crystal_violet_datawrangle.R', knit = FALSE), output_dir = 'analysis_html')";
	rm crystal_violet_datawrangle.Rmd

crystal_violet_stats_vis.html: Data/crystal_violet.Rdata
	Rscript -e "rmarkdown::render(knitr::spin('crystal_violet_stats_vis.R', knit = FALSE), output_dir = 'analysis_html')";
	rm crystal_violet_stats_vis.Rmd

crystal_violet_figures:
	Rscript crystal_violet_stats_vis.R
	rm Rplots.pdf
	echo "$@ created, it is a large file..."

crystal_violet_agr_figures:
	Rscript crystal_violet_agr_stats_vis.R
	rm Rplots.pdf
	echo "$@ created, it is a large file..."

Data/hla_tidy.RData: 
	Rscript hla_tidy.R
	rm Rplots.pdf

Data/merged_raft_cfu.RData:
	Rscript rafts_mutants_datawrangle.R
	rm Rplots.pdf

rafts_mutants_datawrangle.html:
	Rscript -e "rmarkdown::render(knitr::spin('rafts_mutants_datawrangle.R', knit = FALSE), output_dir = 'analysis_html')";
	rm rafts_mutants_datawrangle.Rmd

rafts_mutants_stats_vis.html: Data/hla_tidy.RData Data/merged_raft_cfu.RData
	Rscript -e "rmarkdown::render(knitr::spin('rafts_mutants_stats_vis.R', knit = FALSE), output_dir = 'analysis_html')";
	rm rafts_mutants_stats_vis.Rmd

rafts_agr_stats_vis.html:
	Rscript -e "rmarkdown::render(knitr::spin('rafts_agr_stats_vis.R', knit = FALSE), output_dir = 'analysis_html')";
	rm rafts_agr_stats_vis.Rmd

rafts_mutants_figures:
	Rscript rafts_mutants_stats_vis.R;
	rm Rplots.pdf;
	echo "$@ created, it is a large file..."

rafts_agr_figures:
	Rscript rafts_agr_stats_vis.R;
	rm Rplots.pdf;
	echo "$@ created, it is a large file..."


cfu_calibration_datawrangle.html:
	Rscript -e "rmarkdown::render(knitr::spin('cfu_calibration_datawrangle.R', knit = FALSE), output_dir = 'analysis_html')";
	rm cfu_calibration_datawrangle.Rmd

cfu_calibration_figures:
	Rscript cfu_calibration_stats_vis.R
	rm Rplots.pdf
	echo "$@ created, it is a large file..."

Figures/fig1.pdf: cfu_calibration_figures
	inkscape Figures/Figure1.svg -A $@;
	rm Figures/*.tiff

Figures/fig2.pdf:
	inkscape Figures/Figure2.svg -A $@;
	rm Figures/*.tiff

Figures/fig3.pdf: rafts_agr_figures
	inkscape Figures/Figure3.svg -A $@;
	rm Figures/*.tiff

Figures/fig4.pdf: crystal_violet_agr_figures
	inkscape Figures/Figure4.svg -A $@;
	rm Figures/*.tiff

Figures/figS1.pdf:
	inkscape Figures/Supplement_figure1.svg -A $@

Figures/figS2.pdf:
	inkscape Figures/Supplement_figure2.svg -A $@

Figures/figS3.pdf: rafts_mutants_figures
	inkscape Figures/Supplement_figure3.svg -A $@;
	rm Figures/*.tiff

Figures/figS4.pdf:
	inkscape Figures/Supplement_figure4.svg -A $@

Figures/figS5.pdf:
	inkscape Figures/Supplement_figure5.svg -A $@

Figures/figS6.pdf: crystal_violet_figures
	inkscape Figures/Supplement_figure6.svg -A $@;
	rm Figures/*.tiff

Figures/%.pdfcrop: Figures/%.pdf
	pdfcrop $< $<;
	touch $@;
	echo "$@ has created and cropped the figure, .pdfcrop token created"

$(TEXFILE).pdf: $(TEXFILE).tex
	pdflatex -interaction nonstopmode -halt-on-error -file-line-error $(TEXFILE).tex;
	bibtex $(TEXFILE).aux;
	pdflatex -interaction nonstopmode -halt-on-error -file-line-error $(TEXFILE).tex;
	pdflatex -interaction nonstopmode -halt-on-error -file-line-error $(TEXFILE).tex;
	rm -v $(TEXFILE).blg $(TEXFILE).out $(TEXFILE).aux $(TEXFILE).bbl $(TEXFILE).log


#export svg inkscape files as pdfs for latex import
# why am I exporting my figure as .tiff files then deleting them?:
# because inkscape imports vector based things as a collection of objects that do not update dynamically in the inkscape file
# since I want dynamiclly updating figures and I want high quality I must export graphs as
allfigs: mainfigs suppfigs

mainfigs: fig1_characterize_model.pdf fig2_SEM.pdf fig3_growth.pdf fig4_biofilm.pdf

suppfigs: figS1_experimental_setup.pdf figS2_sonicated.pdf figS3_22251.pdf figS4_newman.pdf figS5_PBS.pdf
