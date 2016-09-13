TEXFILE = second_draft

#this the the rendered and knit R script to html for browsing intermediate figures and looking at analysis decisions

analysis_html_documentation: *.html

Data/crystal_violet.RData: Data/crystal_violet_biofilm/*
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

Data/hla_tidy.RData: Data/hla/*
	Rscript hla_tidy.R
	rm Rplots.pdf

Data/merged_raft_cfu.RData: Data/hla_tidy.RData Data/psm_module/* Data/rafts_cfu_merged/*
	Rscript rafts_mutants_datawrangle.R
	rm Rplots.pdf

rafts_mutants_datawrangle.html:
	Rscript -e "rmarkdown::render(knitr::spin('rafts_mutants_datawrangle.R', knit = FALSE), output_dir = 'analysis_html')";
	rm rafts_mutants_datawrangle.Rmd

rafts_mutants_stats_vis.html: Data/merged_raft_cfu.RData
	Rscript -e "rmarkdown::render(knitr::spin('rafts_mutants_stats_vis.R', knit = FALSE), output_dir = 'analysis_html')";
	rm rafts_mutants_stats_vis.Rmd

rafts_agr_stats_vis.html: Data/merged_raft_cfu.RData
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

Data/cfu_calibration_merged.RData: Data/calibration/*
	Rscript calibration_datamerge.R;
	rm Rplots.pdf;
	echo "$@ created, it is a large file..."

Data/cfu_calibration.RData: Data/cfu_calibration_merged.RData Data/calibration/*
	Rscript calibration_datamerge.R;
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

allfigs: mainfigs suppfigs
	touch $@;

mainfigs: Figures/fig1.pdf Figures/fig2.pdf Figures/fig3.pdf Figures/fig4.pdf
	touch $@;

suppfigs: Figures/figS1.pdf Figures/figS2.pdf Figures/figS3.pdf Figures/figS4.pdf Figures/figS5.pdf Figures/figS6.pdf
	touch $@;

$(TEXFILE).pdf: $(TEXFILE).tex
	pdflatex -interaction nonstopmode -halt-on-error -file-line-error $(TEXFILE).tex;
	bibtex $(TEXFILE).aux;
	pdflatex -interaction nonstopmode -halt-on-error -file-line-error $(TEXFILE).tex;
	pdflatex -interaction nonstopmode -halt-on-error -file-line-error $(TEXFILE).tex;
	rm -v $(TEXFILE).blg $(TEXFILE).out $(TEXFILE).aux $(TEXFILE).bbl $(TEXFILE).log