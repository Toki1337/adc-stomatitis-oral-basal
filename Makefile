ROOT := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
PYTHON = PYTHONPATH=$(ROOT)/scripts PYTHONUNBUFFERED=1 python3 -u
LOGDIR = $(ROOT)/logs

.PHONY: all part1 part3 part5 part4 numbers figures supp docx verify clean_logs

all: part1 part3 part5 part4 numbers figures supp docx verify

part1:
	@mkdir -p $(LOGDIR)
	$(PYTHON) $(ROOT)/scripts/part1_dual_extraction.py > $(LOGDIR)/part1_dual_extraction.log 2>&1
	$(PYTHON) $(ROOT)/scripts/part1_hpa_dab.py > $(LOGDIR)/part1_hpa_dab.log 2>&1
	$(PYTHON) $(ROOT)/scripts/part1_model.py > $(LOGDIR)/part1_model.log 2>&1
	$(PYTHON) $(ROOT)/scripts/part1_robust.py > $(LOGDIR)/part1_robust.log 2>&1
	$(PYTHON) $(ROOT)/scripts/part1_celltype_scan.py > $(LOGDIR)/part1_celltype_scan.log 2>&1
	$(PYTHON) $(ROOT)/scripts/part1_faers2.py > $(LOGDIR)/part1_faers2.log 2>&1

part3:
	@mkdir -p $(LOGDIR)
	$(PYTHON) $(ROOT)/scripts/part3_absolute_vs_ratio.py > $(LOGDIR)/part3_absolute_vs_ratio.log 2>&1

part5:
	@mkdir -p $(LOGDIR)
	$(PYTHON) $(ROOT)/scripts/part5_calib.py > $(LOGDIR)/part5_calib.log 2>&1

part4:
	@mkdir -p $(LOGDIR)
	$(PYTHON) $(ROOT)/scripts/surfaceome_build.py > $(LOGDIR)/surfaceome_build.log 2>&1
	$(PYTHON) $(ROOT)/scripts/part4_final.py > $(LOGDIR)/part4_final.log 2>&1
	$(PYTHON) $(ROOT)/scripts/part6_madan_and_vitalorgan.py > $(LOGDIR)/part6_madan_and_vitalorgan.log 2>&1

numbers:
	$(PYTHON) $(ROOT)/scripts/export_fig_data.py
	$(PYTHON) $(ROOT)/scripts/export_numbers.py

figures:
	Rscript $(ROOT)/scripts/make_figures.R

supp:
	$(PYTHON) $(ROOT)/scripts/build_supp_tables.py

docx:
	$(PYTHON) $(ROOT)/scripts/build_docx.py

verify:
	$(PYTHON) $(ROOT)/scripts/verify_manuscript_numbers.py
	Rscript -e "suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(readr);library(ggrepel);library(cowplot);library(scales)}); capture.output(sessionInfo(), file='$(ROOT)/R_sessionInfo.txt')"
