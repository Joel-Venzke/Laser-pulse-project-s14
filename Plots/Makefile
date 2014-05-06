FIGURES	=	overlap_440_1 overlap_444_1 overlap_440_4 overlap_444_4\
	        overlap_410_1 overlap_430_1 overlap_450_1 \
	         overlap_350_4 overlap_350_1 overlap_350_025 \
				approx_ionization \
				overlap_375_4 overlap_375_1 overlap_375_025
# 				approx_ionization_0375 \
# 				overlap_3935_4 overlap_3799_1 overlap_37625_025 \

all: figs
	rm *.dvi *.log *.aux
	cp *.ps Plots/
	@echo ""
	@echo "Compiling figures complete"

figs: 
	for fig in $(FIGURES); do \
		gnuplot -e "datadir='Data'" $$fig.gp; \
		latex c_$$fig.tex; \
		dvips c_$$fig.dvi; \
	done

clean:
	rm -f *.aux *.bbl *.blg *.d *.dvi *.out *.ps *.log *.table
	for fig in $(FIGURES); do \
		rm -f $$fig.tex; \
	done