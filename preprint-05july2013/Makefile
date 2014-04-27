TARGET = emmons_drake_caps
FIGURES = c_fourier_350_2_36_2 c_fourier_350_10_20_10 c_fourier_350_20_0_20 \
		  c_pulse_350_2_36_2 c_pulse_350_10_20_10 c_pulse_350_20_0_20 c_spectrum_350_2_36_2_highin \
		  c_spectrum_350_2_36_2_lowin c_spectrum_350_3_34_3_highin c_spectrum_350_4_32_4_highin \
		  c_spectrum_375_2_36_2_highin c_spectrum_375_2_36_2_lowin
BADFIGURES = c-detune c-ramp c-shape-intens

all: figs badfigs
	pdflatex -shell-escape $(TARGET).tex
	# bibtex $(TARGET).aux
	pdflatex -shell-escape $(TARGET).tex
	sh wordcount.sh $(TARGET).tex

figs:
	for fig in $(FIGURES); do \
		ps2pdf $$fig.ps; \
	done

badfigs:
	for fig in $(BADFIGURES); do \
		ps2pdf $$fig.eps; \
		pdfcrop $$fig.pdf $$fig.pdf; \
	done

clean:
	rm -rf *.aux *.bbl *.blg *.log *.pdf