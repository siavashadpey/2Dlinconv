run:
	find . -maxdepth 1 -type f -name \*.gnu -exec gnuplot {} \;
	find . -maxdepth 1 -type f -name \*.eps -exec ps2pdf {} \;
	rm -f *.eps
	
