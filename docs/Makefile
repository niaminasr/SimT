

all:
	make clean
	make doc

pdf:
	cd latex && make pdf

doc:
	doxygen Doxyfile
	make pdf
	cp latex/refman.pdf ./
#ln -s html/index.html index.html
#rm -r html/ latex/ 

clean:
	rm -rf html/ latex/ 
	rm -rf *.pdf *.html
