bin/coloring: bin/cmain.o bin/probdata_coloring.o bin/reader_coloring.o bin/pricer_coloring.o bin/vardata_coloring.o
	gcc -o bin/coloring bin/cmain.o bin/probdata_coloring.o bin/reader_coloring.o bin/pricer_coloring.o bin/vardata_coloring.o -lscip

bin/cmain.o: src/cmain.c
	gcc -c -o bin/cmain.o src/cmain.c

bin/probdata_coloring.o: src/probdata_coloring.c src/probdata_coloring.h
	gcc -c -o bin/probdata_coloring.o src/probdata_coloring.c

bin/vardata_coloring.o: src/vardata_coloring.c src/vardata_coloring.h
	gcc -c -o bin/vardata_coloring.o src/vardata_coloring.c

bin/reader_coloring.o: src/reader_coloring.c src/reader_coloring.h
	gcc -c -o bin/reader_coloring.o src/reader_coloring.c

bin/pricer_coloring.o: src/pricer_coloring.c src/pricer_coloring.h
	gcc -c -o bin/pricer_coloring.o src/pricer_coloring.c

.PHONY: clean

clean:
	rm -f bin/*.o bin/coloring

