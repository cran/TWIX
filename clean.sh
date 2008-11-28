#!/bin/bash
for k in ./.DS*; do
	rm -f $k
done

for k in ./data/.DS*; do
	rm -f $k
done

for k in ./src/.DS*; do
	rm -f $k
done

for k in ./man/.DS*; do
	rm -f $k
done

for k in ./R/.DS*; do
	rm -f $k
done

for k in ./src/*.o; do
	rm -f $k
done

for k in ./src/*.so; do
	rm -f $k
done

exit 0