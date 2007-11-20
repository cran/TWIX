#!/bin/bash
for f in ./src/*.o; do
	rm -f $f
done

for f in ./src/*.so; do
	rm -f $f
done

exit 0