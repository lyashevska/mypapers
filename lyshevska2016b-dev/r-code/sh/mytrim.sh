#to trim whitespace around the figure
for i in *.png ; do convert $i -trim $i ; done
