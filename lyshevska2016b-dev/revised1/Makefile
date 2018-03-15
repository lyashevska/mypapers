# Makefile for dvips
# now you want to do:
# make PAPER=paper.tex

# $^  eh tudo que esta a direita do : ou seja as dependencias
# a linha $(addsuffix .ps,$(basename $(PAPER)))
# sao comandos do make, mas que voce tem tb na sua shell
# man basename
# eu primeiro extraio o basename do conteudo da variavel $PAPER
# e depois eu adiciono o sufixo ".ps"

all:  $(addsuffix .ps,$(basename $(PAPER))) # text.ps
	ps2pdf $^

%.ps:  $(addsuffix .dvi,$(basename $(PAPER))) # text.dvi 
	dvips $^
    
%.dvi: $(PAPER)
	latex $^ && bibtex $(basename $^) && latex $^ && latex $^ 


#alternatively, if we do not want to use a specific file name
# Makefile for dvips
# build my pdf based on tex
# to build it:
# make 

#all: paper.ps
#	ps2pdf paper.ps

#paper.ps: paper.dvi 
#	dvips paper.dvi
    
#paper.dvi: paper.tex 
#	latex paper.tex && bibtex paper && latex paper.tex && latex paper.tex 
