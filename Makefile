


all: rrcg

rrcg: rrcg.c
	gcc rrcg.c -o rrcg -lm

install: rrcg
	cp rrcg $(HOME)/.local/bin
