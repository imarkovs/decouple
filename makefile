# makefile: SLRA makefile
all: 
	./jemdoc.py index.jemdoc
	./jemdoc.py publications.jemdoc
	./jemdoc.py software.jemdoc
	./jemdoc.py presentations.jemdoc
pull:
	git pull origin master	
push:
	git push origin master	

