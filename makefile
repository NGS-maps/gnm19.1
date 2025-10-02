gnm  : gnm.o readargs.o plain_search.o sequence_handle.o thread_search.o event_loop.o
	gcc -O3 -lpthread -o gnm gnm.o plain_search.o readargs.o sequence_handle.o thread_search.o event_loop.o -lcurses -lm

gnm.o : gnm.c
	gcc -O3 -c gnm.c
plain_search.o : plain_search.c
	gcc -O3 -c plain_search.c
thread_search.o : thread_search.c
	gcc -O3 -c thread_search.c
sequence_handle.o : sequence_handle.c
	gcc -O3 -c sequence_handle.c
event_loop.o : event_loop.c
	gcc -O3 -c event_loop.c
readargs.o : readargs.c
	gcc -O3 -c readargs.c

clean : 
	/bin/rm -f core gnm *.o
