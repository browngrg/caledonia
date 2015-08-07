# The default rule comes first
all:
	@echo "This is a toolset, not a library. Look around."
clean:
	rm -f *.d *.o
veryclean:
	rm -f *.d *.o *.x *.exe

