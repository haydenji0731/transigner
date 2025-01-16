ifeq ($(shell command -v gcc 2>/dev/null),)
  CC=clang
else
  CC=gcc
endif

LDFLAGS=-lm
SOURCES=transigner/em.c transigner/xxhash.c
OBJECTS=$(SOURCES:.c=.o)
TARGET=em

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
