PROG = BARP

all : $(PROG) clean


HTSDIR = third_party/htslib-1.16
HTSLIBDIR = $(HTSDIR)/htslib
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSVERSION = $(HTSDIR)/version.h

LIBBARP = src/libbarp.a

ZLIBDIR = third_party/zlib-1.2.11
include $(ZLIBDIR)/zlib.mk
LIBZ = $(ZLIBDIR)/libz.a


CC		= gcc
CFLAGS	= -Wall -O0 -g -D_FILE_OFFSET_BITS=64
DEBUGFLAGS = -fsanitize=address -fno-omit-frame-pointer -O0 -g 
INCLUDES = -I src -I $(HTSLIBDIR) -I . -I $(ZLIBDIR)
LIBS = -lm -pthread

PACKAGE_VERSION = 0.1.1 

.SUFFIXES: .c.o 

.c.o: 
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@ 

LIB_OBJ = src/read_json.o \
	src/config_hash.o \
	src/kson.o \
	src/parse_config.o \
	src/sam_opts.o \
	src/parse_read.o \
	src/parse_samfield.o \
	src/merge_interval.o

libbarp.a: $(LIB_OBJ)
	@-rm -f src/$@
	@$(AR) -rcs src/$@ $(LIB_OBJ)

BARP: libbarp.a $(LIBZ) $(HTSLIB)
	@$(CC) -I third_party/curl-8.1.2/include $(CFLAGS) $(INCLUDES) src/main.c  src/libbarp.a -o $@ $(HTSLIB) $(LIBZ) $(LIBS) 

clean:
	@-rm -f src/*.o src/*.a 

.PHONY: clean 
