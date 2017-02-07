bin_PROGRAMS=gsamples

gsamples_CPPFLAGS = -std=c++11 -I$(HTSLIB) -I$(BOOST_HOME) 
gsamples_LDFLAGS = -L$(HTSLIB) -pthread
gsamples_LDADD = -lhts -lz

gsamples_SOURCES = interval_map.h interval_map.cc \
				   fscore.h fscore.cc \
				   binomial.h binomial.cc \
				   config.h config.cc \
				   hit.h hit.cc \
				   boundary.h boundary.cc \
				   region.h region.cc \
				   block.h block.cc \
				   junction.h junction.cc \
				   bundle_base.h bundle_base.cc \
				   bundle.h bundle.cc \
				   assembler.h assembler.cc \
				   gtf.h gtf.cc \
				   fasta.h fasta.cc \
				   util.h util.cc \
				   item.h item.cc \
				   gene.h gene.cc \
				   genome.h genome.cc \
				   transcript.h transcript.cc \
				   main.cc
