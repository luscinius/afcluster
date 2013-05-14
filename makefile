CXXFLAGS=-O3
afcluster: afcluster.cc freqs.cc  freqs.hh dists.hh dists.cc centroid.hh centroid.cc hc.cc  hc.hh  em.cc  em.hh  seqio.cc  seqio.hh
	g++ ${CXXFLAGS} -Wall afcluster.cc em.cc hc.cc dists.cc centroid.cc freqs.cc seqio.cc -o afcluster
