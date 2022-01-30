#ifndef CRONOS_OSTREAM_H
#define CRONOS_OSTREAM_H

#include <iostream>
#ifdef parallel
#include <mpi.h>
#endif


class NullBuffer : public std::streambuf
{
public:
  int overflow(int c);
};

class cronos_ostream: public std::ostream {
//
// Ostream wrapper to restrict the output to certain MPI ranks
//
public:
	cronos_ostream(
			std::ostream & stream, 
			int active_rank
	):
		cronos_ostream(stream.rdbuf(), active_rank)
	{}

	cronos_ostream(
			std::streambuf * streambuff,
			int active_rank
	):
		std::ostream(streambuff),
		buffer(streambuff),
		active_rank(active_rank)
	{
#ifdef parallel
		int isinit;
		MPI_Initialized(&isinit);
		if(isinit == 1) {
			checkRank();
		}
#endif
	};

	cronos_ostream(
			const cronos_ostream & costr
	);

	~cronos_ostream();

	cronos_ostream operator()(
			int active_rank
	);


#ifdef parallel
	/**
	* Check if the rank criterion is met
	* This has to be executed after MPI_Init
	*/
	void checkRank()
	{
		int this_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
		if (this_rank == active_rank) {
			rdbuf(buffer);
		} else {
			rdbuf(&nullBuffer);
		}
	};
#endif


private:
	std::streambuf * buffer;
	int active_rank;

	NullBuffer nullBuffer;
};


// Useful streams
// Have to be initialized in main
extern cronos_ostream ccout;
extern cronos_ostream ccerr;

#endif
