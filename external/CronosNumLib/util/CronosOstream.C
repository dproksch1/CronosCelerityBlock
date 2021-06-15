#include "CronosOstream.H"

int NullBuffer::overflow(int c)
{
	return c;
}


cronos_ostream::cronos_ostream(
		const cronos_ostream & costr
):
	std::ostream(costr.rdbuf()),
	buffer(costr.buffer),
	active_rank(costr.active_rank)
{
};

cronos_ostream::~cronos_ostream() {
};

cronos_ostream cronos_ostream::operator()(
		int active_rank
)
{
	return cronos_ostream(buffer, active_rank);
};


