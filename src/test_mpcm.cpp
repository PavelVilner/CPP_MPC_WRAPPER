#include <iostream>
#include <chrono>

#include "mpw.h"

using namespace std;
using namespace std::chrono;
using namespace mpw_ns;

int main()
{
	mpw_defs::set_d_prec(20);

	mpw_defs::set_mpcw_str_length(21);

	cout << "BEGIN TESTING - Wrap Wrap Wrap!!!" << endl << endl;

	for (uint N = 10; N <= 150; N += 10)
	{
		mpcm qqq = crand(N,N);

		high_resolution_clock::time_point t1 = high_resolution_clock::now();

			mpcm www = qqq / qqq;

		high_resolution_clock::time_point t2 = high_resolution_clock::now();

		cout << N << " --> " << duration_cast<milliseconds>(t2 - t1).count() << endl;
	}

	cout << "FINISH TESTING - Wrap Wrap Wrap!!!" << endl;
}



