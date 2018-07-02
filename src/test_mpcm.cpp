#include <iostream>
#include <chrono>

#include "mpw.h"

using namespace std;
using namespace std::chrono;
using namespace mpw_ns;

int main()
{
	mpw_defs::set_d_prec(200);

	mpw_defs::set_mpcw_str_length(21);

	mpw_defs::set_use_Householder_prec(false);

	cout << "BEGIN TESTING - Wrap Wrap Wrap!!!" << endl << endl;


//  TEST MATRIX INVERSAL
//	for (uint N = 10; N <= 150; N += 10)
//	{
//		mpcm qqq = crand(N,N);
//
//		high_resolution_clock::time_point t1 = high_resolution_clock::now();
//
//			mpcm www = qqq / qqq;
//
//		high_resolution_clock::time_point t2 = high_resolution_clock::now();
//
//		cout << N << " --> " << duration_cast<milliseconds>(t2 - t1).count() << endl;
//	}

// TEST HOUSEHOLDER TRANSFORMATION
	/*cout << "HOUSEHOLDER TESTING" << endl;
	for (uint N = 10; N <= 50; N += 10)
	{
		mpcm qqq = crand(N,N);

		high_resolution_clock::time_point t1 = high_resolution_clock::now();

			std::tuple<mpcm, mpcm> www = qqq.Householder_double_sided();

		high_resolution_clock::time_point t2 = high_resolution_clock::now();

		cout << "Max diff: " << str(max_diff(std::get<1>(www)*std::get<0>(www)*std::get<1>(www).H(), qqq)) << endl;

		cout << "N = " << N << ": --> " << duration_cast<milliseconds>(t2 - t1).count() << "ms" << endl;
	}*/

// TEST HERMITIAN EIG AND JACOBI
	cout << "HERMITIAN EIG TESTING" << endl;
	for (uint N = 10; N <= 50; N += 10)
	{
		mpcm qqq = crand(N, N);

		qqq = qqq + qqq.H();

		high_resolution_clock::time_point t1 = high_resolution_clock::now();

			std::tuple<mpcm, mpcm> www = qqq.Hermitian_eig();

		high_resolution_clock::time_point t2 = high_resolution_clock::now();

		mpcm L = std::get<0>(www);

		mpcm V = std::get<1>(www);

		cout << "Max diff: " << str(max_diff(V*L*V.H(), qqq)) << endl;

		cout << "N = " << N << ": --> " << duration_cast<milliseconds>(t2 - t1).count() << "ms" << endl;
	}

	cout << endl << "FINISH TESTING - Wrap Wrap Wrap!!!" << endl;
}
