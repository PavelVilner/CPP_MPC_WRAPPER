#include <iostream>

#include "mpw.h"

using namespace std;
using namespace mpw_ns;

int main()
{
	mpcw::set_b_prec(1000);

	cout << "BEGIN TESTING - Wrap Wrap Wrap!!!" << endl;

	mpcw qqq = mpcw(-0.02,-0.05);

	cout << qqq.str(21) << endl;

	cout << (-0.6 + qqq).atanh().str(50) << endl;

	cout << (-0.6 + qqq).atanh().real().str(50) << endl;

	cout << (-0.6 + qqq).atanh().imag().str(50) << endl;

	cout << "FINISH TESTING - Wrap Wrap Wrap!!!" << endl;
}



