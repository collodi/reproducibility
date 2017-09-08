#include <iostream>
using namespace std;

int f(const int& value)
{
	static int result = 0;
	return result += value;
}

int main()
{
	cout << f(10) << ", " << f(f(10)) << ", " << f(f(f(10))) << endl;
	return 0;
}
