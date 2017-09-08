#include <iostream>
#include <sstream>

class Thing : public std::ostringstream
{
public:
	Thing() : std::ostringstream() {}
	virtual ~Thing() { std::cerr << str(); }
};

int main(int argc, const char * argv[]) {
	Thing() << "Hello" << std::endl;
	return 0;
}
