

// Declarations

#include <vector>
#include <string>     	// std::string, std::to_string

void initUniformSpacedIappValues( std::vector<long double> &iapp, long double deli, long double i0);
void initRandDistIappValues( std::vector<long double> &iapp, long double deli, long double i0);

void initRandIappValues(std::vector<long double> &iapp);
void initRandIappValues_100(std::vector<long double> &iapp);
void initRandIappValues_200(std::vector<long double> &iapp);
void initRandIappValues_400(std::vector<long double> &iapp);

int writeIappToFile(std::vector<long double> &iapp, std::string const fileNameStr);
