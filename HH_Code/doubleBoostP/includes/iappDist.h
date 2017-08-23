

// Declarations
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <vector>
#include <string>     	// std::string, std::to_string

using namespace boost::multiprecision;

void initUniformSpacedIappValues( std::vector<double> &iapp, double deli, double i0);
void initRandDistIappValues( std::vector<double> &iapp, double deli, double i0);

void initRandIappValues(std::vector<cpp_dec_float_100> &iapp);
void initRandIappValues_100(std::vector<cpp_dec_float_100> &iapp);
void initRandIappValues_200(std::vector<cpp_dec_float_100> &iapp);
void initRandIappValues_400(std::vector<cpp_dec_float_100> &iapp);

int writeIappToFile(std::vector<cpp_dec_float_100> &iapp, std::string const fileNameStr);
