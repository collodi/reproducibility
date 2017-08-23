

// Declarations

#include <vector>
#include <string>     	// std::string, std::to_string

void initUniformSpacedIappValues( std::vector<double> &iapp, double deli, double i0);
void initRandDistIappValues( std::vector<double> &iapp, double deli, double i0);

void initRandIappValues(std::vector<double> &iapp);
void initRandIappValues_100(std::vector<double> &iapp);
void initRandIappValues_200(std::vector<double> &iapp);
void initRandIappValues_400(std::vector<double> &iapp);

int writeIappToFile(std::vector<double> &iapp, std::string const fileNameStr);
