/*
 * SpikeTrain.h
 *
 *  Created on: May 3, 2016
 *      Author: Wilfredo Blanco
 */

#ifndef SPIKETRAIN_H_
#define SPIKETRAIN_H_

#include <vector>
#include <string>     // std::string, std::to_string

using namespace std;

class SpikeTrain {
	unsigned int nNeurons;
	vector< vector<double> > mSpikeTrain; // Spike Train matrix

public:
	SpikeTrain(const unsigned int n);
	virtual ~SpikeTrain();

	void addSpikeTimeToNeuron(const unsigned int n, const double t);

	std::string toString(bool labels) const;
	std::string toString(const unsigned int n, bool labels) const;

	void printToTXTFile(string const strFileName);
	void printToMatlabFile (string const strFileName);

	void printToConsole();
};

#endif /* SPIKETRAIN_H_ */
