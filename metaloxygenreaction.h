#ifndef METALOXYGENREACTION_H
#define METALOXYGENREACTION_H

#include <numeric>
#include "chemicalreaction.h"
#include "elementsfileo.h"

class MetalOxygenReaction : public ChemicalReaction {
	public:
		MetalOxygenReaction();
		MetalOxygenReaction(std::vector<Substance>);
		~MetalOxygenReaction();
		void calculateProducts();
		std::vector<int> getIonicCharges();
		
		
	private:
		//int oxidationNumber;
		//static const std::string O2;
		std::vector<int> ionicCharges;		
};

//const std::string O2 = "Oxygen";

#endif