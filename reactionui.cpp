#include "reactionui.h"

ReactionUI::ReactionUI() {}


ReactionUI::~ReactionUI() {	
	delete pReaction;
}


void ReactionUI::displaySelection() {
	std::string inputStr;
	std::list<std::string> reactants_strings;
	std::vector<Substance> reactants;
	int inputChoice;
	Atom atom;
	Substance sub;
	
	elementsFileO.selectListOfAtoms("metal");

	cout << "Select metal from list by entering the atomic number (Z), the element name or element symbol:" << endl << endl;

	cout << "Z\t" << "Element name" << endl;
	cout << "--------------------" << endl;

	for(auto a : elementsFileO.getSelectedAtoms()) {
		cout << a.getAtomicNumber() << "\t" << a.getElementName() << " (" << a.getElementSymbol() << ")" << endl;
	}	
	cout << endl;

	elementsFileO.deleteSelectedAtoms();

	cout << ">";	
	cin >> inputStr;	

	reactants_strings.push_back(inputStr);	
	elementsFileO.selectListOfAtoms(reactants_strings);	
	//cout << "Anzahl Atome: " << elementsFileO.getSelectedAtoms().size() << endl;
	sub.addAtoms(elementsFileO.getSelectedAtoms());	
	sub.addNumberOfAtoms(1);
	//cout << "sub Inhalt: " << sub.getAtoms().size() << endl;
	reactants.push_back(sub);
	
	
	cout << endl << endl;
	cout << "Select substance from list from list by entering the number:" << endl << endl;
	
	cout << "No.\t" << "Substance" << endl;
	cout << "--------------------" << endl;
	cout << "1\t" << "Oxygen" << endl; 
	cout << "2\t" << "HCl (Hydrochloric acid)" << endl;
	cout << endl;

	cout << ">";
	cin >> inputChoice;

	switch(inputChoice) {
		case 1:
			pReaction = new MetalOxygenReaction(reactants);
			pReaction->calculateProducts();
			pReaction->balanceEquation();
			break;

		case 2:
			pReaction = new MetalHydrochloricAcidReaction(reactants);
			pReaction->calculateProducts();
			pReaction->balanceEquation();
			break;
	}	

}


void ReactionUI::displayChemEq() {
	unsigned int x_index = 0;

	for(unsigned int i = 0; i < pReaction->getReactants().size(); i++) {
		if(pReaction->getX().at(x_index) != 1) {
			cout << pReaction->getX().at(x_index) << " ";
		}

		for(unsigned int j = 0; j < pReaction->getReactants().at(i).getAtoms().size(); j++) {
			cout << pReaction->getReactants().at(i).getAtoms().at(j).getElementSymbol();
			if(pReaction->getReactants().at(i).getNumberOfAtoms().at(j) != 1) {
				cout << "_{" << pReaction->getReactants().at(i).getNumberOfAtoms().at(j) << "}";
			}			
		}
		
		if(i != pReaction->getReactants().size() - 1) {
			cout << " + ";
		}
		else {
			cout << " -> ";
		}

		x_index++;
	}

	for(unsigned int i = 0; i < pReaction->getProducts().size(); i++) {
		if((x_index < pReaction->getX().size()) && (pReaction->getX().at(x_index) != 1)) {
			cout << pReaction->getX().at(x_index) << " ";
		}
		else if((x_index >= pReaction->getX().size()) && (pReaction->getX0() != 1)) {
			cout << pReaction->getX0() << " ";
		}

		for(unsigned int j = 0; j < pReaction->getProducts().at(i).getAtoms().size(); j++) {
			cout << pReaction->getProducts().at(i).getAtoms().at(j).getElementSymbol();
			if(pReaction->getProducts().at(i).getNumberOfAtoms().at(j) != 1) {
				cout << "_{" << pReaction->getProducts().at(i).getNumberOfAtoms().at(j) << "}";
			}			
		}
		
		if(i != pReaction->getProducts().size() - 1) {
			cout << " + ";
		}

		x_index++;
	}
	cout << endl;
}


