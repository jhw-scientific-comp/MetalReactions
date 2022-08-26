#include <iostream>
#include "metaloxygenreaction.h"


#include "reactionui.h"

using namespace std;



void display(Matrix<int>&);



int main (int argc, char** argv) {

	std::cout << "Text aus Datei:" << std::endl;
	ElementsFileO fileObject;
	
	Atom a_test;	
	cout << "Element mit Ordnungszahl 19:" << endl;
	cout << "gefunden: " << fileObject.selectAtom(19, a_test) << endl;
	cout << "Elementsymbol: " << a_test.getElementSymbol() << endl;
	cout << "Elementname: " << a_test.getElementName() << endl;
	cout << "Valenzelektronen: " << a_test.getValenceElectrons() << endl;
	cout << "Typ: " << a_test.getType() << endl;
	cout << endl;

	cout << "Element mit dem Elementsymbol \"Ni\":" << endl;
	cout << "gefunden: " << fileObject.selectAtom("Ni", a_test) << endl;
	cout << "Elementsymbol: " << a_test.getElementSymbol() << endl;
	cout << "Elementname: " << a_test.getElementName() << endl;
	cout << "Valenzelektronen: " << a_test.getValenceElectrons() << endl;
	cout << "Typ: " << a_test.getType() << endl;
	cout << endl;

	cout << "Element mit dem Namen \"Barium\":" << endl;
	cout << "gefunden: " << fileObject.selectAtom("Barium", a_test) << endl;
	cout << "Elementsymbol: " << a_test.getElementSymbol() << endl;
	cout << "Elementname: " << a_test.getElementName() << endl;
	cout << "Valenzelektronen: " << a_test.getValenceElectrons() << endl;
	cout << "Typ: " << a_test.getType() << endl;
	cout << endl;
	
	cout << "Element mit Ordnungszahl 100:" << endl;
	cout << "gefunden: " << fileObject.selectAtom(100, a_test) << endl;
	cout << endl;

	
	cout << "Element mit dem Namen \"magNesIum\":" << endl;
	cout << "gefunden: " << fileObject.selectAtom("magNesIum", a_test) << endl;
	cout << "Elementsymbol: " << a_test.getElementSymbol() << endl;
	cout << "Elementname: " << a_test.getElementName() << endl;
	cout << "Valenzelektronen: " << a_test.getValenceElectrons() << endl;
	cout << "Typ: " << a_test.getType() << endl;
	cout << endl;


	list<string> atomlist = {"Mg", "3", "Zn", "Calcium"};
	vector<Atom> selAtoms;

	fileObject.selectListOfAtoms(atomlist);
	selAtoms = fileObject.getSelectedAtoms();

	for(auto a : selAtoms) {
		cout << a.getElementSymbol() << ", " << a.getElementName() << ", " << a.getValenceElectrons() << endl;
	}


	cout << "*************** MAIN ***************" << endl;

	ReactionUI ui;
	ui.displaySelection();
	ui.displayChemEq();
	

	cout << "*************** ENDE ***************" << endl;

	return 0;
}



void display(Matrix<int>& mat) {
	const vector<int>* v;	
	unsigned int i {0};
	
	
	v = mat.getEntries();
	
	
	for(auto e : *v) {	
        cout << e << '\t';	
		
		if(i == (mat.getCols() - 1)){
			i = 0;
			cout << endl;
		}
		else{		
		  i++;
		}		
	}

	cout << endl;	
}

