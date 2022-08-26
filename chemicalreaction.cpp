#include "chemicalreaction.h"


ChemicalReaction::ChemicalReaction() {}

ChemicalReaction::~ChemicalReaction() {
	reactants.clear();
	products.clear();
}



//TODO: Diese Methode muss die Faktoren zusammen mit den zugehörigen Substanzen separat speichern!!!!

void ChemicalReaction::balanceEquation() {
	std::vector<std::string> elemSymbSet;
	Matrix<int>* reactionCoefficientMatrix;
	std::vector<std::string>::iterator iter;
	std::vector<Substance>::iterator iter_sub;
	unsigned int index;

	unsigned int indexCombi;

	std::vector<Substance> productVectorSubs;
	std::vector<int>* productAtomsVector;

	unsigned int start = 0; 
	unsigned int lengthCombi;
	std::vector<Substance> combination; //TODO: vielleicht zu list aendern?????
	std::vector<std::vector<Substance>> combinationsVec; 
	unsigned int fixedElements = 0;
	unsigned int level = 0;

	int det_reactionCoefficientMatrix = 0;

	std::vector<std::pair<int,int>> bpg_edges;
	std::vector<std::pair<int,int>> edges_temp;
	BalancedBipartiteGraph* pBBPG;
	std::vector<int> rowRearrangement;

	/* std::vector<double> x;
	double x_0 = 0.0; */

	
	for(auto s : reactants) {
		for(auto a : s.getAtoms()) {
			iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), a.getElementSymbol());

			if(iter == elemSymbSet.end()) {
				elemSymbSet.push_back(a.getElementSymbol());
			}
		}
	}

	for(auto s : products) {
		for(auto a : s.getAtoms()) {
			iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), a.getElementSymbol());

			if(iter == elemSymbSet.end()) {
				elemSymbSet.push_back(a.getElementSymbol());
			}
		}
	}

	std::cout << "Elemente: " << elemSymbSet.size() << std::endl;

	reactionCoefficientMatrix = new Matrix<int> (elemSymbSet.size(), elemSymbSet.size());

	productAtomsVector = new std::vector<int> (elemSymbSet.size(), 0);

	for(auto s : elemSymbSet) {
		std::cout << s << std::endl;
	}



	std::cout << "Matrix aufstellen:" << std::endl;

	for(unsigned int k = 0; k < reactants.size(); k++) {			
		for(unsigned int l = 0; l < reactants.at(k).getAtoms().size(); l++) {
			iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), reactants.at(k).getAtoms().at(l).getElementSymbol());
			
			index = std::distance(elemSymbSet.begin(), iter);

			reactionCoefficientMatrix->setEntry(index, k, reactants.at(k).getNumberOfAtoms().at(l));

			bpg_edges.push_back(std::make_pair(index, k));
		}
	}


	
	if(reactants.size() < elemSymbSet.size()) {			

		lengthCombi = elemSymbSet.size() - reactants.size();
		std::cout << "lengthCombi: " << lengthCombi << std::endl;
		findCombinations(start, lengthCombi, combination, combinationsVec, fixedElements, level);
		
		for(unsigned int i = 0; (i < combinationsVec.size()) && (det_reactionCoefficientMatrix == 0); i++) {
			for(unsigned int l = 0; l < lengthCombi; l++) {
				reactionCoefficientMatrix->setCol(elemSymbSet.size() - lengthCombi + l, 0);			
			}
	
			edges_temp.clear();

			for(unsigned int j = 0; j < combinationsVec.at(i).size(); j++) {				
				for(unsigned int k = 0; k < combinationsVec.at(i).at(j).getAtoms().size(); k++) {
					
					iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), combinationsVec.at(i).at(j).getAtoms().at(k).getElementSymbol());
			
					index = std::distance(elemSymbSet.begin(), iter);

					reactionCoefficientMatrix->setEntry(index, elemSymbSet.size() - combinationsVec.at(i).size() + j, -combinationsVec.at(i).at(j).getNumberOfAtoms().at(k));
					
					edges_temp.push_back(std::make_pair(index, elemSymbSet.size() - combinationsVec.at(i).size() + j));
				}
				
			}
			
			det_reactionCoefficientMatrix = det(*reactionCoefficientMatrix);
			indexCombi = i;

			

			if(det_reactionCoefficientMatrix != 0) {
				for(auto p : edges_temp) {
					bpg_edges.push_back(p);
				}

				combination = combinationsVec.at(i);
			}
		}		


		//verbleibende Substanzen fuer die rechte Seite ermitteln
		for(auto s : products) {			
			iter_sub = std::find(combinationsVec.at(indexCombi).begin(), combinationsVec.at(indexCombi).end(), s);
			
			if(iter_sub == combinationsVec.at(indexCombi).end()) {
				productVectorSubs.push_back(s);
			}
		}

		
		//rechte Seite des GLS aufstellen
		for(auto s : productVectorSubs) {		
			for(unsigned int i = 0; i < s.getAtoms().size(); i++) {			
				iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), s.getAtoms().at(i).getElementSymbol());
				
				index = std::distance(elemSymbSet.begin(), iter);				

				productAtomsVector->at(index) += s.getNumberOfAtoms().at(i);
			}
		}
		
	}
	else {		
		//rechte Seite des GLS aufstellen
		for(auto s : products) {		
			for(unsigned int i = 0; i < s.getAtoms().size(); i++) {			
				iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), s.getAtoms().at(i).getElementSymbol());
				
				index = std::distance(elemSymbSet.begin(), iter);				

				productAtomsVector->at(index) += s.getNumberOfAtoms().at(i);
			}
		}
	}


//TODO: Folgender Code darf nur ausgefuehrt werden, falls die Determninante der Reaktionsmatrix ungleich 0 ist.
//Ansonsten wäre nicht sichergestellt, dass das Maximum-Matching ein perfektes Matching wäre.

	//Bipartiten Graph aus der Matrix aufstellen und Maximum-Matching finden (ist hier automatisch ein perfektes Matching)
	pBBPG = new BalancedBipartiteGraph(elemSymbSet.size(), bpg_edges);
	rowRearrangement = pBBPG->findMaxMatch();	
	

	//Matrix und rechte Seite entsprechend des Matchings umordnen
	std::vector<unsigned int> newIndex;
	std::vector<unsigned int>::iterator it_newIndex;
	unsigned int product_temp;
	
	for(unsigned int i = 0; i < elemSymbSet.size(); i++) {
		newIndex.push_back(i);
	}
	
	
	
	for(unsigned int k = 0; k < rowRearrangement.size(); k++) {
		reactionCoefficientMatrix->swapRows(newIndex.at(k), rowRearrangement.at(k));

		product_temp = productAtomsVector->at(newIndex.at(k));
		productAtomsVector->at(newIndex.at(k)) = productAtomsVector->at(rowRearrangement.at(k));
		productAtomsVector->at(rowRearrangement.at(k)) = product_temp;

		it_newIndex = std::find(newIndex.begin(), newIndex.end(), rowRearrangement.at(k));
		index = std::distance(newIndex.begin(), it_newIndex);
		newIndex.at(k) = rowRearrangement.at(k);
		newIndex.at(index) = k;
	}


	//Lösen des LGS

	x = solveSLE_ICRB(*reactionCoefficientMatrix, *productAtomsVector, x_0);
	

	std::cout << "Reaktionsmatrix:" << std::endl;
	for(unsigned int i = 0; i < reactionCoefficientMatrix->getRows(); i++) {
		for(unsigned int j = 0; j < reactionCoefficientMatrix->getCols(); j++) {
			std::cout << reactionCoefficientMatrix->getEntry(i,j) << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "Combinations:" << std::endl;
	for(auto s : combination) {
		for(unsigned int i = 0; i < s.getAtoms().size(); i++) {
			std::cout << s.getAtoms().at(i).getElementSymbol() << "_{" << s.getNumberOfAtoms().at(i) << "}";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;



	//Umordnen des Produktvektors damit die Reihenfolge dem Loesungsvektor x angeglichen wird
	std::vector<Substance> products_temp(products);
	
	for(auto s : combination) {
		iter_sub = std::find(products_temp.begin(), products_temp.end(), s);
		products_temp.erase(iter_sub);
	}
	
	products.clear();
	
	for(auto s : combination) {
		products.push_back(s);
	}

	for(auto s : products_temp) {
		products.push_back(s);
	}





	std::cout << "products:" << products.size() << std::endl;
	for(unsigned int i = 0; i < products.size(); i++) {
		for(unsigned int j = 0; j < products.at(i).getAtoms().size(); j++){
			std::cout << products.at(i).getAtoms().at(j).getElementSymbol() << "_{" << products.at(i).getNumberOfAtoms().at(j) << "}";			
		}
		std::cout << std::endl;
	}
	
	
	std::cout << "Produktvektor:" << std::endl;
	for(unsigned int i = 0; i < productAtomsVector->size(); i++) {
		std::cout << productAtomsVector->at(i) << std::endl;
	}	


	std::cout << "Size of combinationsVec: " << combinationsVec.size() << std::endl;
	for(auto v : combinationsVec) {
		for(auto s : v) {			
			for(unsigned int k = 0; k < s.getAtoms().size(); k++) {
				std::cout << "(" << s.getAtoms().at(k).getElementSymbol() << "," << s.getNumberOfAtoms().at(k) << ")";
			}
			std::cout << "  ";
		}
		std::cout << std::endl;
	}


	std::cout << "BPG edges:" << std::endl;
	for(auto v : bpg_edges) {
		std::cout << "(" << v.first << ", " << v.second << ")" << std::endl;
	}

	
	std::cout << "Loesung x:" << std::endl;
	for(auto e : x) {
		std::cout << e << std::endl;
	}
	std::cout << "Loesung x_0: " << x_0 << std::endl;


	std::cout << "-------------------------" << std::endl;	

	
	delete pBBPG;
	delete reactionCoefficientMatrix;
	delete productAtomsVector;
}





std::vector<Substance> ChemicalReaction::getReactants() {
	return reactants;
}


std::vector<Substance> ChemicalReaction::getProducts() {
	return products;
}


void ChemicalReaction::addReactants(std::vector<Substance> reactants_) {
	for(auto r : reactants_) {
		reactants.push_back(r);
	}
}


void ChemicalReaction::addReactants(Substance reactant_) {
	reactants.push_back(reactant_);
}


void ChemicalReaction::addProducts(std::vector<Substance> products_) {
	for(auto p : products_) {
		products.push_back(p);
	}
}


void ChemicalReaction::addProducts(Substance product_) {
	products.push_back(product_);
}



//TODO: evtl combination als queue??????
//finde alle möglichen Kombinationen der Produktsubstanzen der Laenge lengthCombi
void ChemicalReaction::findCombinations(unsigned int startElement, unsigned int lengthCombi, std::vector<Substance>& combination, std::vector<std::vector<Substance>>& combinationsVec, unsigned int& fixedElements, unsigned int& level) {
		
	for(unsigned int i = startElement; i < products.size() - lengthCombi + level + 1; i++) {		
		combination.push_back(products.at(i));	
		fixedElements++;

		//std::cout << "level: " << level << ", schleife: " << i << ", fixedElements: " << fixedElements << ", combination: " << combination.size() << std::endl;

		if((lengthCombi - fixedElements) > 0) {	
			level++;
			findCombinations(i + 1, lengthCombi, combination, combinationsVec, fixedElements, level);
		}
		else if(i < products.size() - 1) {			 
			combinationsVec.push_back(combination);
			//std::cout << "----------------" << std::endl;
			combination.pop_back();
			fixedElements--;
		}
		else {
			combinationsVec.push_back(combination);
			//std::cout << "----------------" << std::endl;						
		}
	}
	for(unsigned int k = 0; k < level; k++) {
		combination.pop_back();
		fixedElements--;
	}

	level--;

}


//FEHLER: IRGENDWAS FUNKTIONIERT HIER NICHT!!!!!!!!!!!!!!!!!!
//PROBLEM ERKANNT?: 
//- ein Vektor als Rueckgabewert gibt keine Kopie zurueck, sondern der Pointer des oben zugewiesenen Vektors wird auf die Adresse des zurueckgewiesenen Vektors gesetzt.
//- Speicherlecks bei der Parameteruebergabe der Matrix oder des Vektors
std::vector<double> ChemicalReaction::solveSLE_ICRB(const Matrix<int>& rcMatrix, std::vector<int> prodVector, double& x_0_var) {
	unsigned int n = prodVector.size();
	std::vector<double> x_vec(n, 0);
	std::vector<double> r_vec(n, 0);
	double s = 0.0;
	int stop = 0;
	unsigned int i_max = 100;
	double tol = 0.0001;	
	double r = 0.0;

	

	for(unsigned int k = 0; (k < i_max) && (stop == 0); k++) {
		for(unsigned int i = 0; i < n; i++) {
			s = 0.0;
			r_vec = x_vec;

			for(unsigned int j = 0; j < n; j++) {
				if(j != i) {
					s = s + static_cast<double>(rcMatrix.getEntry(i, j) * x_vec.at(j));
				}
			}
			x_vec.at(i) = (static_cast<double>(prodVector.at(i)) - s)/static_cast<double>(rcMatrix.getEntry(i, i));

			if(std::fabs((r_vec.at(i) - x_vec.at(i))/r_vec.at(i)) <= tol) {
				stop = 1;
			}
		}
	}


	std::cout << "Zwischenergebnis x:" << std::endl;
	for(auto e : x_vec) {
		std::cout << e << std::endl;
	}


	r_vec = x_vec;
	stop = 0;

	s = 0.0;
	for(auto e : r_vec) {
		s += e;
	}

	for(unsigned int k = 1; (k < i_max) && (stop == 0); k++) {
		x_0_var = static_cast<double>(k);
		r = std::fabs((k * s - std::round(k * s))/s);
		
		if(r <= tol) {
			stop = 1;
		}
	}

	std::cout << "x_0_var: " << x_0_var << std::endl;

	for(unsigned int i = 0; i < x_vec.size(); i++) {
		x_vec.at(i) = x_0_var * x_vec.at(i);
	}

	return x_vec;
}


std::vector<double> ChemicalReaction::getX() const {
	return x;
}
	

double ChemicalReaction::getX0() const {
	return x_0;
}
