/*
 * Panagiotis Skrimponis
 */

#include "cs_dc_sweep.h"
#include "cs_ac_sweep.h"
#include "cs_component.h"
#include "cs_transient.h"
#include "cs_headerdef.h"

int main(int argc, char * argv[])
{
	std::ifstream iFile;
	std::ofstream oFile;
	std::stringstream iLine;
	std::string token, line, component_name, component_plus_node, component_minus_node, analysis_sweep;
	std::vector < std::unique_ptr<Component> > ComponentVector;
	std::vector <TransientComponent> TransientComponentVector;
	std::vector <ACComponent> ACComponentVector;
	std::vector < std::unique_ptr<DCSweepAnalysis>> DCSweepAnalysisVector;
	std::vector < std::unique_ptr<ACSweepAnalysis>> ACSweepAnalysisVector;
	std::vector < std::unique_ptr<TransientAnalysis>> TransientAnalysisVector;
	std::unordered_map < std::string, int > nodeHash, ctrlHash, componentHash;
	std::unordered_map < int, std::string > nodeHash_reverse;
	std::size_t last, found;
	char component_type;
	unsigned int i=0, j=0, k=0, nodeID=1, componentID=0, component_branch=0, branchID=0, nonZero=0, vectorSize=0, tranZero=0;
	bool spd=false, iter=false, sparse=false, method=false, ac=false;
	double analysis_start=0.0, analysis_stop=0.0, analysis_step=0.0, component_value=0.0, component_mag=0.0, component_phase=0.0, itol = 1e-3, g=0.0;
	
	// Initialize Hash-Tables
	nodeHash["0"] = nodeHash["GND"] = -1;
	ctrlHash[".DC"]		 = DC;
	ctrlHash[".AC"]		 = AC;
	ctrlHash[".TRAN"]	 = TRAN;
	ctrlHash[".OPTIONS"] = OPTIONS;
	ctrlHash[".PLOT"] 	 = PLOT;
	ctrlHash[".PRINT"]	 = PRINT;
	ctrlHash["SPD"] 	 = SPD;
	ctrlHash["ITER"]	 = ITER;
	ctrlHash["SPARSE"] 	 = SPARSE;
	ctrlHash["METHOD"]	 = METHOD;
	ctrlHash["ITOL"]	 = ITOL;

	// Check if there are enough arguments
	if (argc != 2)
	{
		std::cout << "\033[1;31mERROR: not enough arguments\033[0m\n";
		return -1;
	}

	// Check if you can open the input file
	iFile.open(argv[1]);
	if (!iFile.is_open())
	{
		std::cout << "\033[1;31mERROR: cannot open the file\033[0m" << argv[1] << "for reading\n";
		return -1;
	}

	iFile >> token;
	str_toupper(token)
	do {
		switch(token[0])
		{
			case('V'):
			case('I'):
			case('R'):
			case('C'):
			case('L'):
				component_type = token[0];
				component_name = token;
				iFile >> component_plus_node >> component_minus_node;
				str_toupper(component_plus_node)
				str_toupper(component_minus_node)
				componentHash[component_name] = componentID++;
				if (!nodeHash[component_plus_node])
				{
					nodeHash[component_plus_node] = nodeID;
					nodeHash_reverse[nodeID++] = component_plus_node;
				}
				if (!nodeHash[component_minus_node])
				{
					nodeHash[component_minus_node] = nodeID;
					nodeHash_reverse[nodeID++] = component_minus_node;
				}		
				if(token[0]=='R') 
				{
					nonZero += (((nodeHash[component_plus_node] != -1) & (nodeHash[component_minus_node] != -1)) << 1) + (nodeHash[component_plus_node] != -1) + (nodeHash[component_minus_node] != -1);
				}
				else if ((token[0]=='V') || (token[0]=='L'))
				{
					nonZero += ((nodeHash[component_plus_node] != -1) << 1) + ((nodeHash[component_minus_node] != -1) << 1);
					component_branch = branchID++;
					tranZero += (token[0]=='L');
				}
				else if (token[0]=='C')
				{
					tranZero  += (((nodeHash[component_plus_node] != -1) & (nodeHash[component_minus_node] != -1)) << 1) + (nodeHash[component_plus_node] != -1) + (nodeHash[component_minus_node] != -1);
				}
				getline(iFile, token);
				break;
			case('.'):
				switch(ctrlHash[token])
				{
					case(DC): // Create a new DC Anaylsis
						iFile >> token >> analysis_start >> analysis_stop >> analysis_step;
						str_toupper(token)
						component_type = token[0];
						component_name = token;
						DCSweepAnalysisVector.push_back(std::unique_ptr<DCSweepAnalysis> (new DCSweepAnalysis(component_type, component_name, analysis_start, analysis_step, analysis_stop)));
						break;
					case(AC): // Create a new AC Anaylsis
						iFile >> analysis_sweep >> analysis_step >> analysis_start >> analysis_stop;
						ACSweepAnalysisVector.push_back(std::unique_ptr <ACSweepAnalysis> (new ACSweepAnalysis(analysis_sweep, analysis_start, analysis_step, analysis_stop)));
						break;
					case(TRAN): // Create a new Transient Anaylsis
						iFile >> analysis_step >> analysis_stop;
						TransientAnalysisVector.push_back(std::unique_ptr <TransientAnalysis> (new TransientAnalysis(analysis_step, analysis_stop)));
						break;
					case(OPTIONS):
						iFile.seekg(- token.size(), iFile.cur);
						read_line()
						iLine.str(token); iLine.clear();
						while(!iLine.eof())
						{
							iLine >> token;
							if(ctrlHash[token] == SPD) spd = true;
							else if(ctrlHash[token] == SPARSE) sparse = true;
							else if(ctrlHash[token] == ITER) iter = true;
							else if(ctrlHash[token] == ITOL) iLine >> itol;
							else if(ctrlHash[token] == METHOD) { iLine >> token; method = (token == "BE"); }
						}
						break;
					case(PRINT):
					case(PLOT):
						iFile >> token;
						switch(token[0])
						{
							case('D'):
								if (DCSweepAnalysisVector.size() == 0) break;
								read_line();
								iLine.str(token); iLine.clear();
								do {
									iLine >> token >> component_name;
									DCSweepAnalysisVector[DCSweepAnalysisVector.size()-1]->add_node(component_name);
								} while(!iLine.eof());
								break;
							case('T'):
								if (TransientAnalysisVector.size() == 0) break;
								read_line();
								iLine.str(token); iLine.clear();
								do {
									iLine >> token >> component_name;
									TransientAnalysisVector[TransientAnalysisVector.size()-1]->add_node(component_name);
								} while(!iLine.eof());
								break;
							case('A'):
								if (ACSweepAnalysisVector.size() == 0) break;
								read_line();
								iLine.str(token); iLine.clear();
								do {
									iLine >> token >> component_name;
									ACSweepAnalysisVector[ACSweepAnalysisVector.size()-1]->add_node(component_name);
								} while(!iLine.eof());
								break;
							default:
								read_line();
								iLine.str(token); iLine.clear();
								
								do {
									iLine >> token >> component_name;
									if (  DCSweepAnalysisVector.size() != 0) DCSweepAnalysisVector[DCSweepAnalysisVector.size()-1]->add_node(component_name);
									if (TransientAnalysisVector.size() != 0) TransientAnalysisVector[TransientAnalysisVector.size()-1]->add_node(component_name);
									if (  ACSweepAnalysisVector.size() != 0) ACSweepAnalysisVector[ACSweepAnalysisVector.size()-1]->add_node(component_name);
								} while(!iLine.eof());
								
								break;
						}
						break;
				}
				break;
			default: getline(iFile, token);
		}
		iFile >> token;
		str_toupper(token)
	} while(!iFile.eof());
	iFile.clear();
	iFile.seekg(0, iFile.beg);
	nodeHash["0"] = nodeHash["GND"] = 0;
	nodeHash_reverse[0] = "0";
	vectorSize = nodeID + branchID;
	if(sparse)
	{
		int idx=0, tdx=0;
		double * CS_b,* CS_x;
		CS_b = new double[vectorSize-1]; for(i = 0; i < (vectorSize - 1); ++i) CS_b[i] = 0.0;
		CS_x = new double[vectorSize-1]; for(i = 0; i < (vectorSize - 1); ++i) CS_x[i] = 0.0;
		std::complex<double>* MNA_AC_b = new std::complex<double>[vectorSize];
		SparseMatrix<double> CS_A("Sparse MNA DC Matrix", nodeID, branchID, nonZero);
		SparseMatrix<double> CS_C("Sparse MNA Transient Matrix", nodeID, branchID, tranZero);
		branchID = 0;
		iFile >> token;
		str_toupper(token)
		while (!iFile.eof())
		{			
			switch (token[0]) 
			{
				case('V'):
				case('I'):
				case('R'):
				case('C'):
				case('L'):
					component_type = token[0];
					component_name = token.substr(1);
					read_line()
					ac = ((found = token.find(" AC ")) != std::string::npos);
					iLine.str(token.substr(0, found)); iLine.clear();

					iLine >> component_plus_node >> component_minus_node >> component_value;
					str_toupper(component_plus_node)
					str_toupper(component_minus_node)

					if(ac)
					{
						iLine.str(token.substr(found, std::string::npos));
						iLine >> token >> component_mag >> component_phase;
						ACComponentVector.push_back(ACComponent(component_type, component_name, component_value, component_mag, component_phase, nodeHash[component_plus_node], nodeHash[component_minus_node], branchID));
					}

					i = nodeHash[component_plus_node] - 1;
					j = nodeHash[component_minus_node] - 1;
					if (!iLine.eof())
					{
						TransientComponentVector.push_back(TransientComponent(component_type, component_name, iLine, component_value, nodeHash[component_plus_node], nodeHash[component_minus_node], branchID));
					}
					switch (component_type) 
					{
						case('V'):
							ComponentVector.push_back(std::unique_ptr <Component> (new Component(component_type, component_name, nodeHash[component_plus_node], nodeHash[component_minus_node], component_value)));
							k = (branchID + nodeID - 1);
							ComponentVector[ComponentVector.size()-1]->set_branch(branchID++);
							CS_b[k] += component_value;
							if((i != -1) && (j != -1))
							{
								CS_A(idx++, k, i,  1);
								CS_A(idx++, i, k,  1);
								CS_A(idx++, k, j, -1);
								CS_A(idx++, j, k, -1);
							}
							else if(i != -1)
							{
								CS_A(idx++, k, i,  1);
								CS_A(idx++, i, k,  1);
							}
							else if(j != -1)
							{
								CS_A(idx++, k, j, -1);
								CS_A(idx++, j, k, -1);
							}
							break;
						case('I'):
							ComponentVector.push_back(std::unique_ptr <Component> (new Component(component_type, component_name, nodeHash[component_plus_node], nodeHash[component_minus_node], component_value)));
							if (i != -1) CS_b[i] -= component_value;
							if (j != -1) CS_b[j] += component_value;
							break;
						case('R'):
							g = (double) 1.0/component_value;
							if((i != -1) && (j != -1))
							{
								CS_A(idx++, j, j,  g);
								CS_A(idx++, i, i,  g);
								CS_A(idx++, i, j, -g);
								CS_A(idx++, j, i, -g);
							}
							else if(i != -1)
							{
								CS_A(idx++, i, i,  g);
							}
							else if(j != -1)
							{
								CS_A(idx++, j, j,  g);
							}
							break;
						case('C'):
							if((i != -1) && (j != -1))
							{
								CS_C(tdx++, j, j,  component_value);
								CS_C(tdx++, i, i,  component_value);
								CS_C(tdx++, i, j, -component_value);
								CS_C(tdx++, j, i, -component_value);
							}
							else if(i != -1)
							{
								CS_C(tdx++, i, i,  component_value);
							}
							else if(j != -1)
							{
								CS_C(tdx++, j, j,  component_value);
							}
							break;						
						case('L'):
							k = (branchID + nodeID - 1); branchID++;
							if((i != -1) && (j != -1))
							{
								CS_C(tdx++, k, k,  -component_value);
								CS_A(idx++, k, i,  1);
								CS_A(idx++, i, k,  1);
								CS_A(idx++, k, j, -1);
								CS_A(idx++, j, k, -1);
							}
							else if(i != -1)
							{
								CS_C(tdx++, k, k,  -component_value);
								CS_A(idx++, k, i,  1);
								CS_A(idx++, i, k,  1);
							}
							else if(j != -1)
							{
								CS_C(tdx++, k, k,  -component_value);
								CS_A(idx++, k, j, -1);
								CS_A(idx++, j, k, -1);
							}
							break;
					}
					break;
				default: getline(iFile, token);
			}
			iFile >> token;
			str_toupper(token)
		}
		CS_A.compressMatrix();
		CS_C.compressMatrix();
		if(ac)
		{
			cs_ci * csi_miami = cs_ci_add(cs_i_complex(CS_A.matrix(), true), cs_i_complex(CS_C.matrix(), false), 1.0, 1.0);
			ComplexSparseMatrix<std::complex<double>> CS_AC ("MNA AC Matrix", nodeID, branchID, csi_miami);	
			for(auto &it : ACComponentVector)
			{
				if (it.type() == 'V')
					MNA_AC_b[it.branch() + nodeID] = it.value();
				else
				{
					if(it.plus () != 0)
						MNA_AC_b[it.plus ()-1] -= it.value();
					if(it.minus () != 0)
						MNA_AC_b[it.minus()-1] += it.value();
				}
			}
			for(auto &it : ACSweepAnalysisVector) 
			{
				it->update(nodeHash);
				it->sp_analyse(CS_AC, MNA_AC_b, itol, iter, spd, ACComponentVector);
			}
		}

		for(auto &it : TransientAnalysisVector)
		{
			it->update(nodeHash);
			it->sp_analyse(CS_A, CS_C, CS_x, CS_b, itol, iter, spd, method, TransientComponentVector);
		}

		if(iter)
		{
			spd ? CS_A.CG(CS_x, CS_b, itol) : CS_A.BiCG(CS_x, CS_b, itol);
		}
		else
		{
			spd ? CS_A.Cholesky() : CS_A.LU();
			spd ? CS_A.solveSPD(CS_x, CS_b) : CS_A.solve(CS_x, CS_b);
		}
		oFile.open("DC_Operating_Point_Analysis.txt");
		for(i=0;i<(vectorSize-1-branchID);++i) oFile << nodeHash_reverse[i+1] << " " << CS_x[i] << "\n";
		oFile.close();

		for(auto &it : DCSweepAnalysisVector) 
		{
			it->update(ComponentVector[componentHash[it->name()]]->plus(), ComponentVector[componentHash[it->name()]]->minus(), ComponentVector[componentHash[it->name()]]->branch(), ComponentVector[componentHash[it->name()]]->value(), nodeHash);
			it->sp_analyse(CS_A, CS_b, itol, iter, spd);
		}
	}
	else
	{
		Matrix <std::complex<double>> MNA_AC("MNA AC Matrix", nodeID, branchID);
		Matrix <double> MNA_DC ("MNA DC Matrix", nodeID, branchID);
		Matrix <double> MNA_TRAN ("MNA Transient Matrix", nodeID, branchID);
		double * MNA_b = new double[vectorSize];
		double * MNA_x = new double[vectorSize];
		std::complex<double>* MNA_AC_b = new std::complex<double>[vectorSize];

		for(i = 0; i < vectorSize; ++i) MNA_x[i] = 0.0;
		branchID = 0;
		iFile >> token;
		str_toupper(token)
		do {			
			switch (token[0]) 
			{
				case('V'):
				case('I'):
				case('R'):
				case('C'):
				case('L'):
					component_type = token[0];
					component_name = token.substr(1);
					read_line()
					ac = ((found = token.find(" AC ")) != std::string::npos);
					iLine.str(token.substr(0, found)); iLine.clear();

					iLine >> component_plus_node >> component_minus_node >> component_value;
					str_toupper(component_plus_node)
					str_toupper(component_minus_node)

					if(ac)
					{
						iLine.str(token.substr(found, std::string::npos));
						iLine >> token >> component_mag >> component_phase;
						ACComponentVector.push_back(ACComponent(component_type, component_name, component_value, component_mag, component_phase, nodeHash[component_plus_node], nodeHash[component_minus_node], branchID));
					}

					ComponentVector.push_back(std::unique_ptr <Component> (new Component(component_type, component_name, nodeHash[component_plus_node], nodeHash[component_minus_node], component_value)));
					if (!iLine.eof())
					{
						TransientComponentVector.push_back(TransientComponent(component_type, component_name, iLine, component_value, nodeHash[component_plus_node], nodeHash[component_minus_node], branchID));
					}
					if ((component_type=='V') || (component_type=='L'))
					{
  						ComponentVector[ComponentVector.size()-1]->set_branch(branchID++);
					}
					break;
				default: getline(iFile, token);
			}
			iFile >> token;
			str_toupper(token)
		} while (!iFile.eof());

		for (auto &it : ComponentVector)
		{
			i 				= it->plus();
			j				= it->minus();
			component_value = it->value();
			switch ( it->type() ) {
				case ('R'):
					g = (double) 1.0/component_value;
					MNA_DC(j, j) += g;
					MNA_DC(i, i) += g;
					MNA_DC(i, j) -= g;
					MNA_DC(j, i) -= g;
					break;
				case('V'):
					k = ( it->branch() + nodeID);
					MNA_b[k] 	 += component_value;
					MNA_DC(j, k) -= 1;
					MNA_DC(k, j) -= 1;
					MNA_DC(i, k) += 1;
					MNA_DC(k, i) += 1;
					break;
				case('I'):
					MNA_b[j] += component_value;
					MNA_b[i] -= component_value;
					break;
				case('L'):
					k = ( it->branch() + nodeID);
					MNA_DC(j, k) -= 1;
					MNA_DC(k, j) -= 1;
					MNA_DC(i, k) += 1;
					MNA_DC(k, i) += 1;
					MNA_TRAN(k, k) -= component_value;
					break;
				case('C'):
					MNA_TRAN(j, j) += component_value;
					MNA_TRAN(i, i) += component_value;
					MNA_TRAN(i, j) -= component_value;
					MNA_TRAN(j, i) -= component_value;
				default: break; 
			}
		}

		for(i=1; i < vectorSize; i++)
			for(j=1; j < vectorSize; j++)
				MNA_AC(i, j) = MNA_DC(i,j) + std::complex<double>(0.0, MNA_TRAN(i,j));

		for(auto &it : ACComponentVector)
		{
			if (it.type() == 'V')
				MNA_AC_b[it.branch() + MNA_DC.nodeID()] = it.value();
			else
			{
				MNA_AC_b[it.plus ()] -= it.value();
				MNA_AC_b[it.minus()] += it.value();
			}
		}

		for(auto &it : ACSweepAnalysisVector) 
		{
			it->update(nodeHash);
			it->analyse(MNA_AC, MNA_AC_b, itol, iter, spd, ACComponentVector);
		}

		for(auto &it : TransientAnalysisVector)
		{
			it->update(nodeHash);
			it->analyse(MNA_DC, MNA_TRAN, MNA_x, MNA_b, itol, iter, spd, method, TransientComponentVector);
		}

		if(iter)
		{
			spd ? MNA_DC.CG(MNA_x, MNA_b, itol) : MNA_DC.BiCG(MNA_x, MNA_b, itol);
		}
		else
		{
			spd ? MNA_DC.Cholesky(): MNA_DC.LU();
			spd ? MNA_DC.solveSPD(MNA_x, MNA_b) : MNA_DC.solve(MNA_x, MNA_b);
		}

		oFile.open("DC_Operating_Point_Analysis.txt");
		for(i=1 ; i < nodeID; ++i) oFile << nodeHash_reverse[i] << " " << MNA_x[i] << "\n";
		oFile.close();

		for(auto &it : DCSweepAnalysisVector) 
		{
			it->update(ComponentVector[componentHash[it->name()]]->plus(), ComponentVector[componentHash[it->name()]]->minus(), ComponentVector[componentHash[it->name()]]->branch(), ComponentVector[componentHash[it->name()]]->value(), nodeHash);
			it->analyse(MNA_DC, MNA_b, itol, iter, spd);
		}
		delete [] MNA_b, MNA_x;
	}
	return 0;
}