//LQ
#include <stdio.h>
#include <string>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <strstream>
#include <cmath>
#include "math.hpp"
#include "math-quat.hpp"
#include "readPDB.hpp"

#define THIRD 0.333333333333333
#define ROOTTHREE 1.73205080756888       
#define PI 3.141592653589793238462643383279502884197
 
using namespace std;

bool BothAreSpaces(char lhs, char rhs) { 
  // From the web
  return (lhs == rhs) && (lhs == ' '); 
}

int main(int argc, char *argv[])
{
  double const qtoC = 1.60e-19; // 1 elementary charge to Coulomb
  double const CoulCost = 8.99e9; // Coulomb's constant (N m^2)/C^2
  double const Avogadro = 6.022e23; // Avogadro number
  double const dielWat = 78.4; // Water dielectric constant at 298K 
  double const Jtokcal = 0.000239; // Joule to kcal
  double const JtoeV = 6.242e+18; // Joule to electronVolt
  double const AAtom = 1e-10; // Angstrom to m

  if (argc < 4 || argc > 5)
    {
      //Check number of parameters
      std::cout << "[ERROR] Wrong number of input parameters: ./ComputeCoulOnAtom <pdbfile> <psffile> <atomid> (<dielectric>)" << std::endl;
      return 0;
    }

  if (string(argv[1]) == "--help") 
    {
      std::cout << " " << std::endl;
      std::cout << " Welcome to the ComputeCoulOnAtom code!" << std::endl;
      std::cout << " This code will compute Coulomb interactions between one chosen atom and **all** other atoms of the sample." << std::endl;
      std::cout << " A single configuration (pdb) and its topology (psf) are needed. " << std::endl;
      std::cout << " The code assumes dielectric constant to be that of water at 298K (i.e. 78.4), if not otherwise specified." << std::endl;
      std::cout << " " << std::endl;
      std::cout << " ./ComputeCoulOnAtom <pdbfile> <psffile> <atomid> (<dielectric>)" << std::endl;
    }
  else 
    {
      ifstream infile (argv[1]) ;
      if ( !infile )
	{
	  std::cout <<"ERROR: Could not open configuration file " << argv[1] << std::endl;
	  return 0;
	}
      
      ifstream topfile (argv[2]) ;
      if ( !topfile )
	{
	  std::cout <<"ERROR: Could not open topology file " << argv[2] << std::endl;
	return 0;
      }

    // Read chosen atom id
    int atomsel = atoi(argv[3]);

    // Assign dielectric constant
    double dielecost;
    if ( argc == 5 )
      {
	dielecost = atof(argv[4]);
      }
    else 
      {
	//	std::cout << "Assuming dielectric constant of water at 298K (78.4)." << std::endl;
	dielecost = dielWat;
      }

    // Read PDB
    string str;
    vector<string> vec; 
    // Get PDB size
    while (getline(infile, str))
      {
	if (str.size() > 0)
	  {
	    vec.push_back(str);
	  }
      }

    infile.close();
    //std::cout << "Total number of lines: " << vec.size() << " " << std::endl;
  
    // Save PDB info
    double sidea, sideb, sidec;
    double alpha, beta,  gamma;
    string sgroup;
    double zvalue;

    vector<int> id;
    vector<string> atom;
    vector<string> altloc;
    vector<string> residue;
    vector<string> chain;
    vector<string> resseq;
    vector<string> icode;
    vector<double> coordx;
    vector<double> coordy;
    vector<double> coordz;
    vector<double> occup; //Occupancy
    vector<double> tempf; //TempFactor
    vector<string> element; 
    vector<string> charge; 

    string cryst = "CRYS";
    string at = "ATOM";
    string ha = "HETA";
    string cnc = "CONN";
    string end = "END";

    int totatom = 0;

    // Save selected atom info
    double selx, sely, selz;
    string atsel, ressel, chsel;


    for (int i=0; i<vec.size(); i++)
      {
	string start (vec[i], 0, 4);
	//std::cout << "First 4 chars: "  << start << " (" << i << ") " << std::endl;
	//std::cout << " " << std::endl;	

	if (start == cryst)
	  {
	    //std::cout << "header section" << std::endl;

	    ReadCryst(vec[i], sidea, sideb, sidec, alpha, beta, gamma, sgroup, zvalue);

	    //std::cout << "HEADER " <<  setprecision(6)  << setw(8) << start << " " << sidea << " - "  << sideb << " - " << sidec << " - " << alpha << " - " << beta << " - " << gamma << " - " << sgroup << " - "  << zvalue << std::endl;
	    //std::cout << " " << std::endl;	
	  }
	else if (start == at || start == ha || start == cnc)
	  {

	    //std::cout << "atom section" << std::endl;
	    
	    ReadAtom(vec[i], id, atom, altloc, residue, resseq, chain, icode, coordx, coordy, coordz, occup, tempf, element, charge);

	    if (id.back() == atomsel)
	      {
		selx = coordx.back();
		sely = coordy.back();
		selz = coordz.back();
		atsel = atom.back();
		ressel = residue.back();
		chsel = chain.back(); //chain
	      }
	    
	    //std::cout << "ATOM " << std::setprecision(6)  << std::setw(8) << " - " << id[0] << " - " << atom[0] << " - " << residue[0]  << " - " << resseq[0] << " - " << coordx[0] << " - " <<  coordy[0] << " - " << coordz[0] << " - " << occup[0] << " - " << tempf[0] << " - " << element[0] << " - " <<  charge[0] << " - " << std::endl;

	    totatom++;
	  }
	else if (start == end)
	  {
	    //ADD???
	  }
	else
	  {
	    std::cout << "error at line " << i << std::endl;
	    return 0;
	  }

      }

    // Print selected atom info
    //    std::cout << "Selected info: id  " << atomsel << ", type " << atsel << ", residue " <<ressel << ", at "  << selx << " " << sely << " " << selz << ", file " << argv[1] << std::endl;
    //    std::cout << " " << std::endl;

    // Reading PSF (topology)
    //    std::cout << "Reading topology from file " << argv[2] << std::endl;
    //    std::cout << " " << std::endl;
    string natom = "NATOM"; // Keyword! Look for it in the psf to start saving charges

    vector<int> idpsf;
    vector<string> respsf;
    vector<string> atnamepsf;
    vector<string> attypepsf;
    vector<double> chargepsf;
    vector<double> masspsf;

    while ( !topfile.eof() )
      {
	string check;
        getline(topfile, check);
	// Start reading PSF from the line containing "NATOM"
	if ( check.find(natom) != std::string::npos )
	  {

	    while ( !topfile.eof() )
	      {
		string ch;
		getline(topfile, ch);
		
		if (ch.find("NBOND") != std::string::npos || ch.empty() )
		  {
		    // Stop reading PSF file at the end of "atom" section
		    break;
		  }
		else
		  {
		    string null;
		    int idt; //temporary storage 
		    double charget, masst; //temporary storage

		    // Remove double spaces
		    std::string::iterator new_end = std::unique(ch.begin(), ch.end(), BothAreSpaces);
		    ch.erase(new_end, ch.end());   

		    //Tokenizing a string
		    // from: https://www.geeksforgeeks.org/tokenizing-a-string-cpp/
		    // stringstream class check1 
		    stringstream check1(ch); 
		    string intermediate; 
      		    vector<string> tokens;
		    // Tokenizing w.r.t. space ' ' 
		    while(getline(check1, intermediate, ' ')) 
		      { 
			tokens.push_back(intermediate); 
		      } 
      
		    // Printing the token vector 
		    //for(int i = 0; i < tokens.size(); i++) 
		    //  std::cout << tokens[i] << std::endl; 
		    //		    std::cout << "size tokens " << tokens.size() << " mao " << tokens[1] << std::endl;

		    //Fill vectors of PSF info
		    idt = atoi(tokens[1].c_str());
		    idpsf.push_back(idt);
		    
		    respsf.push_back(tokens[4].c_str());
		    atnamepsf.push_back(tokens[5].c_str());
		    attypepsf.push_back(tokens[6].c_str());

		    charget = atof(tokens[7].c_str());
		    chargepsf.push_back(charget);

		    masst = atof(tokens[8].c_str());
		    masspsf.push_back(masst);

		  }
	      }
	  }

      }

    topfile.close();
 
    // Check that both psf and pdb vectors have the same dimension
    //std::cout << idpsf.size() << " vs " << id.size() << " i.e. " << totatom << std::endl;

    // Assign correct charge to selected atom
    double qsel; //charge 
    for (int i=0; i<totatom; i++)
      {
	if (idpsf[i] == atomsel)
	    qsel = chargepsf[i];
      }

    // Print info of selected atom
    std::cout << "Atom selected: id  " << atomsel << ", type " << atsel << ", residue " <<ressel << ", position "  << selx << " " << sely << " " << selz << ", charge " << qsel << ", from file " << argv[1] << std::endl;
    std::cout << "Dielectric constant: " << dielecost << std::endl;
    std::cout << " " << std::endl;
 
    // Compute Coulomb energy
    double uelec = 0.0;
    for (int i=1; i<totatom+1; i++)
      {
	if (idpsf[i] != atomsel)
	  {
	    double q2 = chargepsf[i];
	    double crd2x, crd2y, crd2z;

	    // Check if the order of atoms in pdb and psf are the same
	    if (idpsf[i] == id[i])
	      {
		crd2x = coordx[i];
		crd2y = coordy[i];
		crd2z = coordz[i];
	      }
	    else
	      {
		for (int i=1; i<totatom+1; i++)
		  {
		    if (idpsf[i] == id[i])
		      {
			crd2x = coordx[i];
			crd2y = coordy[i];
			crd2z = coordz[i];
		      }
		  }
	      }

	    double distx = selx-crd2x;
	    double disty = sely-crd2y;
	    double distz = selz-crd2z;
	    double dist = distx*distx + disty*disty + distz*distz;
	    dist = sqrt(dist);
	    
	    //	    std::cout << "charges " << qsel <<  " and " << q2 << " at distance " << dist << " uelec-pre "  << uelec << " id-atom1 " << atomsel << " id-atom2 " << id[i] << std::endl;
	    uelec = uelec + (qsel*q2)/dist;
	  }
      }

    uelec = uelec*(qtoC*qtoC)/AAtom; // Convert to N m
    uelec = uelec*CoulCost; // Multiply for Coulomb's constant
    uelec = uelec/dielecost; // Divide for dielectric constant

    // Print electrostatic constant with different units
    std::cout << "The electrostatic energy of your selected atom interacting with the whole sample is: " << std::endl;
    std::cout << std::setw(15) << uelec << " J (N m) " << std::endl;
    std::cout << std::setw(15) << uelec*JtoeV << " eV "  << std::endl;
    std::cout << std::setw(15) << uelec*Jtokcal*Avogadro << " kcal/mol " << std::endl;


  }  
}

