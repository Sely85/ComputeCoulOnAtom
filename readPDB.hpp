#include <string>
#include <vector>
#include <sstream>
#include <iostream>

void ReadCryst(std::string str, double &sidea, double &sideb, double &sidec, double &alpha, double &beta,  double &gamma, std::string &sgroup, double &zvalue)
{

  sidea   = atof(str.substr(6,10).c_str());  
  sideb   = atof(str.substr(15,10).c_str());
  sidec   = atof(str.substr(24,10).c_str());
  alpha   = atof(str.substr(33,8).c_str());  
  beta    = atof(str.substr(40,8).c_str()); 
  gamma   = atof(str.substr(47,8).c_str()); 
  sgroup   = (str.substr(55,12).c_str()); 
  zvalue  = atof(str.substr(66,5).c_str()); 

  // std::cout << std::setprecision(6)  << std::setw(8) << " " << sidea << " - "  << sideb << " - " << sidec << " - " << alpha << " - " << beta << " - " << gamma << " - " << sgroup << " - "  << zvalue << std::endl;
  
}



void ReadAtom(std::string str, std::vector<int> &id, std::vector<std::string> &atom, std::vector<std::string> &altloc, std::vector<std::string> &residue, std::vector<std::string> &resseq, std::vector<std::string> &chain, std::vector<std::string> &icode, std::vector<double> &coordx, std::vector<double> &coordy, std::vector<double> &coordz, std::vector<double> &occup, std::vector<double> &tempf, std::vector<std::string> &element, std::vector<std::string> &charge)
{

  id.push_back(atoi(str.substr(6,6).c_str()));  
  atom.push_back(str.substr(12,5).c_str());  
  altloc.push_back(str.substr(16,2).c_str());  
  residue.push_back(str.substr(17,4).c_str());  
  chain.push_back(str.substr(21,2).c_str());  
  resseq.push_back(str.substr(22,5).c_str());  
  icode.push_back(str.substr(26,2).c_str());  
  coordx.push_back(atof(str.substr(30,9).c_str()));  
  coordy.push_back(atof(str.substr(38,9).c_str()));  
  coordz.push_back(atof(str.substr(46,9).c_str()));  
  occup.push_back(atof(str.substr(54,7).c_str()));  
  tempf.push_back(atof(str.substr(60,7).c_str()));  
  element.push_back(str.substr(76,3).c_str());  
  charge.push_back(str.substr(78,3).c_str());  

}
