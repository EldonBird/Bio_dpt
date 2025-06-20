#include <pybind11/pybind11.h>
#include <string>
#include <unordered_map>


    

std::string reverse_complement(std::string s) {
  
  std::unordered_map<char, char> complement = {
      
    {'A', 'T'}, {'T', 'A'},
    {'C', 'G'}, {'G', 'C'},
    {'a', 't'}, {'t', 'a'},
    {'c', 'g'}, {'g', 'c'}
    
    };
  
    std::string result;
    
    for(int i = 0; i < s.length(); i++) {
      
      char cur_complement = complement[s[i]];
      result += cur_complement;
      
    }
    
    
    return result;
    
  
}



PYBIND11_MODULE(pcr_lib, m) {
    m.def("reverse_complement", &reverse_complement, "This function takes a string, and returns the reverse complement of the string or sequence");
}
