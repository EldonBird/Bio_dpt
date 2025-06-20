#include <pybind11/pybind11.h>
#include <stdio.h>
#include <string>
#include <unordered_map>


class snp_data {
    public:
        std::string rsid;
        std::string allele;
        std::string primer_sequence;
        std::string direction;
        int length;
};



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
    
    py::class_<snp_data>(m, "snp_data")
        .def(py::init<const std::string&>())
        .def("rsid", &snp_data::rsid);


}
