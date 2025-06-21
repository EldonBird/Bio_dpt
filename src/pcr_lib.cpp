#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace py = pybind11;

class snp_data {
    public:
        std::string rsid;
        std::string allele;
        std::string primer_sequence;
        std::string direction;
        int length;
		std::string gc_content;
		std::string hairpin;
		std::string homodimer;


        std::string get(const std::string& key) const {
   			if (key == "rsid") return rsid;
    		if (key == "allele") return allele;
    		if (key == "primer_sequence") return primer_sequence;
    		if (key == "direction") return direction;
    		if (key == "gc_content") return gc_content;
    		if (key == "hairpin") return hairpin;
    		if (key == "homodimer") return homodimer;
    		if (key == "length") return std::to_string(length);
            throw std::invalid_argument("Invalid key name");
		}
};

class snp_array{
	  public:
			std::vector<snp_data> data;

			snp_data get(const int index) const {
                          if( index < 0){ throw std::out_of_range("the index is out of raeng"); }

	                  return data[index];
			}
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
        .def(py::init<>())
        .def_readwrite("rsid", &snp_data::rsid)
        .def_readwrite("allele", &snp_data::allele)
        .def_readwrite("primer_sequence", &snp_data::primer_sequence)
        .def_readwrite("direction", &snp_data::direction)
        .def_readwrite("length", &snp_data::length)
        .def_readwrite("gc_content", &snp_data::gc_content)
        .def_readwrite("hairpin", &snp_data::hairpin)
        .def_readwrite("homodimer", &snp_data::homodimer)
    	.def("__getitem__", &snp_data::get);

    py::class_<snp_array>(m, "snp_array")
		.def(py::init<>())
    	.def_readwrite("data", &snp_array::data)
    	.def("__getitem__", &snp_array::get);


}
