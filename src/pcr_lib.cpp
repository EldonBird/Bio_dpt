#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace py = pybind11;

class primer {
    public:
        std::string rsid;
        std::string allele;
        std::string primer_sequence;
        std::string direction;
        int length;
		std::string gc_content;
        int tm;
		int hairpin;
		int homodimer;


        std::string get(const std::string& key) const {
   			if (key == "rsid") return rsid;
    		if (key == "allele") return allele;
    		if (key == "primer_sequence") return primer_sequence;
    		if (key == "direction") return direction;
    		if (key == "gc_content") return gc_content;
    		if (key == "tm") return std::to_string(tm);
    		if (key == "hairpin") return std::to_string(hairpin);
    		if (key == "homodimer") return std::to_string(homodimer);
    		if (key == "length") return std::to_string(length);
            throw std::invalid_argument("Invalid key name");
		}
};

class primer_group{
	  public:
			std::vector<primer> data;

			primer get(const int index) const {
                          if( index < 0){ throw std::out_of_range("the index is out of raeng"); }

	                  return data[index];
			}
			void add(primer p) {
				data.push_back(p);
			}
};




primer_group generate_allele_spesific_primers(primer data, int min_length, int max_length) {

	primer_group result;



	return result;
}

std::vector<std::string> introduce_missmatch(std::string str, int pos) {

	std::vector<std::string> primers;


	return primers;
}

std::unordered_map<std::string, std::string> evaluate_primer(const std::string & string) {

	std::unordered_map<std::string, std::string> result;




	return result;
}


primer_group filter_primers(primer_group evaluated_primers, int tm_min, int tm_max, double hairpin_max, double homodimer_max) {
  	primer_group result;

    for (int i = 0; i < evaluated_primers.data.size(); i++) {
      	primer current_primer = evaluated_primers.get(i);


    	if (current_primer.tm < tm_min) continue;
    	if (current_primer.tm > tm_max) continue;
    	if (current_primer.hairpin > hairpin_max) continue;
    	if (current_primer.homodimer > homodimer_max) continue;

		result.add(current_primer);
    }

	return result;
}

primer_group rank_primers(primer_group data){
	primer_group result;


	return result;
}

primer_group generate_matching_primers(primer_group data, primer_group allele_spesific_primers, int min_distance, int max_distance) {

	primer_group result;



    return result;
}

primer_group check_multiplex_compatibility(primer_group data, double heterodimer_max){
  primer_group result;



  return result;
}





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
    m.def("reverse_complement", &reverse_complement, "This function takes a string, and returns the reverse complement of the string or sequence.");
    m.def("generate_allele_spesific_primers", &generate_allele_spesific_primers, "This function takes in a an array of primers or snp data and a minimum and maximun length, and returns a generated list of allele spesific primers.");
    m.def("introduce_missmatch", &introduce_missmatch, "Takes in a sequence as a string, and a position, and returns a *list* of missmatches.");
    m.def("filter_primers", &filter_primers, "Takes in an array of primers, a tm min & max, a hairpin max, and a homodimer max, and returns the filtered list of primers.");
    m.def("rank_primers", &rank_primers, "Takes in an array of primers, and returns an ranked array of primers.");
    m.def("generate_matching_primers", &generate_matching_primers, "Takes in 2 arrays of primers, and a min & max distance, and returns an array of matching primers.");
    m.def("check_multiplex_compatibility", &check_multiplex_compatibility, "Takes in a list of primers, and a heterodimer maximun, and returns an array of all possible primer combinations? REDO THIS ONE!!!!!!!!!!!!!!!");

    py::class_<primer>(m, "primer")
        .def(py::init<>())
        .def_readwrite("rsid", &primer::rsid)
        .def_readwrite("allele", &primer::allele)
        .def_readwrite("primer_sequence", &primer::primer_sequence)
        .def_readwrite("direction", &primer::direction)
        .def_readwrite("length", &primer::length)
        .def_readwrite("gc_content", &primer::gc_content)
        .def_readwrite("hairpin", &primer::hairpin)
        .def_readwrite("homodimer", &primer::homodimer)
    	.def("__getitem__", &primer::get);

    py::class_<primer_group>(m, "primer_group")
		.def(py::init<>())
    	.def_readwrite("data", &primer_group::data)
    	.def("__getitem__", &primer_group::get);


}
