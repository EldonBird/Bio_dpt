#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <bits/fs_ops.h>
#include <pstl/execution_defs.h>

namespace py = pybind11;

class primer {
    public:
		// ?
        std::string rsid;
		std::string snp_id;

		std::string allele;
		std::string primer_sequence;
		std::string direction;

		float score;
		float gc;
		float tm;
		float hairpin;
		float homodimer;
		int position;
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

std::string introduce_mismatch(std::string str) {

	std::vector<std::string> result;

	std::unordered_map<char, char> complement = {
      
		{'A', 'T'}, {'T', 'A'},
		{'C', 'G'}, {'G', 'C'},
		{'a', 't'}, {'t', 'a'},
		{'c', 'g'}, {'g', 'c'}
    
	};

	if (str.length() < 3) {
		throw std::invalid_argument("Invalid Sequence, Not larger than 3.");
	}

	int position = str.length() - 3;
	char base = str[position];
	char missmatch = complement[base];

	if (missmatch == ' ') {
		throw std::invalid_argument("Warning: No mismatch rule for base" + str);
	}
	
	std::string new_str = "";

	for (int i = 0; i < str.length(); i++) {
		
		if (i == position) {
			new_str += missmatch;
		}
		else {
			str[i];
		}
	}
	return new_str;	
}

std::vector<primer> generate_allele_specific_primers(std::vector<primer> data, int min_length, int max_length) {
	std::vector<primer> result;

	for (int i = 0; i < data.data.size(); i++) {

		std::string snp_id = data[i].snp_id;
		std::string allele = data[i].allele;
		std::string sequence = data[i].primer_sequence;
		int center = data[i].position;

		std::string forward = sequence.substr(center - max_length + 1, center + 1);
		std::string forward_mismatch = introduce_mismatch(forward);

		if (forward.length() >= min_length) {

			for (int length = 0; length < max_length - min_length - 1; length++) {

				std::string trimmed = forward_mismatch.substr(length);
				primer p;

				// might have a constructor later, please keep this in consideration

				p.snp_id = snp_id;
				p.allele = allele;
				p.primer_sequence = trimmed;
				p.direction = "forward";
				p.length = forward.length() - length;

				result.push_back(p);
			}
		}
		else {
			throw std::invalid_argument("Error: thte forward length was not greater than or equal to the min length for item, " +
				std::to_string(i) + " : " + snp_id + " . This error happened in the generate_allele_specific_primers");
		}


		// I have no idea if this is what was intended, please have someone double check this ;)
		std::string reverse = reverse_complement(sequence.substr(center + 1, center + max_length - 1));
		std::string reverse_mismatch = introduce_mismatch(reverse);

		if (reverse.length() >= min_length) {

			for (int length = 0; length < max_length - min_length - 1; length++) {

				// double check that this is doing what I think that it is doing
				std::string trimmed = reverse_mismatch.substr(length);
				primer p;

				p.snp_id = snp_id;
				p.allele = allele;
				p.primer_sequence = trimmed;
				p.direction = "reverse";
				p.length = reverse.length() - length;

				result.push_back(p);
			}
		}
		else {
			throw std::invalid_argument("Error: the reverse length was not greater than or equal to the min length for item, " +
				std::to_string(i) + " : " + snp_id + " . This error happened in the generate_allele_specific_primers");
		}
	}

	
	return result;
}



std::unordered_map<std::string, std::string> evaluate_primer(const std::string seq) {

	std::unordered_map<std::string, std::string> result;




	return result;
}

std::vector<primer> filter_primers(std::vector<primer> evaluated_primers, int tm_min, int tm_max, double hairpin_max, double homodimer_max) {
  	std::vector<primer> result;

    for (int i = 0; i < evaluated_primers.data.size(); i++) {
      	primer current_primer = evaluated_primers[i];


    	if (current_primer.tm < tm_min) continue;
    	if (current_primer.tm > tm_max) continue;
    	if (current_primer.hairpin > hairpin_max) continue;
    	if (current_primer.homodimer > homodimer_max) continue;

		result.push_back(current_primer);
    }

	return result;
}

std::vector<primer> rank_primers(std::vector<primer> data, float target_tm, float target_gc){

	for (int i = 0; i < data.size(); i++) {

		float tm_score = abs(data[i].tm - target_tm);
		float gc_Score = abs(data[i].gc - target_gc);
		data[i].score = tm_score + gc_Score + data[i].hairpin + data[i].homodimer;
		
	}

	// I think that this is sorted correctly...?

	std::sort(data.begin(), data.end(), [](primer a, primer b) {
		return a.score > b.score;
	});
	return data;
}

std::vector<primer> generate_matching_primers(std::vector<primer> data, std::vector<primer> allele_spesific_primers, int min_distance, int max_distance) {

	std::vector<primer> result;

	for (int i = 0; i < allele_spesific_primers.size(); i++) {
		
		int snp_id = allele_spesific_primers[i].snp_id;
		std::string allele = allele_spesific_primers[i].allele;
		std::string direction = allele_spesific_primers[i].direction;
		std::string sequence = allele_spesific_primers[i].primer_sequence;
		int center = allele_spesific_primers[i].center;



		for (int dist = min_distance; dist <= max_distance + 1; dist += 10) {

			for (int length = 18; length < 29; length++) {

				if (direction == "forward") {
					int start = center + dist;
					std::string intermediate = sequence.substr(start, start + length);
					std::string primer_seq = reverse_complement(intermediate);
				}
				else {
					int start = center - dist - length;
					std::string primer_Sec = sequence.substr(start, start+length);
				}
				std::unordered_map<std::string, std::string> metrics = evaluate_primer(sequence);

				if (metric.tm >= 60.0 && metric.tm < 65.0 &&
					metrics.gc >= 40.0 && metrics.gc <= 60.0 &&
					metrics.hairpin < 45.0 &&
					metrics.homodimer < 45.0) {

					result.data.push_back(allele_spesific_primers.get(i));

					
				}

				
			}
			
		}

		
	}



    return result;
}

std::vector<primer> check_multiplex_compatibility(std::vector<primer> data, double heterodimer_max){
  std::vector<primer> result;



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
        .def_readwrite("homodimer", &primer::homodimer);

}
