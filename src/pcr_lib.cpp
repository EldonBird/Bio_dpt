#include <algorithm>
#include <ranges>
#include <filesystem>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <pybind11/stl.h>
#include <bits/fs_ops.h> // enable these if on linux...
#include <pstl/execution_defs.h> // other one
// #include <filesystem> // and this for windows...


namespace fs = std::filesystem;
namespace py = pybind11;

class primer {
    public:
		// ?
        std::string rsid;
		std::string snp_id;

		std::string allele;
		std::string sequence;
		std::string direction;

		float score;
		float gc;
		float tm;
		float hairpin;
		float homodimer;
		int position;
        int length;
		int best; 
};

std::vector<primer> df_to_listprimers(const py::object& df) {
	std::vector<primer> result;

	py::list ids = df.attr("__getitem__")("snp_id");
	py::list allele = df.attr("__getitem__")("allele");
	py::list sequence = df.attr("__getitem__")("sequence");
	py::list positions = df.attr("__getitem__")("position");

	auto length = py::len(df);

	for (size_t i = 0; i < length; i++) {

		primer p;
		p.snp_id = ids[i].cast<std::string>();
		p.allele = allele[i].cast<std::string>();
		p.sequence = sequence[i].cast<std::string>();
		p.position = positions[i].cast<int>();
		result.push_back(p);
		
	}

	
	return result;
}

py::object& df listprimers_to_df(vector<primer> data) {




	
}


std::string reverse_complement(std::string s) {
  
	std::unordered_map<char, char> complement = {
      
		{'A', 'T'}, {'T', 'A'},
		{'C', 'G'}, {'G', 'C'},
		{'a', 't'}, {'t', 'a'},
		{'c', 'g'}, {'g', 'c'}
    
	};
  
	std::string result;
    
	for(std::size_t i = 0; i < s.length(); i++) {
      
		char cur_complement = complement[s[i]];
		result += cur_complement;
      
	}
    
    
	return result;
}

std::string introduce_mismatch(std::string primer_sequence) {

	std::unordered_map<char, char> complement = {
      
		{'A', 'T'}, {'T', 'A'},
		{'C', 'G'}, {'G', 'C'},
		{'a', 't'}, {'t', 'a'},
		{'c', 'g'}, {'g', 'c'}
    
	};

	if (primer_sequence.length() < 3) {
		throw std::invalid_argument("Invalid Sequence, Not larger than 3 : " + primer_sequence);
	}

	size_t position = primer_sequence.length() - 3;
	char base = primer_sequence[position];
	char mismatch = complement[base];

	if (mismatch == ' ') {
		throw std::invalid_argument("Warning: No mismatch rule for base" + primer_sequence + " : " + mismatch);
	}

	// this might need an improvement, I don't know if the substring is inclusive or if the position replacement is correct
	return primer_sequence.substr(0, position) + mismatch + primer_sequence.substr(position + 1, primer_sequence.length());
	
}

std::vector<primer> generate_allele_specific_primers(std::vector<primer> data, std::size_t min_length, int max_length) {
	std::vector<primer> result;

	for (std::size_t i = 0; i < data.size(); i++) {

		std::string snp_id = data[i].snp_id;
		std::string allele = data[i].allele;
		std::string sequence = data[i].sequence;
		int center = data[i].position;

		// took out the + 1 from center - max_length + 1 double check this
		std::string forward = sequence.substr(center - max_length, center + 1);
		std::string forward_mismatch = introduce_mismatch(forward);

		if (forward.length() >= min_length) {

			for (std::size_t length = 0; length < max_length - min_length + 1; length++) {

				std::string trimmed = forward_mismatch.substr(length);
				primer p;

				// might have a constructor later, please keep this in consideration

				p.snp_id = snp_id;
				p.allele = allele;
				p.sequence = trimmed;
				p.direction = "forward";
				p.length = forward.length() - length;

				result.push_back(p);
			}
		}
		else {
			throw std::invalid_argument("Error: the forward length was not greater than or equal to the min length for item, " +
				std::to_string(i) + " : " + snp_id + " . This error happened in the generate_allele_specific_primers");
		}


		// I have no idea if this is what was intended, please have someone double check this ;)
		std::string reverse = reverse_complement(sequence.substr(center, center + max_length + 1));
		std::string reverse_mismatch = introduce_mismatch(reverse);

		if (reverse.length() >= min_length) {

			for (std::size_t length = 0; length < max_length - min_length + 1; length++) {

				// double check that this is doing what I think that it is doing
				std::string trimmed = reverse_mismatch.substr(length);
				primer p;

				p.snp_id = snp_id;
				p.allele = allele;
				p.sequence = trimmed;
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

std::vector<primer> filter_primers(std::vector<primer> evaluated_primers, int tm_min, int tm_max, double hairpin_max, double homodimer_max) {
	std::vector<primer> result;

	for (std::size_t i = 0; i < evaluated_primers.size(); i++) {
		primer current_primer = evaluated_primers[i];


		if (current_primer.tm <= tm_min) continue;
		if (current_primer.tm > tm_max) continue;
		if (current_primer.hairpin >= hairpin_max) continue;
		if (current_primer.homodimer >= homodimer_max) continue;

		result.push_back(current_primer);
	}

	return result;
}

std::vector<primer> rank_primers(std::vector<primer> data, float target_tm, float target_gc){

	for (std::size_t i = 0; i < data.size(); i++) {

		float tm_score = abs(data[i].tm - target_tm);
		float gc_Score = abs(data[i].gc - target_gc);
		data[i].score = tm_score + gc_Score + data[i].hairpin + data[i].homodimer + data[i].best; // maybe include length of sequences here in the score
		
	}

	// I think that this is sorted correctly...?

	std::sort(data.begin(), data.end(), [](primer a, primer b) {
		return a.score > b.score;
	});
	return data;
}

std::unordered_map<std::string, float> evaluate_primer(const std::string seq) {

	std::unordered_map<std::string, float> result;

	
	


	return result;
}



std::vector<primer> generate_matching_primers(std::vector<primer> snp_data, std::vector<primer> allele_specific_primers, std::size_t min_distance, std::size_t max_distance) {

	std::vector<primer> matching_primers;

	for (std::size_t i = 0; i < allele_specific_primers.size(); i++) {

		std::string snp_id = allele_specific_primers[i].snp_id;
		std::string allele = allele_specific_primers[i].allele;
		std::string direction = allele_specific_primers[i].direction;

		std::string sequence = snp_data[i].sequence;
		int center = snp_data[i].position;

		for (std::size_t dist = min_distance; dist <= max_distance + 1; dist += 10) {

			for (std::size_t len = 18; len < 29; len++) {

				int start;
				std::string primer_seq;

				if (direction == "forward") {
					start = center + dist;
					primer_seq = reverse_complement(sequence.substr(start, start + len));
				}
				else {
					start = center - dist - len;
					primer_seq = reverse_complement(sequence.substr(start, start + len));
				}

				// this doesn't currently actually do anything, I am not quite sure if I should work on this further I have to find out what the,
				// primer 3 call is doing, this might take too long, or in fact, be impossible. Got rid of the evaluate primers, we decided that we didn't
				// really understand what it was doing in comparison to filter/rank. We can re add it when needed when the logic is understood.
				
				//std::unordered_map<std::string, float> metrics = evaluate_primer(primer_seq);


				// if (metrics["tm"] >= 60.0 and metrics["tm"] <= 65.0 and
				// 	metrics["gc"] >= 40.0 and metrics["gc"] <= 60.0 and
				// 	metrics["hairpin"] < 45.0 and
				// 	metrics["homodimer"] < 45.0 ) {
					
					primer p;
					p.snp_id = snp_id;
					p.allele = allele;
					p.sequence = sequence;
					p.direction = direction;
					

					
				//}
			}
		}
	}
	return matching_primers;
}

std::vector<primer> check_multiplex_compatibility(std::vector<primer> data, double heterodimer_max){
	std::vector<primer> result;

	
	
	

	return result;
}





PYBIND11_MODULE(pcr_lib, m) {
    m.def("reverse_complement", &reverse_complement, "This function takes a string, and returns the reverse complement of the string or sequence.");
    m.def("generate_allele_specific_primers", &generate_allele_specific_primers, "This function takes in a an array of primers or snp data and a minimum and maximun length, and returns a generated list of allele specific primers.");
    m.def("introduce_mismatch", &introduce_mismatch, "Takes in a sequence as a string, and a position, and returns a *list* of missmatches.");
    m.def("filter_primers", &filter_primers, "Takes in an array of primers, a tm min & max, a hairpin max, and a homodimer max, and returns the filtered list of primers.");
    m.def("rank_primers", &rank_primers, "Takes in an array of primers, and returns an ranked array of primers.");
    m.def("generate_matching_primers", &generate_matching_primers, "Takes in 2 arrays of primers, and a min & max distance, and returns an array of matching primers.");
    m.def("check_multiplex_compatibility", &check_multiplex_compatibility, "Takes in a list of primers, and a heterodimer maximun, and returns an array of all possible primer combinations? REDO THIS ONE!!!!!!!!!!!!!!!");
    m.def("df_to_listprimers", &df_to_listprimers, "df to list of primers :)");
    m.def("listprimers_to_df", &listprimers_to_df, "list of primers to pandas df ;)");

    py::class_<primer>(m, "primer")
        .def(py::init<>())
        .def_readwrite("rsid", &primer::rsid)
        .def_readwrite("snp_id", &primer::snp_id)
        .def_readwrite("allele", &primer::allele)
        .def_readwrite("sequence", &primer::sequence)
        .def_readwrite("direction", &primer::direction)
        .def_readwrite("score", &primer::score)
        .def_readwrite("gc", &primer::gc)
        .def_readwrite("tm", &primer::tm)
        .def_readwrite("hairpin", &primer::hairpin)
        .def_readwrite("homodimer", &primer::homodimer)
        .def_readwrite("position", &primer::position)
        .def_readwrite("length", &primer::length);
}
