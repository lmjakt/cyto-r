#ifndef FCS3_H
#define FCS3_H

#include <vector>
#include <string>
#include <map>
#include <fstream>

// an object representing the data in an FCS file.
// keeps everything as simple as possible
// Keeps all values as floats.

class FCS3 {
 public :
  FCS3(const char* fname);
  ~FCS3();

  float* event(unsigned int i);
  unsigned int event_no();
  unsigned int par_no();
  std::vector<std::string> pars();
  std::vector<std::string> par_names();
  std::string par(unsigned int i);
  std::string par_name(unsigned int i);

  void write_text(const char* ofile);

 private:
  float** events;
  unsigned int events_size;
  std::vector<std::string> parameters;
  std::vector<std::string> parameter_names;
  std::vector<unsigned int> par_bit_sizes;
  std::vector<float> par_ranges;
  
  bool byte_big_endian;  // byte order
  
  std::map<std::string, std::string> key_values;
  
  size_t text_begin;
  size_t text_end;
  size_t data_start;
  size_t data_end;

  bool read_header(std::ifstream& in);
  bool read_text(std::ifstream& in);
  bool read_data(std::ifstream& in);
  void reorder_bytes(float* f, unsigned int l);

  std::vector<std::string> read_key_value_pair(std::ifstream& in, char sep);
  bool parse_pairs();
};

#endif
