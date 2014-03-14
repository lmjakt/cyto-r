#include "fcs3.h"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <sstream>

FCS3::FCS3(const char* fname)
{
  // set pointers to null
  events = 0;
  events_size = 0;
  data_start = 0;
  data_end = 0;
  byte_big_endian = true; // ? our machines seem to do this...

  std::ifstream in(fname);
  if(!in.is_open()){
    std::cerr << "FCS3 Constructor unable to open fcs 3 file" << std::endl;
    return;
  }
  if(!read_header(in)){
    std::cerr << "FCS3 Constructor error in reading header information" << std:: endl;
    return;
  }
  if(!read_text(in)){
    std::cerr << "FCS3 Constructor unable to read text information" << std::endl;
    return;
  }
  if(!read_data(in))
    std::cerr << "FCS3 Constructor unable to read data" << std::endl;
  in.close();
}

FCS3::~FCS3()
{
  for(unsigned int i=0; i < events_size; ++i)
    delete []events[i];
  delete []events;
}

bool FCS3::read_header(std::ifstream& in)
{
  char id[7];
  id[6] = 0;
  size_t t_s = 8;     // the size of the field indicating file pos
  size_t t_b = 10;
  size_t t_e = 18;
  char txt_beg[8];
  char txt_end[8];    // read only these.
  
  in.read(id, 6);
  if(in.fail() || strncmp(id, "FCS3.0", 6)){
    std::cerr << "read_header, incorrect header info: " << id << std::endl;
    return(false);
  }
  in.seekg(t_b);
  in.read(txt_beg, t_s);
  if(in.fail()){
    std::cerr << "read_header, unable to read txt_beg" << std::endl;
    return(false);
  }
  in.seekg(t_e);
  in.read(txt_end, t_s);
  if(in.fail()){
    std::cerr << "read_header, unable to read txt_end" << std::endl;
    return(false);
  }
  ssize_t td_beg = atol(txt_beg);
  ssize_t td_end = atol(txt_end);

  // some sanity checking
  if(!td_beg || !td_end || td_end <= td_beg){
    std::cerr << "Values obtained for beginning and end of text don't make sense: "
	      << td_beg << " -> " << td_end << std::endl;
    return(false);
  }
  text_begin = td_beg;
  text_end = td_end;
  return(true);
}

bool FCS3::read_text(std::ifstream& in)
{
  in.seekg(text_begin);
  if(in.fail()){
    std::cerr << "Unable to move to beginning of text " << text_begin << std::endl;
    return(false);
  }
  char sep;
  in.read(&sep, 1);
  if(in.fail())
    return(false);
  std::vector<std::string> kv_pair = read_key_value_pair(in, sep);
  while(kv_pair.size() == 2){
    key_values.insert(make_pair(kv_pair[0], kv_pair[1]));
    kv_pair = read_key_value_pair(in, sep);
  }
  std::cout << "Parsing pairs: " << std::endl;
  if(!parse_pairs()){
    std::cerr << "Unable to obtain all necessary information" << std::endl;
    return(false);
  }
  std::cout << "events : " << events_size << " parno: " << parameters.size() << std::endl;
  for(unsigned int i=0; i < parameters.size(); ++i){
    std::cout << i << "\t" << parameters[i] << "\t" << parameter_names[i] 
	 << "\t" << par_bit_sizes[i] << "\t" << par_ranges[i] << std::endl;
  }
  return(true);
}

float* FCS3::event(unsigned int i)
{
  if(i >= events_size)
    return(0);
  return(events[i]);
}

unsigned int FCS3::event_no()
{
  return(events_size);
}

unsigned int FCS3::par_no()
{
  return(parameters.size());
}

std::vector<std::string> FCS3::pars()
{
  return(parameters);
}

std::vector<std::string> FCS3::par_names()
{
  return(parameter_names);
}

std::string FCS3::par(unsigned int i)
{
  std::string blank;
  if(i >= parameters.size())
    return(blank);
  return(parameters[i]);
}

std::string FCS3::par_name(unsigned int i)
{
  std::string blank;
  if(i >= parameter_names.size())
    return(blank);
  return(parameter_names[i]);
}


void FCS3::write_text(const char* ofile)
{
  std::ofstream out(ofile);
  if(!out.is_open()){
    std::cerr << "Unable to open file for writing" << std::endl;
    return;
  }
  out << "index";
  for(unsigned int i=0; i < parameters.size(); ++i){
    if(parameters[i].length()){
      out << "\t" << parameters[i];
    }else{
      out << "\t" << parameter_names[i];
    }
  }
  out << std::endl;
  for(unsigned int i=0; i < events_size; ++i){
    out << i + 1;
    for(unsigned int j=0; j < parameters.size(); ++j)
      out << "\t" << events[i][j];
    out << std::endl;
  }
  out.close();
}

// can only read data from float fields.. 
bool FCS3::read_data(std::ifstream& in)
{
  in.seekg(data_start);
  if(in.fail()){
    std::cerr << "read_data unable to seek to beginning of data : " << data_start << std::endl;
    return(false);
  }
  events = new float*[events_size];
  for(unsigned int i=0; i < events_size; ++i){
    events[i] = new float[ parameters.size() ];
    in.read((char*)events[i], sizeof(float) * parameters.size());
    if(byte_big_endian)   /// bad, NON-PORTABLE CODE assumes small endian byte order for float
      reorder_bytes(events[i], parameters.size());
  }
  if(in.fail())
    return(false);
  return(true);
}

void FCS3::reorder_bytes(float* f, unsigned int l){
  float r;
  for(unsigned int i=0; i < l; ++i){
    char* s = (char*)(f + i);
    char* d = (char*)&r;
    d[0] = s[3];
    d[1] = s[2];
    d[2] = s[1];
    d[3] = s[0];
    f[i] = r;
  }
}

// this is kind of a slow way, but who cares..?
std::vector<std::string> FCS3::read_key_value_pair(std::ifstream& in, char sep)
{
  std::vector<std::string> value_pair;
  value_pair.reserve(2);
  char c;
  std::string key, value;
  while( in.tellg() <= (ssize_t)text_end ){
    in.read(&c, 1);
    if(c == sep)
      break;
    key.append(1, c);
  }
  while( in.tellg() <= (ssize_t)text_end ){
    in.read(&c, 1);
    if(c == sep)
      break;
    value.append(1, c);
  }
  if(in.fail() || in.tellg() >= (ssize_t)text_end)
    return(value_pair); // empty;
  boost::to_upper(key);
  boost::to_upper(value);  // should be case insensitive.. hence.. 
  value_pair.push_back(key);
  value_pair.push_back(value);
  return(value_pair);
}

bool FCS3::parse_pairs(){
  std::map<int, std::string> par_names;
  std::map<int, std::string> par_snames;

  std::string data_type;
  std::string byte_order; // not sure how to use this, but we need somehow  
  unsigned int par_no=0;

  // this is a stupid way of doing it. Just ask directly using
  // map.count(). Not sure why I ended up doing it like this.
  // kind of strange in the mind.
  for(std::map<std::string, std::string>::iterator it=key_values.begin();
      it != key_values.end(); ++it){
    if(it->first == "$TOT")
      events_size = atoi(it->second.c_str());
    if(it->first == "$PAR")
      par_no = atoi(it->second.c_str());
    if(it->first == "$DATATYPE")
      data_type = it->second;
    if(it->first == "$BEGINDATA")
      data_start = atol(it->second.c_str());
    if(it->first == "$DATATYPE")
      data_type = it->second;
    if(it->first == "$BYTEORD")
      byte_order = it->second;
    //std::cout << (*it).first << " -> " << (*it).second << std::endl;
  }
  if(!events_size || !par_no || !data_start){
    std::cerr << "parse_pairs unable to find necessary pars: "
	      << events_size << "\t" << par_no << "\t" << data_start << std::endl;
    return(false);
  }
  if(data_type != "F"){
    std::cerr << "parse_pairs: currently we only know how to handle floats" << std::endl;
    return(false);
  }
  if(byte_order != "4,3,2,1" && byte_order != "1,2,3,4"){
    std::cerr << "Unable to interpret byte order from " << byte_order << std::endl;
    return(false);
  }
  byte_big_endian = (byte_order == "4,3,2,1");
  parameters.resize(par_no);
  parameter_names.resize(par_no);
  par_bit_sizes.resize(par_no);
  par_ranges.resize(par_no);

  for(unsigned int i=0; i < parameters.size(); ++i){
    std::ostringstream par_b;
    std::ostringstream par_r;  // the first two are mandatory
    std::ostringstream par_s;
    std::ostringstream par_n;
    
    unsigned int j = i + 1;
    par_b << "$P" << j << "B";
    par_r << "$P" << j << "R";
    par_s << "$P" << j << "S";
    par_n << "$P" << j << "N";

    if(!key_values.count(par_b.str()) || !key_values.count(par_r.str())){
      std::cerr << "Unable to find paramter : " << par_b.str()
		<< " or " << par_r.str() << std::endl;
      return(false);
    }
    par_bit_sizes[i] = atoi(key_values[par_b.str()].c_str());
    par_ranges[i] = atof(key_values[par_r.str()].c_str());
    
    if(key_values.count(par_s.str()))
      parameters[i] = key_values[ par_s.str() ];
    if(key_values.count(par_n.str()))
      parameter_names[i] = key_values[ par_n.str() ];
  }
	
  return(true);
}
