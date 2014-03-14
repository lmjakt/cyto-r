#include "fcs3.h"
#include <iostream>
#include <stdlib.h>

int main(int argc, char** argv){
  if(argc < 2){
    std::cerr << "read_fcs filne_name out.txt" << std::endl;
    exit(1);
  }
  FCS3 fcs(argv[1]);
  if(argc > 2)
    fcs.write_text(argv[2]);
}
